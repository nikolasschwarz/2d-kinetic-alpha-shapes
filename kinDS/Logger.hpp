// C++ program to implement a basic logging system.
// Based on https://www.geeksforgeeks.org/cpp/logging-system-in-cpp/

#pragma once

#include <atomic>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace kinDS
{
using namespace std;
// Enum to represent log levels
enum class LogLevel : unsigned int
{
  Debug = 1,
  Info = 2,
  Warning = 4,
  Error = 8,
  Critical = 16,
  /// Targeted diagnostic tracing (flip/crossing/radius monitors). Off by default; enable with --log-add monitor.
  Monitor = 32
};

inline LogLevel operator|(LogLevel a, LogLevel b)
{
  return static_cast<LogLevel>(static_cast<unsigned>(a) | static_cast<unsigned>(b));
}

inline LogLevel operator&(LogLevel a, LogLevel b)
{
  return static_cast<LogLevel>(static_cast<unsigned>(a) & static_cast<unsigned>(b));
}

inline LogLevel& operator|=(LogLevel& a, LogLevel b)
{
  a = a | b;
  return a;
}

inline LogLevel operator&=(LogLevel& a, LogLevel b)
{
  a = a & b;
  return a;
}

inline LogLevel operator~(LogLevel a) { return static_cast<LogLevel>(~static_cast<unsigned>(a)); }

class Logger
{
 private:
  // Default: everything except Debug. Atomic so disabled levels are rejected without locking (see isEnabled).
  std::atomic<unsigned> log_level { static_cast<unsigned>(
    LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical) };
  ofstream log_file; // File stream for the log file (optional)
  bool console_enabled_ = true;
  std::function<void(LogLevel, const std::string&)> external_sink_;
  mutable std::mutex mutex_;

 public:
  // Default constructor: no log file
  Logger() = default;

  // Constructor: Opens the log file in append mode (for backward compatibility)
  Logger(const string& filename) { setLogFile(filename); }

  // Destructor: Closes the log file if open
  ~Logger()
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (log_file.is_open())
    {
      log_file.close();
    }
  }

  // Set or change the log file path. Returns true if the file is open for writing.
  bool setLogFile(const string& filename)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    // Close existing file if open
    if (log_file.is_open())
    {
      log_file.close();
    }

    if (filename.empty())
    {
      return false;
    }

    // Open new file
    log_file.open(filename, ios::app);
    if (!log_file.is_open())
    {
      cerr << "Warning: Could not open log file: " << filename << endl;
      return false;
    }
    return true;
  }

  /// When false, @ref log writes only to the log file (if open), not stdout.
  void setConsoleEnabled(bool enabled)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    console_enabled_ = enabled;
  }

  bool getConsoleEnabled() const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    return console_enabled_;
  }

  /// Optional host sink (EcoSysLab console, etc.). Invoked after the local log write.
  void setExternalSink(std::function<void(LogLevel, const std::string&)> sink)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    external_sink_ = std::move(sink);
  }

  // Explicitly set the full log level bitmask (overwrites existing mask).
  void setLogLevelMask(LogLevel mask)
  {
    log_level.store(static_cast<unsigned>(mask));
  }

  void setLogLevel(LogLevel level, bool set = true)
  {
    if (set)
    {
      log_level.fetch_or(static_cast<unsigned>(level));
    }
    else
    {
      log_level.fetch_and(~static_cast<unsigned>(level));
    }
  }

  LogLevel getLogLevel() const { return static_cast<LogLevel>(log_level.load(std::memory_order_relaxed)); }

  /// Lock-free check used by the KINDS_* macros to skip building messages for disabled levels.
  bool isEnabled(LogLevel level) const
  {
    return (log_level.load(std::memory_order_relaxed) & static_cast<unsigned>(level)) != 0;
  }

  // Logs a message with a given log level
  void log(LogLevel level, const string& message)
  {
    std::function<void(LogLevel, const std::string&)> sink;
    {
      if (!isEnabled(level))
      {
        return; // Log level not enabled
      }
      std::lock_guard<std::mutex> lock(mutex_);
      // Get current timestamp (thread-safe)
      time_t now = time(nullptr);
      tm timeinfo{};
#ifdef _WIN32
      localtime_s(&timeinfo, &now);
#else
      localtime_r(&now, &timeinfo);
#endif
      char timestamp[20];
      strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", &timeinfo);

      // Create log entry
      ostringstream logEntry;
      logEntry << "[kinDS] [" << timestamp << "] " << levelToString(level) << ": " << message << endl;

      // Output to console (flush so lines appear before a subsequent throw)
      if (console_enabled_)
      {
        if (level == LogLevel::Error || level == LogLevel::Critical)
        {
          cerr << logEntry.str();
          cerr.flush();
          fflush(stderr);
        }
        else
        {
          cout << logEntry.str();
          cout.flush();
          fflush(stdout);
        }
      }

      // Output to log file
      if (log_file.is_open())
      {
        log_file << logEntry.str();
        log_file.flush(); // Ensure immediate write to file
      }
      sink = external_sink_;
    }
    if (sink)
    {
      sink(level, message);
    }
  }

  void log(LogLevel level, const char* fmt, ...)
  {
    if (!isEnabled(level))
    {
      return; // Log level not enabled
    }

    va_list args;
    va_start(args, fmt);

    // Get the size needed
    va_list args_copy;
    va_copy(args_copy, args);
    int size = std::vsnprintf(nullptr, 0, fmt, args_copy);
    va_end(args_copy);

    if (size < 0)
    {
      va_end(args);
      log(level, std::string("formatting error in log"));
      return;
    }

    std::vector<char> buffer(static_cast<size_t>(size) + 1);
    std::vsnprintf(buffer.data(), buffer.size(), fmt, args);
    va_end(args);

    std::string log_string(buffer.data(), static_cast<size_t>(size));
    log(level, log_string);
  }

 private:
  // Converts log level to a string for output
  string levelToString(LogLevel level)
  {
    switch (level)
    {
    case LogLevel::Debug:
      return "DEBUG";
    case LogLevel::Info:
      return "INFO";
    case LogLevel::Warning:
      return "WARNING";
    case LogLevel::Error:
      return "ERROR";
    case LogLevel::Critical:
      return "CRITICAL";
    case LogLevel::Monitor:
      return "MONITOR";
    default:
      return "UNKNOWN";
    }
  }
};

// Global logger instance (no log file by default)
inline Logger logger;

}  // namespace kinDS

#define KINDS_DEBUG(msg)                                                                                               \
  if (::kinDS::logger.isEnabled(::kinDS::LogLevel::Debug))                                                            \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n";                                                         \
    ::kinDS::logger.log(::kinDS::LogLevel::Debug, ss.str());                                                           \
  }

#define KINDS_MONITOR(msg)                                                                                             \
  if (::kinDS::logger.isEnabled(::kinDS::LogLevel::Monitor))                                                          \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n";                                                         \
    ::kinDS::logger.log(::kinDS::LogLevel::Monitor, ss.str());                                                         \
  }

#define KINDS_INFO(msg)                                                                                                \
  if (::kinDS::logger.isEnabled(::kinDS::LogLevel::Info))                                                             \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n";                                                         \
    ::kinDS::logger.log(::kinDS::LogLevel::Info, ss.str());                                                            \
  }

#define KINDS_WARNING(msg)                                                                                             \
  if (::kinDS::logger.isEnabled(::kinDS::LogLevel::Warning))                                                          \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n";                                                         \
    ::kinDS::logger.log(::kinDS::LogLevel::Warning, ss.str());                                                         \
  }

#define KINDS_ERROR(msg)                                                                                               \
  if (::kinDS::logger.isEnabled(::kinDS::LogLevel::Error))                                                            \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n";                                                         \
    ::kinDS::logger.log(::kinDS::LogLevel::Error, ss.str());                                                           \
  }
