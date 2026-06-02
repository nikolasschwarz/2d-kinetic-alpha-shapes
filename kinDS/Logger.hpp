// C++ program to implement a basic logging system.
// Based on https://www.geeksforgeeks.org/cpp/logging-system-in-cpp/

#pragma once

#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
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
  Critical = 16
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
  // Default: everything except Debug
  LogLevel log_level = LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical;
  ofstream log_file; // File stream for the log file (optional)

 public:
  // Default constructor: no log file
  Logger() = default;

  // Constructor: Opens the log file in append mode (for backward compatibility)
  Logger(const string& filename) { setLogFile(filename); }

  // Destructor: Closes the log file if open
  ~Logger()
  {
    if (log_file.is_open())
    {
      log_file.close();
    }
  }

  // Set or change the log file path
  void setLogFile(const string& filename)
  {
    // Close existing file if open
    if (log_file.is_open())
    {
      log_file.close();
    }

    // Open new file
    log_file.open(filename, ios::app);
    if (!log_file.is_open())
    {
      cerr << "Warning: Could not open log file: " << filename << endl;
    }
  }

  // Explicitly set the full log level bitmask (overwrites existing mask).
  void setLogLevelMask(LogLevel mask) { log_level = mask; }

  void setLogLevel(LogLevel log_level, bool set = true)
  {
    if (set)
    {
      this->log_level |= log_level;
    }
    else
    {
      log_level = ~log_level;
      this->log_level &= log_level;
    }
  }

  LogLevel getLogLevel() const { return log_level; }

  // Logs a message with a given log level
  void log(LogLevel level, const string& message)
  {
    if (!static_cast<unsigned>(log_level & level))
    {
      return; // Log level not enabled
    }
    // Get current timestamp
    time_t now = time(0);
    tm* timeinfo = localtime(&now);
    char timestamp[20];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);

    // Create log entry
    ostringstream logEntry;
    logEntry << "[kinDS] [" << timestamp << "] " << levelToString(level) << ": " << message << endl;

    // Output to console
    cout << logEntry.str();

    // Output to log file
    if (log_file.is_open())
    {
      log_file << logEntry.str();
      log_file.flush(); // Ensure immediate write to file
    }
  }

  void log(LogLevel level, const char* fmt, ...)
  {

    if (!static_cast<unsigned>(log_level & level))
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
    }

    std::vector<char> buffer(size + 1);
    std::vsnprintf(buffer.data(), buffer.size(), fmt, args);
    va_end(args);

    std::string log_string(buffer.data(), size);
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
    default:
      return "UNKNOWN";
    }
  }
};

// Global logger instance (no log file by default)
inline Logger logger;

#define KINDS_DEBUG(msg)                                                                                               \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::Debug, ss.str());                                                                             \
  }

#define KINDS_INFO(msg)                                                                                                 \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::Info, ss.str());                                                                              \
  }

#define KINDS_WARNING(msg)                                                                                             \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::Warning, ss.str());                                                                           \
  }

#define KINDS_ERROR(msg)                                                                                               \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::Error, ss.str());                                                                             \
  }
}
