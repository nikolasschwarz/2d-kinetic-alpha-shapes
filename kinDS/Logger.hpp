// C++ program to implement a basic logging system.
// Source: https://www.geeksforgeeks.org/cpp/logging-system-in-cpp/

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
enum LogLevel : unsigned int
{
  DEBUG = 1,
  INFO = 2,
  WARNING = 4,
  ERROR = 8,
  CRITICAL = 16
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
  LogLevel log_level = DEBUG | INFO | WARNING | ERROR | CRITICAL;

 public:
  // Constructor: Opens the log file in append mode
  Logger(const string& filename)
  {
    logFile.open(filename, ios::app);
    if (!logFile.is_open())
    {
      cerr << "Error opening log file." << endl;
    }
  }

  // Destructor: Closes the log file
  ~Logger() { logFile.close(); }

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
    if (!(log_level & level))
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
    logEntry << "[" << timestamp << "] " << levelToString(level) << ": " << message << endl;

    // Output to console
    cout << logEntry.str();

    // Output to log file
    if (logFile.is_open())
    {
      logFile << logEntry.str();
      logFile.flush(); // Ensure immediate write to file
    }
  }

  void log(LogLevel level, const char* fmt, ...)
  {

    if (!(log_level & level))
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
  ofstream logFile; // File stream for the log file

  // Converts log level to a string for output
  string levelToString(LogLevel level)
  {
    switch (level)
    {
    case DEBUG:
      return "DEBUG";
    case INFO:
      return "INFO";
    case WARNING:
      return "WARNING";
    case ERROR:
      return "ERROR";
    case CRITICAL:
      return "CRITICAL";
    default:
      return "UNKNOWN";
    }
  }
};

static Logger logger("kinDS_logfile.txt");

#define KINDS_INFO(msg)                                                                                                \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::INFO, ss.str());                                                                              \
  }

#define KINDS_ERROR(msg)                                                                                               \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::ERROR, ss.str());                                                                             \
  }

#define KINDS_WARNING(msg)                                                                                             \
  {                                                                                                                    \
    std::stringstream ss;                                                                                              \
    ss << msg << " (" << __FILE__ << ": line " << __LINE__ << ")\n";                                                   \
    logger.log(LogLevel::WARNING, ss.str());                                                                           \
  }
}
