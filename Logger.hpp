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
using namespace std;

// Enum to represent log levels
enum LogLevel
{
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    CRITICAL
};

class Logger
{
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

    // Logs a message with a given log level
    void log(LogLevel level, const string& message)
    {
        // Get current timestamp
        time_t now = time(0);
        tm* timeinfo = localtime(&now);
        char timestamp[20];
        strftime(timestamp, sizeof(timestamp),
            "%Y-%m-%d %H:%M:%S", timeinfo);

        // Create log entry
        ostringstream logEntry;
        logEntry << "[" << timestamp << "] "
                 << levelToString(level) << ": " << message
                 << endl;

        // Output to console
        cout << logEntry.str();

        // Output to log file
        if (logFile.is_open())
        {
            logFile << logEntry.str();
            logFile
                .flush(); // Ensure immediate write to file
        }
    }

    void log(LogLevel level, const char* fmt, ...)
    {
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

static Logger logger("logfile.txt");
