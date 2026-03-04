#include "kinDS/Logger.hpp"

using namespace kinDS;

namespace
{
struct EnableDebugLogLevel
{
  EnableDebugLogLevel()
  {
    logger.setLogLevelMask(
      LogLevel::Debug | LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical);
  }
};

static EnableDebugLogLevel enableDebugLogLevel;
} // namespace

