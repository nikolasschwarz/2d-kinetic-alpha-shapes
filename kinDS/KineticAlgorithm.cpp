#include "KineticAlgorithm.hpp"

using namespace kinDS;

void KineticAlgorithm::clear()
{
  while (!events_.empty())
  {
    events_.pop();
  }
  next_queue_sequence_ = 0;
}

void KineticAlgorithm::processEvents()
{
  while (!events_.empty())
  {
    auto event = events_.top();
    events_.pop();
    event->handleEvent();
  }
}
