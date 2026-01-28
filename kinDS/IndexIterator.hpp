#pragma once
#include <iterator>

struct IndexIterator
{
  using value_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = size_t;
  using pointer = void;
  using iterator_category = std::forward_iterator_tag;

  size_t i;

  size_t operator*() const { return i; }

  IndexIterator& operator++()
  {
    ++i;
    return *this;
  }
  IndexIterator operator+(difference_type n) const { return { i + n }; }
  difference_type operator-(const IndexIterator& other) const { return i - other.i; }

  bool operator==(const IndexIterator& other) const { return i == other.i; }
  bool operator!=(const IndexIterator& other) const { return i != other.i; }
};
