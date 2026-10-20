#pragma once

#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace kinDS
{

/// Chunked array. A new chunk is published under a mutex; indexing a slot below an observed
/// @c size() takes no lock. Elements do not move when the table grows. Same-slot writes are the caller's problem.
template<typename T>
class SlotTable
{
 public:
  static constexpr size_t chunk_shift = 12;
  static constexpr size_t chunk_size = size_t { 1 } << chunk_shift;
  static constexpr size_t max_chunks = 4096;

  SlotTable() = default;
  SlotTable(const SlotTable&) = delete;
  SlotTable& operator=(const SlotTable&) = delete;
  SlotTable(SlotTable&& other) noexcept { moveFrom(other); }
  SlotTable& operator=(SlotTable&& other) noexcept
  {
    if (this != &other)
    {
      clearUnlocked();
      moveFrom(other);
    }
    return *this;
  }
  ~SlotTable() { clearUnlocked(); }

  size_t size() const { return size_.load(std::memory_order_acquire); }
  bool empty() const { return size() == 0; }

  T& operator[](size_t index) { return element(index); }
  const T& operator[](size_t index) const { return element(index); }

  T& at(size_t index)
  {
    if (index >= size())
    {
      throw std::out_of_range("SlotTable::at");
    }
    return element(index);
  }
  const T& at(size_t index) const
  {
    if (index >= size())
    {
      throw std::out_of_range("SlotTable::at");
    }
    return element(index);
  }

  /// Allocates chunks for @p n slots. Does not change @ref size.
  void reserve(size_t n)
  {
    std::lock_guard<std::recursive_mutex> lock(grow_mutex_);
    ensureChunksUnlocked(n);
  }

  void clear()
  {
    std::lock_guard<std::recursive_mutex> lock(grow_mutex_);
    clearUnlocked();
  }

  void resize(size_t n) { resizeWith(n, [] { return T {}; }); }

  void resize(size_t n, const T& fill)
  {
    resizeWith(n, [&] { return fill; });
  }

  /// Grows or shrinks to @p n. New elements are built by @p factory before the new size is published.
  template<class Factory>
  void resizeWith(size_t n, Factory&& factory)
  {
    std::lock_guard<std::recursive_mutex> lock(grow_mutex_);
    const size_t old = size_.load(std::memory_order_relaxed);
    if (n < old)
    {
      for (size_t i = n; i < old; ++i)
      {
        elementUnlocked(i) = T {};
      }
      size_.store(n, std::memory_order_release);
      return;
    }
    if (n == old)
    {
      return;
    }
    ensureChunksUnlocked(n);
    for (size_t i = old; i < n; ++i)
    {
      elementUnlocked(i) = factory();
    }
    size_.store(n, std::memory_order_release);
  }

  void push_back(const T& value) { append(value); }
  void push_back(T&& value) { append(std::move(value)); }
  void emplace_back() { append(T {}); }

  /// Drops existing chunks and takes ownership of @p values. Single-threaded rebuilds use this.
  void replaceWith(std::vector<T> values)
  {
    std::lock_guard<std::recursive_mutex> lock(grow_mutex_);
    clearUnlocked();
    if (values.empty())
    {
      return;
    }
    ensureChunksUnlocked(values.size());
    for (size_t i = 0; i < values.size(); ++i)
    {
      elementUnlocked(i) = std::move(values[i]);
    }
    size_.store(values.size(), std::memory_order_release);
  }

  template<bool IsConst>
  class BasicIterator
  {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using Table = std::conditional_t<IsConst, const SlotTable, SlotTable>;
    using reference = std::conditional_t<IsConst, const T&, T&>;
    using pointer = std::conditional_t<IsConst, const T*, T*>;

    BasicIterator() = default;
    BasicIterator(Table* table, size_t index)
      : table_(table)
      , index_(index)
    {
    }

    reference operator*() const { return (*table_)[index_]; }
    pointer operator->() const { return &(*table_)[index_]; }
    BasicIterator& operator++()
    {
      ++index_;
      return *this;
    }
    BasicIterator operator++(int)
    {
      BasicIterator copy = *this;
      ++*this;
      return copy;
    }
    bool operator==(const BasicIterator& other) const { return index_ == other.index_; }
    bool operator!=(const BasicIterator& other) const { return index_ != other.index_; }

   private:
    Table* table_ = nullptr;
    size_t index_ = 0;
  };

  using iterator = BasicIterator<false>;
  using const_iterator = BasicIterator<true>;

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size()); }
  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, size()); }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }

 private:
  struct Chunk
  {
    T elements[chunk_size];
  };

  T& element(size_t index)
  {
    const size_t chunk_index = index >> chunk_shift;
    assert(chunk_index < max_chunks);
    Chunk* chunk = chunks_[chunk_index].load(std::memory_order_acquire);
    assert(chunk != nullptr);
    return chunk->elements[index & (chunk_size - 1)];
  }
  const T& element(size_t index) const
  {
    const size_t chunk_index = index >> chunk_shift;
    assert(chunk_index < max_chunks);
    Chunk* chunk = chunks_[chunk_index].load(std::memory_order_acquire);
    assert(chunk != nullptr);
    return chunk->elements[index & (chunk_size - 1)];
  }

  T& elementUnlocked(size_t index)
  {
    Chunk* chunk = chunks_[index >> chunk_shift].load(std::memory_order_relaxed);
    assert(chunk != nullptr);
    return chunk->elements[index & (chunk_size - 1)];
  }

  void ensureChunksUnlocked(size_t n)
  {
    if (n == 0)
    {
      return;
    }
    const size_t needed = (n + chunk_size - 1) >> chunk_shift;
    if (needed > max_chunks)
    {
      throw std::length_error("SlotTable capacity exceeded");
    }
    size_t published = published_.load(std::memory_order_relaxed);
    while (published < needed)
    {
      chunks_[published].store(new Chunk(), std::memory_order_release);
      ++published;
    }
    published_.store(published, std::memory_order_release);
  }

  template<class U>
  void append(U&& value)
  {
    std::lock_guard<std::recursive_mutex> lock(grow_mutex_);
    const size_t index = size_.load(std::memory_order_relaxed);
    ensureChunksUnlocked(index + 1);
    elementUnlocked(index) = std::forward<U>(value);
    size_.store(index + 1, std::memory_order_release);
  }

  void clearUnlocked()
  {
    size_.store(0, std::memory_order_relaxed);
    const size_t published = published_.exchange(0, std::memory_order_relaxed);
    for (size_t i = 0; i < published; ++i)
    {
      delete chunks_[i].exchange(nullptr, std::memory_order_relaxed);
    }
  }

  void moveFrom(SlotTable& other) noexcept
  {
    const size_t published = other.published_.exchange(0, std::memory_order_relaxed);
    size_.store(other.size_.exchange(0, std::memory_order_relaxed), std::memory_order_relaxed);
    published_.store(published, std::memory_order_relaxed);
    for (size_t i = 0; i < published; ++i)
    {
      chunks_[i].store(other.chunks_[i].exchange(nullptr, std::memory_order_relaxed), std::memory_order_relaxed);
    }
  }

  std::recursive_mutex grow_mutex_;
  std::atomic<size_t> size_ { 0 };
  std::atomic<size_t> published_ { 0 };
  std::array<std::atomic<Chunk*>, max_chunks> chunks_ {};
};

} // namespace kinDS
