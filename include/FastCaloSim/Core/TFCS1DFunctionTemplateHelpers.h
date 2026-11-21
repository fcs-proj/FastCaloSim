// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <memory>

#include "FastCaloSim/Definitions/Index.h"
#include "TBuffer.h"

// For the purpose of FastCaloSim, 32bit are sufficient for bin counting
using TFCS1DFunction_size_t = uint32_t;

/// Converter functions to store/retrieve float numbers in the range [0,1] as
/// unsigned integers. T is checked with static_assert instead of a constrained
/// template parameter: the ROOT dictionary forward-declares this template
/// without the constraint, and the interpreter rejects the constrained
/// definition as a different redeclaration.
template<typename T, typename Tfloat>
class TFCS1DFunction_Numeric
{
  static_assert(std::unsigned_integral<T>);

public:
  static constexpr T MaxValue = std::numeric_limits<T>::max();
  static constexpr Tfloat MaxValueFloat = static_cast<Tfloat>(MaxValue);

  static auto MaxCeilOnlyForInt(const Tfloat value) -> T
  {
    return static_cast<T>(std::ceil(MaxValueFloat * value));
  }
  static auto ExpandToMaxRange(const Tfloat value) -> T
  {
    return static_cast<T>(value * (1 + MaxValueFloat));
  }
  static auto ToNormalizedRange(const T value) -> Tfloat
  {
    return static_cast<Tfloat>(value) / (1 + MaxValueFloat);
  }
};

template<typename T>
class TFCS1DFunction_Array
{
public:
  using size_t = TFCS1DFunction_size_t;

  TFCS1DFunction_Array() = default;
  TFCS1DFunction_Array(size_t count) { resize(count); }

  auto MemorySizeArray() const -> std::size_t { return size() * sizeof(T); }
  auto MemorySize() const -> std::size_t
  {
    return sizeof(*this) + MemorySizeArray();
  }

  auto size() const -> size_t { return m_size; }

  /// resize to given count, copying old content
  void resize(size_t count)
  {
    std::unique_ptr<T[]> new_cont;
    if (count > 0) {
      new_cont = std::make_unique_for_overwrite<T[]>(count);
    }
    if (m_content && new_cont) {
      std::copy_n(m_content.get(), std::min(count, size()), new_cont.get());
    }
    m_size = count;
    m_content = std::move(new_cont);
  }

  /// Direct data pointer
  auto data() -> T* { return m_content.get(); }
  auto data() const -> const T* { return m_content.get(); }

  /// Direct access operators. Values are in the range
  /// [0,TFCS1DFunction_Numeric<T>::MaxValue]
  auto operator[](size_t pos) -> T& { return m_content[pos]; }
  auto operator[](size_t pos) const -> const T& { return m_content[pos]; }

  /// begin() iterators
  auto begin() -> T* { return m_content.get(); }
  auto begin() const -> const T* { return m_content.get(); }

  /// end() iterators
  auto end() -> T* { return m_content.get() + size(); }
  auto end() const -> const T* { return m_content.get() + size(); }

private:
  std::unique_ptr<T[]> m_content;  //!
  size_t m_size {0};

  // Use ClassDef without virtual functions. Saves 8 bytes per instance
  ClassDefNV(TFCS1DFunction_Array, 1)  // TFCS1DFunction_Array
};

/// Streamer method to persitify objects of type TFCS1DFunction_Array<T>
template<typename T>
void TFCS1DFunction_Array<T>::Streamer(TBuffer& b)
{
  if (b.IsReading()) {
    size_t count;
    b >> count;
    resize(count);
    if (m_size > 0) {
      b.ReadFastArray(m_content.get(), m_size);
    }
  } else {
    b << m_size;
    if (m_size > 0) {
      b.WriteFastArray(m_content.get(), m_size);
    }
  }
}

// Class to represent histogram content. Trandom should be a type with a
// floating point range [0,1]
template<typename T, typename Trandom = float>
class TFCS1DFunction_HistogramContent
{
public:
  using size_t = TFCS1DFunction_size_t;
  using value_type = T;
  using random_type = Trandom;

  TFCS1DFunction_HistogramContent(size_t nbins = 0)
      : m_array(nbins >= 1 ? nbins - 1 : 0)
  {
  }

  auto MemorySizeArray() const -> std::size_t
  {
    return m_array.MemorySizeArray();
  }
  auto MemorySize() const -> std::size_t
  {
    return sizeof(*this) + MemorySizeArray();
  }

  /// Set the content of bin pos to a given value, where value is in the range
  /// [0,1]
  void set_fraction(size_t pos, Trandom value)
  {
    if (pos >= size()) {
      return;
    }
    m_array[pos] = TFCS1DFunction_Numeric<T, Trandom>::MaxCeilOnlyForInt(value);
  }

  /// Get the cumulative content at bin pos as fraction in the range [0,1]
  auto get_fraction(size_t pos) const -> Trandom
  {
    if (pos >= size()) {
      return 1;
    }
    return static_cast<Trandom>(m_array[pos])
        / TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat;
  }

  /// set number of bins.
  /// The actually allocated size is one smaller than count, as the last bin is
  /// fixed with the range [get_fraction(size()-1,1]
  void set_nbins(size_t nbins) { m_array.resize(nbins >= 1 ? nbins - 1 : 0); }

  /// return number of bins.
  /// This is one larger than size, as the last bin is fixed with the range
  /// [get_fraction(size()-1,1]
  auto get_nbins() const -> size_t { return size() + 1; }

  /// Get the matching bin for a given random value in the range [0,1).
  /// A residual random number to calculate a position inside this bin is
  /// returned in residual_rnd
  auto get_bin(Trandom drnd, Trandom& residual_rnd) const -> size_t
  {
    if (drnd <= 0) {
      drnd = 0;
    }
    if (drnd >= 1) {
      drnd =
          std::nextafter(static_cast<Trandom>(1.0), static_cast<Trandom>(0.0));
    }
    if (size() == 0) {
      residual_rnd = drnd;
      return 0;
    }
    T rnd = static_cast<T>(TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat
                           * drnd);
    auto it = std::upper_bound(m_array.begin(), m_array.end(), rnd);

    T basecont = 0;
    if (it != m_array.begin()) {
      basecont = *(it - 1);
    }

    T fullcont = TFCS1DFunction_Numeric<T, Trandom>::MaxValue;
    if (it != m_array.end()) {
      fullcont = *it;
    }

    T dcont = fullcont - basecont;
    if (dcont > 0) {
      residual_rnd =
          static_cast<Trandom>(rnd - basecont) / static_cast<Trandom>(dcont);
      if (residual_rnd > 1) {
        residual_rnd = std::nextafter(static_cast<Trandom>(1.0),
                                      static_cast<Trandom>(0.0));
      }
    } else {
      residual_rnd = 0.5;
    }

    if (it == m_array.end()) {
      return size();
    }
    // The bin numbers of this class are 32 bit (size_t is its own alias)
    return static_cast<size_t>(FastCaloSim::index_of(m_array, it));
  }

private:
  TFCS1DFunction_Array<T> m_array;
  auto size() const -> size_t { return m_array.size(); }

  // Use ClassDef without virtual functions. Saves 8 bytes per instance
  ClassDefNV(TFCS1DFunction_HistogramContent,
             1)  // TFCS1DFunction_HistogramContent
};

/// Fraction x in [0,1] of the way through a bin for a random number drnd in
/// [0,1], such that histogramming x gives a density with linear slope m (in
/// units of the bin width, limited to [-2,2]). Declared inline so that it is
/// inlined into position_lin in the shared library (per-hit code).
template<typename Trandom>
inline auto TFCS1DFunction_linear_slope_fraction(Trandom m, Trandom drnd)
    -> Trandom
{
  if (m > 2) {
    m = 2;
  }
  if (m < -2) {
    m = -2;
  }
  return std::abs(m) > 0.001
      ? (0.5 * std::sqrt(m * (m + 8 * drnd - 4) + 4) - 1) / m + 0.5
      : drnd;
}

// Class to represent histogram bin edges, where the interval between GetMin()
// and GetMax() can be stored as a different (smaller) data type Tint Trandom
// should be a type with a floating point range [0,1]
template<typename T, typename Tint, typename Trandom = float>
class TFCS1DFunction_HistogramCompactBinEdges
{
public:
  using size_t = TFCS1DFunction_size_t;
  using value_type = T;
  using random_type = Trandom;
  using internal_storage_type = Tint;

  TFCS1DFunction_HistogramCompactBinEdges(size_t nbins = 0)
      : m_array(nbins >= 1 ? nbins - 1 : 0)
  {
  }
  ~TFCS1DFunction_HistogramCompactBinEdges() = default;

  auto MemorySizeArray() const -> std::size_t
  {
    return m_array.MemorySizeArray();
  }
  auto MemorySize() const -> std::size_t
  {
    return sizeof(*this) + MemorySizeArray();
  }

  /// set number of bins
  void set_nbins(size_t nbins) { m_array.resize(nbins >= 1 ? nbins - 1 : 0); }

  /// return number of bins
  auto get_nbins() const -> size_t { return size() + 1; }

  /// set position of lower edge of bins. Requires GetMin() and GetMax() to be
  /// set and may not be changed!
  void SetBinLowEdge(size_t pos, const T& value)
  {
    if (pos == 0) {
      SetMin(value);
      return;
    }
    if (pos >= get_nbins()) {
      SetMax(value);
      return;
    }
    m_array[pos - 1] = TFCS1DFunction_Numeric<Tint, T>::ExpandToMaxRange(
        (value - GetMin()) / Length());
  }

  /// get position of lower edge of bins. Requires GetMin() and GetMax() to be
  /// set and may not be changed!
  auto GetBinLowEdge(size_t pos) const -> T
  {
    if (pos == 0) {
      return GetMin();
    }
    if (pos >= get_nbins()) {
      return GetMax();
    }
    return GetMin()
        + Length()
        * TFCS1DFunction_Numeric<Tint, T>::ToNormalizedRange(m_array[pos - 1]);
  }

  /// get the length of a bin
  auto GetBinLength(size_t pos) const -> T
  {
    return GetBinLowEdge(pos + 1) - GetBinLowEdge(pos);
  }

  /// set and get minimum
  void SetMin(const T& value) { m_Min = value; }
  auto GetMin() const -> const T& { return m_Min; }

  /// set and get maximum
  void SetMax(const T& value) { m_Max = value; }
  auto GetMax() const -> const T& { return m_Max; }

  /// set minimum and maximum
  void SetMinMax(const T& valuemin, const T& valuemax)
  {
    SetMin(valuemin);
    SetMax(valuemax);
  }

  /// Get length of interval of all bins
  auto Length() const -> T { return GetMax() - GetMin(); }

  /// return linear interpolated position for bin pos. Interpolation is done
  /// with a random value in the range [0,1]
  auto position(size_t pos, const Trandom& drnd) const -> T
  {
    T pos1 = GetBinLowEdge(pos);
    T pos2 = GetBinLowEdge(pos + 1);
    return (1 - drnd) * pos1 + drnd * pos2;
  }

  /// return linearly interpolated position for bin pos, such that histograming
  /// the position gives a linear slope m, where m is in units of the bin width
  /// for bin pos. Interpolation is done with a random value in the range [0,1]
  auto position_lin(size_t pos, Trandom m, const Trandom& drnd) const -> T
  {
    const Trandom x = TFCS1DFunction_linear_slope_fraction(m, drnd);
    T pos1 = GetBinLowEdge(pos);
    T pos2 = GetBinLowEdge(pos + 1);
    return (1 - x) * pos1 + x * pos2;
  }

private:
  TFCS1DFunction_Array<Tint> m_array;
  auto size() const -> size_t { return m_array.size(); }
  T m_Min {0};
  T m_Max {0};

  ClassDefNV(TFCS1DFunction_HistogramCompactBinEdges,
             1)  // TFCS1DFunction_HistogramCompactBinEdges
};

class TFCS1DFunction_HistogramInt8BinEdges
    : public TFCS1DFunction_HistogramCompactBinEdges<float, uint8_t, float>
{
public:
  using TFCS1DFunction_HistogramCompactBinEdges::
      TFCS1DFunction_HistogramCompactBinEdges;

  ClassDefNV(TFCS1DFunction_HistogramInt8BinEdges,
             1)  // TFCS1DFunction_HistogramInt8BinEdges
};

class TFCS1DFunction_HistogramInt32BinEdges
    : public TFCS1DFunction_HistogramCompactBinEdges<float, uint32_t, float>
{
public:
  using TFCS1DFunction_HistogramCompactBinEdges::
      TFCS1DFunction_HistogramCompactBinEdges;

  ClassDefNV(TFCS1DFunction_HistogramInt32BinEdges,
             1)  // TFCS1DFunction_HistogramInt32BinEdges
};
