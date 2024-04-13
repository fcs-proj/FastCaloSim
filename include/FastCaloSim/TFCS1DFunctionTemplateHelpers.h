#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCS1DFunctionTemplateHelpers_h
#  define ISF_FASTCALOSIMEVENT_TFCS1DFunctionTemplateHelpers_h

#  include <algorithm>
#  include <cmath>
#  include <cstring>
#  include <vector>

#  include "TBuffer.h"

// For the purpose of FastCaloSim, 32bit are sufficient for bin counting
typedef uint32_t TFCS1DFunction_size_t;

/// Converter functions that does nothing for floats
template<typename T, typename Tfloat = float>
class TFCS1DFunction_Numeric
{
public:
  static constexpr T MaxValue = 1;
  static constexpr Tfloat MaxValueFloat = MaxValue;

  static inline T MaxCeilOnlyForInt(const Tfloat value) { return value; };
  static inline T ExpandToMaxRange(const Tfloat value) { return value; };
  static inline Tfloat ToNormalizedRange(const T value) { return value; };
};

/// Converter functions to store/retrieve float numbers in integers
template<typename Tfloat>
class TFCS1DFunction_Numeric<uint8_t, Tfloat>
{
public:
  static constexpr uint8_t MaxValue = UINT8_MAX;
  static constexpr Tfloat MaxValueFloat = MaxValue;

  static inline uint8_t MaxCeilOnlyForInt(const Tfloat value)
  {
    return std::ceil(MaxValueFloat * value);
  };
  static inline uint8_t ExpandToMaxRange(const Tfloat value)
  {
    return value * (1 + MaxValueFloat);
  };
  static inline Tfloat ToNormalizedRange(const uint8_t value)
  {
    return value / (1 + MaxValueFloat);
  };
};

/// Converter functions to store/retrieve float numbers in integers
template<typename Tfloat>
class TFCS1DFunction_Numeric<uint16_t, Tfloat>
{
public:
  static constexpr uint16_t MaxValue = UINT16_MAX;
  static constexpr Tfloat MaxValueFloat = MaxValue;

  static inline uint16_t MaxCeilOnlyForInt(const Tfloat value)
  {
    return std::ceil(MaxValueFloat * value);
  };
  static inline uint16_t ExpandToMaxRange(const Tfloat value)
  {
    return value * (1 + MaxValueFloat);
  };
  static inline Tfloat ToNormalizedRange(const uint16_t value)
  {
    return value / (1 + MaxValueFloat);
  };
};

/// Converter functions to store/retrieve float numbers in integers
template<typename Tfloat>
class TFCS1DFunction_Numeric<uint32_t, Tfloat>
{
public:
  static constexpr uint32_t MaxValue = UINT32_MAX;
  static constexpr Tfloat MaxValueFloat = MaxValue;

  static inline uint32_t MaxCeilOnlyForInt(const Tfloat value)
  {
    return std::ceil(MaxValueFloat * value);
  };
  static inline uint32_t ExpandToMaxRange(const Tfloat value)
  {
    return value * (1 + MaxValueFloat);
  };
  static inline Tfloat ToNormalizedRange(const uint32_t value)
  {
    return value / (1 + MaxValueFloat);
  };
};

template<typename T>
class TFCS1DFunction_Array
{
public:
  typedef TFCS1DFunction_size_t size_t;

  TFCS1DFunction_Array() {};
  TFCS1DFunction_Array(size_t count) { resize(count); };
  ~TFCS1DFunction_Array()
  {
    if (m_content)
      delete[] m_content;
  };

  std::size_t MemorySizeArray() const { return size() * sizeof(T); };
  std::size_t MemorySize() const { return sizeof(*this) + MemorySizeArray(); };

  inline size_t size() const { return m_size; };

  /// resize to given count, copying old content
  void resize(size_t count)
  {
    T* new_cont = nullptr;
    if (count > 0)
      new_cont = new T[count];
    if (m_content && new_cont) {
      size_t ncopy = count;
      if (size() < ncopy)
        ncopy = size();
      ncopy *= sizeof(T);
      std::memcpy(new_cont, m_content, ncopy);
    }
    if (m_content)
      delete[] m_content;
    m_size = count;
    m_content = new_cont;
  };

  /// Direct data pointer
  inline T* data() { return m_content; };
  inline const T* data() const { return m_content; };

  /// Direct access operators. Values are in the range
  /// [0,TFCS1DFunction_Numeric<T>::MaxValue]
  inline T& operator[](size_t pos) { return m_content[pos]; };
  inline const T& operator[](size_t pos) const { return m_content[pos]; };

  /// begin() iterators
  inline T* begin() { return m_content; };
  inline const T* begin() const { return m_content; };

  /// end() iterators
  inline T* end() { return m_content + size(); };
  inline const T* end() const { return m_content + size(); };

private:
  T* m_content {nullptr};  //!
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
    if (m_size > 0)
      b.ReadFastArray(m_content, m_size);
  } else {
    b << m_size;
    if (m_size > 0)
      b.WriteFastArray(m_content, m_size);
  }
}

// Class to represent histogram content. Trandom should be a type with a
// floating point range [0,1]
template<typename T, typename Trandom = float>
class TFCS1DFunction_HistogramContent
{
public:
  typedef TFCS1DFunction_size_t size_t;
  typedef T value_type;
  typedef Trandom random_type;

  TFCS1DFunction_HistogramContent(size_t nbins = 0)
      : m_array(nbins >= 1 ? nbins - 1 : 0) {};

  std::size_t MemorySizeArray() const { return m_array.MemorySizeArray(); };
  std::size_t MemorySize() const { return sizeof(*this) + MemorySizeArray(); };

  /// Set the content of bin pos to a given value, where value is in the range
  /// [0,1]
  void set_fraction(size_t pos, Trandom value)
  {
    if (pos >= size())
      return;
    m_array[pos] = TFCS1DFunction_Numeric<T, Trandom>::MaxCeilOnlyForInt(value);
  };

  /// Get the cumulative content at bin pos as fraction in the range [0,1]
  inline Trandom get_fraction(size_t pos) const
  {
    if (pos >= size())
      return 1;
    return m_array[pos] / TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat;
  };

  /// Get the content at bin pos as fraction in the range [0,1]
  inline Trandom get_binfraction(size_t pos) const
  {
    if (pos == 0)
      return m_array[pos] / TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat;
    if (pos == size())
      return 1
          - (m_array[size() - 1]
             / TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat);
    return (m_array[pos] - m_array[pos - 1])
        / TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat;
  };

  /// set number of bins.
  /// The actually allocated size is one smaller than count, as the last bin is
  /// fixed with the range [get_fraction(size()-1,1]
  void set_nbins(size_t nbins) { m_array.resize(nbins >= 1 ? nbins - 1 : 0); };

  /// return number of bins.
  /// This is one larger than size, as the last bin is fixed with the range
  /// [get_fraction(size()-1,1]
  inline size_t get_nbins() const { return size() + 1; };

  /// Get the matching bin for a given random value in the range [0,1).
  /// A residual random number to calculate a position inside this bin is
  /// returned in residual_rnd
  size_t get_bin(Trandom drnd, Trandom& residual_rnd) const
  {
    if (drnd <= 0)
      drnd = 0;
    if (drnd >= 1)
      drnd = std::nextafter((Trandom)1.0, (Trandom)0.0);
    if (size() == 0) {
      residual_rnd = drnd;
      return 0;
    }
    T rnd = TFCS1DFunction_Numeric<T, Trandom>::MaxValueFloat * drnd;
    auto it = std::upper_bound(m_array.begin(), m_array.end(), rnd);

    T basecont = 0;
    if (it != m_array.begin())
      basecont = *(it - 1);

    T fullcont = TFCS1DFunction_Numeric<T, Trandom>::MaxValue;
    if (it != m_array.end())
      fullcont = *it;

    T dcont = fullcont - basecont;
    if (dcont > 0) {
      residual_rnd = ((Trandom)(rnd - basecont)) / dcont;
      if (residual_rnd > 1)
        residual_rnd = std::nextafter((Trandom)1.0, (Trandom)0.0);
    } else {
      residual_rnd = 0.5;
    }

    if (it == m_array.end())
      return size();
    return std::distance(m_array.begin(), it);
  };

private:
  TFCS1DFunction_Array<T> m_array;
  inline size_t size() const { return m_array.size(); };

  // Use ClassDef without virtual functions. Saves 8 bytes per instance
  ClassDefNV(TFCS1DFunction_HistogramContent,
             1)  // TFCS1DFunction_HistogramContent
};

// Class to represent histogram bin edges. Trandom should be a type with a
// floating point range [0,1]
template<typename T, typename Trandom = float>
class TFCS1DFunction_HistogramBinEdges
{
public:
  typedef TFCS1DFunction_size_t size_t;
  typedef T value_type;
  typedef Trandom random_type;

  TFCS1DFunction_HistogramBinEdges(size_t nbins = 0)
      : m_array(nbins + 1) {};
  ~TFCS1DFunction_HistogramBinEdges() {};

  std::size_t MemorySizeArray() const { return m_array.MemorySizeArray(); };
  std::size_t MemorySize() const { return sizeof(*this) + MemorySizeArray(); };

  /// set number of bins
  void set_nbins(size_t nbins) { m_array.resize(nbins + 1); };

  /// return number of bins
  inline size_t get_nbins() const { return size() - 1; };

  /// set position of lower edge of bins
  void SetBinLowEdge(size_t pos, const T& value) { m_array[pos] = value; };

  /// get position of lower edge of bins
  inline const T& GetBinLowEdge(size_t pos) const { return m_array[pos]; };

  /// get the length of a bin
  inline const T GetBinLength(size_t pos) const
  {
    return GetBinLowEdge(pos + 1) - GetBinLowEdge(pos);
  };

  /// set and get minimum
  void SetMin(const T& value) { m_array[0] = value; };
  inline const T& GetMin() const { return m_array[0]; };

  /// set and get maximum
  void SetMax(const T& value) { m_array[get_nbins()] = value; };
  inline const T& GetMax() const { return m_array[get_nbins()]; };

  /// set minimum and maximum
  void SetMinMax(const T& valuemin, const T& valuemax)
  {
    SetMin(valuemin);
    SetMax(valuemax);
  };

  /// Get length of interval of all bins
  inline const T Length() const { return GetMax() - GetMin(); };

  /// return linear interpolated position for bin pos. Interpolation is done
  /// with a random value in the range [0,1]
  inline T position(size_t pos, const Trandom& drnd) const
  {
    return (1 - drnd) * m_array[pos] + drnd * m_array[pos + 1];
  };

  /// return linearly interpolated position for bin pos, such that histograming
  /// the position gives a linear slope m, where m is in units of the bin width
  /// for bin pos. Interpolation is done with a random value in the range [0,1]
  inline T position_lin(size_t pos, Trandom m, const Trandom& drnd) const
  {
    if (m > 2)
      m = 2;
    if (m < -2)
      m = -2;
    Trandom x = fabs(m) > 0.001
        ? (0.5 * std::sqrt(m * (m + 8 * drnd - 4) + 4) - 1) / m + 0.5
        : drnd;
    return (1 - x) * m_array[pos] + x * m_array[pos + 1];
  };

  /// return exponentially interpolated position for bin pos, such that
  /// histograming the position gives a linear slope m, where m is in units of
  /// the bin width for bin pos. Interpolation is done with a random value in
  /// the range [0,1]
  inline T position_exp(size_t pos, Trandom beta, const Trandom& drnd) const
  {
    Trandom z = drnd;  ///(m_array[pos+1] - m_array[pos]);
    T pos1 = GetBinLowEdge(pos);
    T pos2 = GetBinLowEdge(pos + 1);
    if (fabs(beta) < 1.0e-8)
      return (1 - z) * pos1 + z * pos2;
    else {
      T deltax = m_array[pos + 1] - m_array[pos];
      if (deltax == 0)
        return m_array[pos];
      else {
        z = 1 / beta * log(1.0 + drnd * (exp(beta * deltax) - 1.0)) / deltax;
      }
    }

    return (1 - z) * pos1 + z * pos2;
  };

private:
  TFCS1DFunction_Array<T> m_array;
  inline size_t size() const { return m_array.size(); };

  // Use ClassDef without virtual functions. Saves 8 bytes per instance
  ClassDefNV(TFCS1DFunction_HistogramBinEdges,
             1)  // TFCS1DFunction_HistogramBinEdges
};

class TFCS1DFunction_HistogramFloatBinEdges
    : public TFCS1DFunction_HistogramBinEdges<float, float>
{
public:
  TFCS1DFunction_HistogramFloatBinEdges(size_t nbins = 0)
      : TFCS1DFunction_HistogramBinEdges<float, float>(nbins) {};

  ClassDefNV(TFCS1DFunction_HistogramFloatBinEdges,
             1)  // TFCS1DFunction_HistogramFloatBinEdges
};

class TFCS1DFunction_HistogramDoubleBinEdges
    : public TFCS1DFunction_HistogramBinEdges<double, double>
{
public:
  TFCS1DFunction_HistogramDoubleBinEdges(size_t nbins = 0)
      : TFCS1DFunction_HistogramBinEdges<double, double>(nbins) {};

  ClassDefNV(TFCS1DFunction_HistogramDoubleBinEdges,
             1)  // TFCS1DFunction_HistogramDoubleBinEdges
};

// Class to represent histogram bin edges, where the interval between GetMin()
// and GetMax() can be stored as a different (smaller) data type Tint Trandom
// should be a type with a floating point range [0,1]
template<typename T, typename Tint, typename Trandom = float>
class TFCS1DFunction_HistogramCompactBinEdges
{
public:
  typedef TFCS1DFunction_size_t size_t;
  typedef T value_type;
  typedef Trandom random_type;
  typedef Tint internal_storage_type;

  TFCS1DFunction_HistogramCompactBinEdges(size_t nbins = 0)
      : m_array(nbins >= 1 ? nbins - 1 : 0) {};
  ~TFCS1DFunction_HistogramCompactBinEdges() {};

  std::size_t MemorySizeArray() const { return m_array.MemorySizeArray(); };
  std::size_t MemorySize() const { return sizeof(*this) + MemorySizeArray(); };

  /// set number of bins
  void set_nbins(size_t nbins) { m_array.resize(nbins >= 1 ? nbins - 1 : 0); };

  /// return number of bins
  inline size_t get_nbins() const { return size() + 1; };

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
  };

  /// get position of lower edge of bins. Requires GetMin() and GetMax() to be
  /// set and may not be changed!
  inline const T GetBinLowEdge(size_t pos) const
  {
    if (pos == 0)
      return GetMin();
    if (pos >= get_nbins())
      return GetMax();
    return GetMin()
        + Length()
        * TFCS1DFunction_Numeric<Tint, T>::ToNormalizedRange(m_array[pos - 1]);
  };

  /// get the length of a bin
  inline const T GetBinLength(size_t pos) const
  {
    return GetBinLowEdge(pos + 1) - GetBinLowEdge(pos);
  };

  /// set and get minimum
  void SetMin(const T& value) { m_Min = value; };
  inline const T& GetMin() const { return m_Min; };

  /// set and get maximum
  void SetMax(const T& value) { m_Max = value; };
  inline const T& GetMax() const { return m_Max; };

  /// set minimum and maximum
  void SetMinMax(const T& valuemin, const T& valuemax)
  {
    SetMin(valuemin);
    SetMax(valuemax);
  };

  /// Get length of interval of all bins
  inline const T Length() const { return GetMax() - GetMin(); };

  /// return linear interpolated position for bin pos. Interpolation is done
  /// with a random value in the range [0,1]
  inline T position(size_t pos, const Trandom& drnd) const
  {
    T pos1 = GetBinLowEdge(pos);
    T pos2 = GetBinLowEdge(pos + 1);
    return (1 - drnd) * pos1 + drnd * pos2;
  };

  /// return interpolated position for bin pos, such that histograming the
  /// position gives a linear slope m, where m is in units of the bin width for
  /// bin pos. Interpolation is done with a random value in the range [0,1]

  /// return linearly interpolated position for bin pos, such that histograming
  /// the position gives a linear slope m, where m is in units of the bin width
  /// for bin pos. Interpolation is done with a random value in the range [0,1]
  inline T position_lin(size_t pos, Trandom m, const Trandom& drnd) const
  {
    if (m > 2)
      m = 2;
    if (m < -2)
      m = -2;
    Trandom x = fabs(m) > 0.001
        ? (0.5 * std::sqrt(m * (m + 8 * drnd - 4) + 4) - 1) / m + 0.5
        : drnd;
    T pos1 = GetBinLowEdge(pos);
    T pos2 = GetBinLowEdge(pos + 1);
    return (1 - x) * pos1 + x * pos2;
  };

  /// return exponentially interpolated position for bin pos, such that
  /// histograming the position gives a linear slope m, where m is in units of
  /// the bin width for bin pos. Interpolation is done with a random value in
  /// the range [0,1]
  inline T position_exp(size_t pos, Trandom beta, const Trandom& drnd) const
  {
    Trandom z = drnd;  ///(m_array[pos+1] - m_array[pos]);
    T pos1 = GetBinLowEdge(pos);
    T pos2 = GetBinLowEdge(pos + 1);
    if (fabs(beta) < 1.0e-8)
      return (1 - z) * pos1 + z * pos2;
    else {
      T deltax = m_array[pos + 1] - m_array[pos];
      if (deltax == 0)
        return m_array[pos];
      else {
        z = 1 / beta * log(1.0 + drnd * (exp(beta * deltax) - 1.0)) / deltax;
      }
    }

    return (1 - z) * pos1 + z * pos2;
  };

private:
  TFCS1DFunction_Array<Tint> m_array;
  inline size_t size() const { return m_array.size(); };
  T m_Min {0};
  T m_Max {0};

  ClassDefNV(TFCS1DFunction_HistogramCompactBinEdges,
             1)  // TFCS1DFunction_HistogramCompactBinEdges
};

class TFCS1DFunction_HistogramInt8BinEdges
    : public TFCS1DFunction_HistogramCompactBinEdges<float, uint8_t, float>
{
public:
  TFCS1DFunction_HistogramInt8BinEdges(size_t nbins = 0)
      : TFCS1DFunction_HistogramCompactBinEdges<float, uint8_t, float>(
          nbins) {};

  ClassDefNV(TFCS1DFunction_HistogramInt8BinEdges,
             1)  // TFCS1DFunction_HistogramInt8BinEdges
};

class TFCS1DFunction_HistogramInt16BinEdges
    : public TFCS1DFunction_HistogramCompactBinEdges<float, uint16_t, float>
{
public:
  TFCS1DFunction_HistogramInt16BinEdges(size_t nbins = 0)
      : TFCS1DFunction_HistogramCompactBinEdges<float, uint16_t, float>(
          nbins) {};

  ClassDefNV(TFCS1DFunction_HistogramInt16BinEdges,
             1)  // TFCS1DFunction_HistogramInt16BinEdges
};

class TFCS1DFunction_HistogramInt32BinEdges
    : public TFCS1DFunction_HistogramCompactBinEdges<float, uint32_t, float>
{
public:
  TFCS1DFunction_HistogramInt32BinEdges(size_t nbins = 0)
      : TFCS1DFunction_HistogramCompactBinEdges<float, uint32_t, float>(
          nbins) {};

  ClassDefNV(TFCS1DFunction_HistogramInt32BinEdges,
             1)  // TFCS1DFunction_HistogramInt32BinEdges
};

#  if defined(__ROOTCLING__) && defined(__FastCaloSimStandAlone__)
#    pragma link C++ class TFCS1DFunction_Numeric < uint8_t, float> + ;
#    pragma link C++ class TFCS1DFunction_Numeric < uint16_t, float> + ;
#    pragma link C++ class TFCS1DFunction_Numeric < uint32_t, float> + ;
#    pragma link C++ class TFCS1DFunction_Numeric < float, float> + ;
#    pragma link C++ class TFCS1DFunction_Numeric < double, float> + ;
#    pragma link C++ class TFCS1DFunction_Numeric < double, double> + ;

#    pragma link C++ class TFCS1DFunction_Array < float> - ;
#    pragma link C++ class TFCS1DFunction_Array < double> - ;
#    pragma link C++ class TFCS1DFunction_Array < uint8_t> - ;
#    pragma link C++ class TFCS1DFunction_Array < uint16_t> - ;
#    pragma link C++ class TFCS1DFunction_Array < uint32_t> - ;

#    pragma link C++ class TFCS1DFunction_HistogramContent < float, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramContent < double, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramContent < double, double> + ;
#    pragma link C++ class TFCS1DFunction_HistogramContent < uint8_t, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramContent < uint16_t, \
        float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramContent < uint32_t, \
        float> + ;

#    pragma link C++ class TFCS1DFunction_HistogramBinEdges < float, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramBinEdges < double, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramBinEdges < double, \
        double> + ;

#    pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float, \
        uint8_t, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float, \
        uint16_t, float> + ;
#    pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float, \
        uint32_t, float> + ;

#    pragma link C++ class TFCS1DFunction_HistogramInt8BinEdges + ;
#    pragma link C++ class TFCS1DFunction_HistogramInt16BinEdges + ;
#    pragma link C++ class TFCS1DFunction_HistogramInt32BinEdges + ;
#    pragma link C++ class TFCS1DFunction_HistogramFloatBinEdges + ;
#    pragma link C++ class TFCS1DFunction_HistogramDoubleBinEdges + ;

#  endif

#endif

#pragma GCC diagnostic pop
