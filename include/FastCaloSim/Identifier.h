/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef IdentifierStandAlone
#define IdentifierStandAlone

// Not really sure why this is needed.
// If it isn't there, then we get 
// error: 'Long64_t' does not name a type
// in Identifier.h, only I haven't changed that.
// Some dependancy isn't defined quite right.
#include <RtypesCore.h>

class Identifier
{
  public:
  typedef Long64_t value_type;
  
  typedef enum bit_defs_enum
  {
    NBITS = sizeof(value_type) * 8, // bits per byte
    MAX_BIT = (static_cast<value_type>(1) << (NBITS - 1)),
    ALL_BITS = ~(static_cast<value_type>(0))
  } bit_defs;
  
  typedef enum max_value_type_enum {
    //max_value = 0xFFFFFFFFFFFFFFFFULL
    max_value = ~(static_cast<value_type>(0))
  } max_value_type;

  Identifier():m_id(max_value) {};
  Identifier(const Identifier& value):m_id(value.m_id) {};
  Identifier(value_type value):m_id(value) {};
  
  operator value_type() const { return m_id; }

  // cppcheck-suppress operatorEqVarError; false positive
  Identifier& operator = (const Identifier& old) {m_id=old;return (*this);};
  Identifier& operator = (value_type value) {m_id=value;return (*this);};
  bool operator == (const Identifier& other) const {return (m_id == other.m_id);}
  bool operator != (const Identifier& other) const {return (m_id != other.m_id);}
  bool operator < (const Identifier& other) const {return (m_id < other.m_id);}
  bool operator > (const Identifier& other) const {return (m_id > other.m_id);}
  bool operator <= (const Identifier& other) const {return (m_id <= other.m_id);}
  bool operator >= (const Identifier& other) const {return (m_id >= other.m_id);}
  bool operator == (Identifier::value_type other) const {return (m_id == other);}
  bool operator != (Identifier::value_type other) const {return (m_id != other);}
  
  protected:
  value_type m_id;
  
};
  
#endif

