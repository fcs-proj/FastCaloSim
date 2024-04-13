/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_IntArray_h
#define ISF_FASTCALOSIMEVENT_IntArray_h

#include "TArrayI.h"
#include "TObject.h"

class IntArray : public TObject, public TArrayI {

public:
  IntArray();
  IntArray(int);
  ~IntArray();

private:
  ClassDef(IntArray, 1)
};

#endif
