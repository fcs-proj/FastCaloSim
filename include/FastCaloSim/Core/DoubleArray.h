// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "TArrayD.h"
#include "TObject.h"

class DoubleArray
    : public TObject
    , public TArrayD
{
public:
  DoubleArray();
  DoubleArray(int);
  ~DoubleArray() override;

private:
  ClassDef(DoubleArray, 1)
};
