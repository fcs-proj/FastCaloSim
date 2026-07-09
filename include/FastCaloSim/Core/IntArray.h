// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include "TArrayI.h"
#include "TObject.h"

class IntArray
    : public TObject
    , public TArrayI
{
public:
  IntArray();
  IntArray(int);
  ~IntArray() override;

private:
  ClassDefOverride(IntArray, 1)
};
