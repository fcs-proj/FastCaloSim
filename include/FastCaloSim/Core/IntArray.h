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
  explicit IntArray(int /*n*/);
  ~IntArray() override;
  IntArray(const IntArray&) = default;
  IntArray(IntArray&&) = default;
  auto operator=(const IntArray&) -> IntArray& = default;
  auto operator=(IntArray&&) -> IntArray& = default;

private:
  ClassDefOverride(IntArray, 1)
};
