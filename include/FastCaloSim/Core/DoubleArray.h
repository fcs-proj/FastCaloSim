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
  explicit DoubleArray(int /*n*/);
  ~DoubleArray() override;
  DoubleArray(const DoubleArray&) = default;
  DoubleArray(DoubleArray&&) = default;
  auto operator=(const DoubleArray&) -> DoubleArray& = default;
  auto operator=(DoubleArray&&) -> DoubleArray& = default;

private:
  ClassDefOverride(DoubleArray, 1)
};
