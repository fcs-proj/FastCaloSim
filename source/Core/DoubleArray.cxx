// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/DoubleArray.h"

#include "TArrayD.h"

//==========================
//======= DoubleArray =========
//==========================
using namespace FastCaloSim::Core;

DoubleArray::DoubleArray(int n)
    : TArrayD(n)
{
}

DoubleArray::DoubleArray() {}

DoubleArray::~DoubleArray() {}
