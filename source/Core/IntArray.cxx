// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/IntArray.h"

#include "TArrayI.h"

//==========================
//======= IntArray =========
//==========================

IntArray::IntArray(int n)
    : TArrayI(n)
{
}

IntArray::IntArray() = default;

IntArray::~IntArray() = default;
