/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/DoubleArray.h"

#include "TArrayD.h"

//==========================
//======= DoubleArray =========
//==========================

DoubleArray::DoubleArray(int n)
    : TArrayD(n)
{
}

DoubleArray::DoubleArray() {}

DoubleArray::~DoubleArray() {}
