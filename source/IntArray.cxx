/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/IntArray.h"

#include "TArrayI.h"

//==========================
//======= IntArray =========
//==========================

IntArray::IntArray(int n)
    : TArrayI(n)
{
}

IntArray::IntArray() {}

IntArray::~IntArray() {}
