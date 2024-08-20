/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @author Elmar.Ritsch -at- cern.ch
 * @author Robert.Langenberg -at- cern.ch
 * @date March 2017
 * @brief Thread-local storage object used by MagField::AtlasFieldSvc
 */

#ifndef MAGFIELDSERVICES_ATLASFIELDSVCTLS_H
#define MAGFIELDSERVICES_ATLASFIELDSVCTLS_H 1

// MagField includes
#include "BFieldCond.h"
#include "BFieldMeshZR.h"
#include "BFieldZone.h"

namespace MagField
{

/** @class AtlasFieldSvcTLS
 *
 *  @brief Thread-local storage object used by MagField::AtlasFieldSvc
 *
 *  @author Elmar.Ritsch -at- cern.ch
 *  @author Robert.Langenberg -at- cern.ch
 */
struct AtlasFieldSvcTLS
{
  /// Constructor
  AtlasFieldSvcTLS()
      : isInitialized(false)
      , cond(nullptr)
      , cache()
      , cacheZR()
  {
    ;
  }

  /// Is the current AtlasFieldSvcTLS object properly initialized
  bool isInitialized;

  /// Pointer to the conductors in the current field zone (to compute
  /// Biot-Savart component)
  const std::vector<BFieldCond>* cond;

  /// Full 3d field
  BFieldCache cache;
  /// Fast 2d field
  BFieldCacheZR cacheZR;
};

}  // namespace MagField

#endif  // MAGFIELDSERVICES_ATLASFIELDSVCTLS_H
