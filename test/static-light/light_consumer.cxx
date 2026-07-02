// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/MLogging.h"

// A minimal consumer of the lightweight logging mixin, deliberately using
// nothing else from FastCaloSim: pulling in the logger must not require the
// heavy dependency closure (see CMakeLists.txt in this directory).
namespace
{
class LightLogger : public ISF_FCS::MLogging
{
public:
  void greet() const { FCS_MSG_INFO("light logging path works"); }
};
}  // namespace

// Called by the driver executable.
auto fcs_emit_light_log_line() -> int
{
  LightLogger logger;
  logger.setLevel(FCS_MSG::DEBUG);
  logger.greet();
  return 0;
}
