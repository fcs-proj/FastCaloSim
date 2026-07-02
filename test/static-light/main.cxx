// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

// Provided by the shared library that absorbed the lightweight logging part
// of the FastCaloSim archive.
auto fcs_emit_light_log_line() -> int;

auto main() -> int
{
  return fcs_emit_light_log_line();
}
