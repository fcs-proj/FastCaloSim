// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

class G4MagneticField;
class MagFieldPlugin
{
public:
  // Constructor
  MagFieldPlugin() = default;

  // Destructor
  virtual ~MagFieldPlugin() = default;

  virtual auto getField() -> G4MagneticField* = 0;

private:
  MagFieldPlugin(const MagFieldPlugin&) = delete;

  auto operator=(const MagFieldPlugin&) -> MagFieldPlugin& = delete;
};
