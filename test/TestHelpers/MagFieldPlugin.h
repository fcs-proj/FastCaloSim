#pragma once

class G4MagneticField;
class MagFieldPlugin
{
public:
  // Constructor
  MagFieldPlugin() = default;

  // Destructor
  virtual ~MagFieldPlugin() = default;

  virtual G4MagneticField* getField() = 0;

private:
  MagFieldPlugin(const MagFieldPlugin&) = delete;

  MagFieldPlugin& operator=(const MagFieldPlugin&) = delete;
};
