#ifndef MAGFIELDPLUGIN_H
#define MAGFIELDPLUGIN_H
class G4MagneticField;
class MagFieldPlugin
{
public:
  // Constructor
  MagFieldPlugin() = default;

  // Destructor
  virtual ~MagFieldPlugin() = default;

  MagFieldPlugin(const MagFieldPlugin&) = delete;
  auto operator=(const MagFieldPlugin&) -> MagFieldPlugin& = delete;
  MagFieldPlugin(MagFieldPlugin&&) = delete;
  auto operator=(MagFieldPlugin&&) -> MagFieldPlugin& = delete;

  virtual auto getField() -> G4MagneticField* = 0;
};

#endif
