#ifndef _MAGFIELDPLUGIN_H_
#define _MAGFIELDPLUGIN_H_
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

#endif
