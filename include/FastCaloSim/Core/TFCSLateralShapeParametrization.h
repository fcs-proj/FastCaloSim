// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSLateralShapeParametrization_h
#define TFCSLateralShapeParametrization_h

#include "FastCaloSim/Core/TFCSParametrization.h"

namespace FastCaloSim::Core
{
class TFCSLateralShapeParametrization : public TFCSParametrization
{
public:
  TFCSLateralShapeParametrization(const char* name = nullptr,
                                  const char* title = nullptr);

  bool is_match_Ekin_bin(int bin) const override
  {
    if (Ekin_bin() == -1)
      return true;
    return bin == Ekin_bin();
  };
  bool is_match_calosample(int calosample) const override
  {
    return calosample == m_calosample;
  };

  virtual bool is_match_all_Ekin_bin() const override
  {
    if (Ekin_bin() == -1)
      return true;
    return false;
  };
  virtual bool is_match_all_calosample() const override { return false; };

  int Ekin_bin() const { return m_Ekin_bin; };
  void set_Ekin_bin(int bin);

  int calosample() const { return m_calosample; };
  void set_calosample(int cs);

  virtual void set_pdgid_Ekin_eta_Ekin_bin_calosample(
      const TFCSLateralShapeParametrization& ref);

  void Print(Option_t* option = "") const override;

protected:
  bool compare(const TFCSParametrizationBase& ref) const;

private:
  int m_Ekin_bin;
  int m_calosample;

  ClassDefOverride(TFCSLateralShapeParametrization,
                   1)  // TFCSLateralShapeParametrization
};
}  // namespace FastCaloSim::Core
#endif
