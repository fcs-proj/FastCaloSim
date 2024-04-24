#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrization_h
#  define ISF_FASTCALOSIMEVENT_TFCSParametrization_h

#  include "FastCaloSim/Core/TFCSParametrizationBase.h"

class TFCSParametrization : public ::TFCSParametrizationBase
{
public:
  TFCSParametrization(const char* name = nullptr, const char* title = nullptr);
  void clear();

  virtual bool is_match_pdgid(int id) const override
  {
    return TestBit(kMatchAllPDGID) || m_pdgid.find(id) != m_pdgid.end();
  };
  virtual bool is_match_Ekin(float Ekin) const override
  {
    return (Ekin >= m_Ekin_min) && (Ekin < m_Ekin_max);
  };
  virtual bool is_match_eta(float eta) const override
  {
    return (eta >= m_eta_min) && (eta < m_eta_max);
  };

  virtual bool is_match_all_Ekin() const override
  {
    return Ekin_min() == init_Ekin_min && Ekin_max() == init_Ekin_max;
  };
  virtual bool is_match_all_eta() const override
  {
    return eta_min() == init_eta_min && eta_max() == init_eta_max;
  };
  virtual bool is_match_all_Ekin_bin() const override { return true; };
  virtual bool is_match_all_calosample() const override { return true; };

  const std::set<int>& pdgid() const override { return m_pdgid; };
  double Ekin_nominal() const override { return m_Ekin_nominal; };
  double Ekin_min() const override { return m_Ekin_min; };
  double Ekin_max() const override { return m_Ekin_max; };
  double eta_nominal() const override { return m_eta_nominal; };
  double eta_min() const override { return m_eta_min; };
  double eta_max() const override { return m_eta_max; };

  virtual void set_pdgid(int id);
  virtual void set_pdgid(const std::set<int>& ids);
  virtual void add_pdgid(int id);
  virtual void clear_pdgid();

  virtual void set_Ekin_nominal(double min);
  virtual void set_Ekin_min(double min);
  virtual void set_Ekin_max(double max);
  virtual void set_eta_nominal(double min);
  virtual void set_eta_min(double min);
  virtual void set_eta_max(double max);

  virtual void set_Ekin(const TFCSParametrizationBase& ref);
  virtual void set_eta(const TFCSParametrizationBase& ref);
  virtual void set_Ekin_eta(const TFCSParametrizationBase& ref);
  virtual void set_pdgid_Ekin_eta(const TFCSParametrizationBase& ref);

protected:
  bool compare(const TFCSParametrizationBase& ref) const;

private:
  std::set<int> m_pdgid;
  double m_Ekin_nominal, m_Ekin_min, m_Ekin_max;
  double m_eta_nominal, m_eta_min, m_eta_max;

  ClassDefOverride(TFCSParametrization, 1)  // TFCSParametrization
};

#endif

#pragma GCC diagnostic pop
