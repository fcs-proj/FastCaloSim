#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSParametrization_h
#  define ISF_FASTCALOSIMEVENT_TFCSParametrization_h

#  include "FastCaloSim/Core/TFCSParametrizationBase.h"

class TFCSParametrization : public ::TFCSParametrizationBase
{
public:
  TFCSParametrization(const char* name = nullptr, const char* title = nullptr);
  void clear();

  auto is_match_pdgid(int id) const -> bool override
  {
    return TestBit(kMatchAllPDGID) || m_pdgid.find(id) != m_pdgid.end();
  };
  auto is_match_Ekin(float Ekin) const -> bool override
  {
    return (Ekin >= m_Ekin_min) && (Ekin < m_Ekin_max);
  };
  auto is_match_eta(float eta) const -> bool override
  {
    return (eta >= m_eta_min) && (eta < m_eta_max);
  };

  auto is_match_all_Ekin() const -> bool override
  {
    return Ekin_min() == init_Ekin_min && Ekin_max() == init_Ekin_max;
  };
  auto is_match_all_eta() const -> bool override
  {
    return eta_min() == init_eta_min && eta_max() == init_eta_max;
  };
  auto is_match_all_Ekin_bin() const -> bool override { return true; };
  auto is_match_all_calosample() const -> bool override { return true; };

  auto pdgid() const -> const std::set<int>& override { return m_pdgid; };
  auto Ekin_nominal() const -> double override { return m_Ekin_nominal; };
  auto Ekin_min() const -> double override { return m_Ekin_min; };
  auto Ekin_max() const -> double override { return m_Ekin_max; };
  auto eta_nominal() const -> double override { return m_eta_nominal; };
  auto eta_min() const -> double override { return m_eta_min; };
  auto eta_max() const -> double override { return m_eta_max; };

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
  auto compare(const TFCSParametrizationBase& ref) const -> bool;

private:
  std::set<int> m_pdgid;
  double m_Ekin_nominal, m_Ekin_min, m_Ekin_max;
  double m_eta_nominal, m_eta_min, m_eta_max;

  ClassDefOverride(TFCSParametrization, 1)  // TFCSParametrization
};

#endif

#pragma GCC diagnostic pop
