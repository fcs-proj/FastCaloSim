#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitChain.h"

#include <TClass.h>

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TMath.h"

//=============================================
//======= TFCSLateralShapeParametrizationHitChain =========
//=============================================

TFCSLateralShapeParametrizationHitChain::
    TFCSLateralShapeParametrizationHitChain(const char* name, const char* title)
    : TFCSLateralShapeParametrization(name, title)
    , m_number_of_hits_simul(nullptr)
{
}

TFCSLateralShapeParametrizationHitChain::
    TFCSLateralShapeParametrizationHitChain(
        TFCSLateralShapeParametrizationHitBase* hitsim)
    : TFCSLateralShapeParametrization(
          TString("hit_chain_") + hitsim->GetName(),
          TString("hit chain for ") + hitsim->GetTitle())
    , m_number_of_hits_simul(nullptr)
{
  set_pdgid_Ekin_eta_Ekin_bin_calosample(*hitsim);

  m_chain.push_back(hitsim);
}

unsigned int TFCSLateralShapeParametrizationHitChain::size() const
{
  if (m_number_of_hits_simul)
    return m_chain.size() + 1;
  else
    return m_chain.size();
}

void TFCSLateralShapeParametrizationHitChain::push_back_init(
    const Chain_t::value_type& value)
{
  if (m_ninit == size()) {
    chain().push_back(value);
  } else {
    Chain_t::iterator it(&chain()[m_ninit]);
    chain().insert(it, value);
  }
  ++m_ninit;
}

const TFCSParametrizationBase*
TFCSLateralShapeParametrizationHitChain::operator[](unsigned int ind) const
{
  if (m_number_of_hits_simul) {
    if (ind == 0)
      return m_number_of_hits_simul;
    return m_chain[ind - 1];
  } else {
    return m_chain[ind];
  }
}

TFCSParametrizationBase* TFCSLateralShapeParametrizationHitChain::operator[](
    unsigned int ind)
{
  if (m_number_of_hits_simul) {
    if (ind == 0)
      return m_number_of_hits_simul;
    return m_chain[ind - 1];
  } else {
    return m_chain[ind];
  }
}

void TFCSLateralShapeParametrizationHitChain::set_daughter(
    unsigned int ind, TFCSParametrizationBase* param)
{
  TFCSLateralShapeParametrizationHitBase* param_typed = nullptr;
  if (param != nullptr) {
    if (!param->InheritsFrom(TFCSLateralShapeParametrizationHitBase::Class())) {
      FCS_MSG_ERROR("Wrong class type " << param->IsA()->GetName());
      return;
    }
    param_typed = static_cast<TFCSLateralShapeParametrizationHitBase*>(param);
  }
  if (m_number_of_hits_simul) {
    if (ind == 0)
      m_number_of_hits_simul = param_typed;
    else
      m_chain.at(ind - 1) = param_typed;
  } else {
    m_chain.at(ind) = param_typed;
  }
}

int TFCSLateralShapeParametrizationHitChain::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  // TODO: should we still do it?
  if (m_number_of_hits_simul) {
    int n =
        m_number_of_hits_simul->get_number_of_hits(simulstate, truth, extrapol);
    if (n < 1)
      n = 1;
    return n;
  }
  for (TFCSLateralShapeParametrizationHitBase* hitsim : m_chain) {
    int n = hitsim->get_number_of_hits(simulstate, truth, extrapol);
    if (n > 0)
      return n;
  }
  return 1;
}

float TFCSLateralShapeParametrizationHitChain::getMinWeight() const
{
  for (TFCSLateralShapeParametrizationHitBase* hitsim : m_chain) {
    float weight = hitsim->getMinWeight();
    if (weight > 0.)
      return weight;
  }
  return -1.;
}

float TFCSLateralShapeParametrizationHitChain::getMaxWeight() const
{
  for (TFCSLateralShapeParametrizationHitBase* hitsim : m_chain) {
    float weight = hitsim->getMaxWeight();
    if (weight > 0.)
      return weight;
  }
  return -1.;
}

float TFCSLateralShapeParametrizationHitChain::get_E_hit(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  const int nhits = get_number_of_hits(simulstate, truth, extrapol);
  const int sample = calosample();
  if (sample < 0)
    return 0;
  if (nhits <= 0)
    return simulstate.E(sample);
  const float maxWeight =
      getMaxWeight();  // maxWeight = -1 if  shapeWeight class is not in m_chain

  if (maxWeight > 0)
    return simulstate.E(sample)
        / (maxWeight * nhits);  // maxWeight is used only if shapeWeight class
                                // is in m_chain
  else
    return simulstate.E(sample)
        / nhits;  // Otherwise, old definition of E_hit is used
}

float TFCSLateralShapeParametrizationHitChain::get_sigma2_fluctuation(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  if (m_number_of_hits_simul) {
    double sigma2 = m_number_of_hits_simul->get_sigma2_fluctuation(
        simulstate, truth, extrapol);
    if (sigma2 > 0)
      return sigma2;
  }
  for (TFCSLateralShapeParametrizationHitBase* hitsim : m_chain) {
    double sigma2 = hitsim->get_sigma2_fluctuation(simulstate, truth, extrapol);
    if (sigma2 > 0)
      return sigma2;
  }
  // Limit to factor s_max_sigma2_fluctuation fluctuations
  return s_max_sigma2_fluctuation;
}

void TFCSLateralShapeParametrizationHitChain::PropagateMSGLevel(
    FCS_MSG::Level level) const
{
  for (TFCSLateralShapeParametrizationHitBase* reset : m_chain)
    reset->setLevel(level);
}

FCSReturnCode TFCSLateralShapeParametrizationHitChain::init_hit(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  hit.reset_center();
  if (get_nr_of_init() > 0) {
    FCS_MSG_DEBUG("E(" << calosample() << ")=" << simulstate.E(calosample())
                       << " before init");

    auto initloopend = m_chain.begin() + get_nr_of_init();
    for (auto hititr = m_chain.begin(); hititr != initloopend; ++hititr) {
      TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;

      FCSReturnCode status =
          hitsim->simulate_hit(hit, simulstate, truth, extrapol);

      if (status != FCSSuccess) {
        FCS_MSG_ERROR(
            "TFCSLateralShapeParametrizationHitChain::simulate(): "
            "simulate_hit init call failed");
        return FCSFatal;
      }
    }
  }
  return FCSSuccess;
}

FCSReturnCode TFCSLateralShapeParametrizationHitChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  FCS_MSG::Level old_level = level();
  const bool debug = msgLvl(FCS_MSG::DEBUG);
  const bool verbose = msgLvl(FCS_MSG::VERBOSE);

  // Execute the first get_nr_of_init() simulate calls only once. Used for
  // example to initialize the center position
  TFCSLateralShapeParametrizationHitBase::Hit hit;
  if (init_hit(hit, simulstate, truth, extrapol) != FCSSuccess) {
    FCS_MSG_ERROR("init_hit() failed");
    return FCSFatal;
  }

  int cs = calosample();
  // Initialize hit energy only now, as init loop above might change the layer
  // energy
  const float Elayer = simulstate.E(cs);
  if (Elayer == 0) {
    FCS_MSG_VERBOSE("Elayer=0, nothing to do");
    return FCSSuccess;
  }
  const float Ehit = get_E_hit(simulstate, truth, extrapol);
  if (Ehit * Elayer <= 0) {
    FCS_MSG_ERROR(
        "TFCSLateralShapeParametrizationHitChain::simulate(): Ehit "
        "and Elayer have different sign! Ehit="
        << Ehit << " Elayer=" << Elayer);
    return FCSFatal;
  }

  // Call get_number_of_hits() only once inside get_E_hit(...),
  // as it could contain a random number
  int nhit = TMath::Nint(Elayer / Ehit);
  if (nhit < 1)
    nhit = 1;

  float sumEhit = 0;

  if (debug) {
    PropagateMSGLevel(old_level);
    FCS_MSG_DEBUG("E(" << cs << ")=" << simulstate.E(cs) << " #hits~" << nhit);
  }

  {
    auto hitloopstart = m_chain.begin() + get_nr_of_init();
    int ihit = 0;
    int ifail = 0;
    int itotalfail = 0;
    int retry_warning = 1;
    int retry = 0;
    do {
      hit.reset();
      hit.set_E(Ehit);
      bool failed = false;
      if (debug)
        if (ihit == 2)
          if (!verbose) {
            // Switch debug output back to INFO to avoid huge logs, but keep
            // full log in verbose
            PropagateMSGLevel(FCS_MSG::INFO);
          }
      for (auto hititr = hitloopstart; hititr != m_chain.end(); ++hititr) {
        TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;

        FCSReturnCode status =
            hitsim->simulate_hit(hit, simulstate, truth, extrapol);

        if (status == FCSSuccess)
          continue;
        if (status == FCSFatal) {
          if (debug)
            PropagateMSGLevel(old_level);
          return FCSFatal;
        }
        failed = true;
        ++ifail;
        ++itotalfail;
        retry = status - FCSRetry;
        retry_warning = retry >> 1;
        if (retry_warning < 1)
          retry_warning = 1;
        break;
      }
      if (!failed) {
        ifail = 0;
        sumEhit += hit.E();
        ++ihit;

        if (((ihit == 20 * nhit) || (ihit == 100 * nhit)) && ihit >= 100) {
          FCS_MSG_DEBUG(
              "TFCSLateralShapeParametrizationHitChain::simulate(): Iterated "
              << ihit << " times, expected " << nhit << " times. Deposited E("
              << cs << ")=" << sumEhit << " expected E=" << Elayer);
        }
        if (ihit >= 1000 * nhit && ihit >= 1000) {
          FCS_MSG_DEBUG(
              "TFCSLateralShapeParametrizationHitChain::simulate():"
              " Aborting hit chain, iterated "
              << ihit << " times, expected " << nhit << " times. Deposited E("
              << cs << ")=" << sumEhit << " expected E=" << Elayer
              << ", caused by:");
          if (debug)
            Print();
          break;
        }
      } else {
        if (ifail >= retry) {
          FCS_MSG_ERROR(
              "TFCSLateralShapeParametrizationHitChain::simulate(): "
              "simulate_hit call failed after "
              << ifail << "/" << retry
              << "retries, total fails=" << itotalfail);
          if (debug)
            PropagateMSGLevel(old_level);
          return FCSFatal;
        }
        if (ifail >= retry_warning) {
          FCS_MSG_WARNING(
              "TFCSLateralShapeParametrizationHitChain::simulate():"
              " retry simulate_hit calls "
              << ifail << "/" << retry << ", total fails=" << itotalfail);
        }
      }
    } while (std::abs(sumEhit) < std::abs(Elayer));
  }

  if (debug)
    PropagateMSGLevel(old_level);
  return FCSSuccess;
}

void TFCSLateralShapeParametrizationHitChain::Print(Option_t* option) const
{
  TFCSLateralShapeParametrization::Print(option);
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");

  if (m_number_of_hits_simul) {
    if (longprint)
      FCS_MSG_INFO(optprint << "#:Number of hits simulation:");
    m_number_of_hits_simul->Print(opt + "#:");
  }
  if (longprint && get_nr_of_init() > 0)
    FCS_MSG_INFO(optprint << "> Simulation init chain:");
  auto hitloopstart = m_chain.begin() + get_nr_of_init();
  for (auto hititr = m_chain.begin(); hititr != hitloopstart; ++hititr) {
    TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;
    hitsim->Print(opt + "> ");
  }
  if (longprint)
    FCS_MSG_INFO(optprint << "- Simulation chain:");
  char count = 'A';
  for (auto hititr = hitloopstart; hititr != m_chain.end(); ++hititr) {
    TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;
    hitsim->Print(opt + count + " ");
    count++;
  }
}

#pragma GCC diagnostic pop
