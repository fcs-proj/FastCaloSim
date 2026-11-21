#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitChain.h"

#include <TClass.h>

#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"

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
          (std::string("hit_chain_") + hitsim->GetName()).c_str(),
          (std::string("hit chain for ") + hitsim->GetTitle()).c_str())
    , m_number_of_hits_simul(nullptr)
{
  set_pdgid_Ekin_eta_Ekin_bin_calosample(*hitsim);

  m_chain.push_back(hitsim);
}

auto TFCSLateralShapeParametrizationHitChain::size() const -> unsigned int
{
  if (m_number_of_hits_simul)
    return static_cast<unsigned int>(m_chain.size() + 1);
  else
    return static_cast<unsigned int>(m_chain.size());
}

void TFCSLateralShapeParametrizationHitChain::push_back_init(
    const Chain_t::value_type& value)
{
  if (m_ninit == size()) {
    chain().push_back(value);
  } else {
    const auto it = chain().begin() + m_ninit;
    chain().insert(it, value);
  }
  ++m_ninit;
}

auto TFCSLateralShapeParametrizationHitChain::operator[](unsigned int ind) const
    -> const TFCSParametrizationBase*
{
  if (m_number_of_hits_simul) {
    if (ind == 0)
      return m_number_of_hits_simul;
    return m_chain[ind - 1];
  } else {
    return m_chain[ind];
  }
}

auto TFCSLateralShapeParametrizationHitChain::operator[](unsigned int ind)
    -> TFCSParametrizationBase*
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

auto TFCSLateralShapeParametrizationHitChain::get_number_of_hits(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> int
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

auto TFCSLateralShapeParametrizationHitChain::getMinWeight() const -> double
{
  for (TFCSLateralShapeParametrizationHitBase* hitsim : m_chain) {
    double weight = hitsim->getMinWeight();
    if (weight > 0.)
      return weight;
  }
  return -1.;
}

auto TFCSLateralShapeParametrizationHitChain::getMaxWeight() const -> double
{
  for (TFCSLateralShapeParametrizationHitBase* hitsim : m_chain) {
    double weight = hitsim->getMaxWeight();
    if (weight > 0.)
      return weight;
  }
  return -1.;
}

auto TFCSLateralShapeParametrizationHitChain::get_E_hit(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> double
{
  const int nhits = get_number_of_hits(simulstate, truth, extrapol);

  // Store nhits in the simulation state as the get function may only be
  // called once per event
  simulstate.setAuxInfo<int>("FCSHitChainNHits"_FCShash, nhits);

  const int sample = calosample();
  if (sample < 0)
    return 0;
  if (nhits <= 0)
    return simulstate.E(sample);
  const double maxWeight =
      getMaxWeight();  // maxWeight = -1 if  shapeWeight class is not in m_chain

  if (maxWeight > 0)
    return simulstate.E(sample)
        / (maxWeight * nhits);  // maxWeight is used only if shapeWeight class
                                // is in m_chain
  else
    return simulstate.E(sample)
        / nhits;  // Otherwise, old definition of E_hit is used
}

auto TFCSLateralShapeParametrizationHitChain::get_sigma2_fluctuation(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> double
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
    FastCaloSim::MSG::Level level) const
{
  for (TFCSLateralShapeParametrizationHitBase* reset : m_chain)
    reset->setLevel(level);
}

auto TFCSLateralShapeParametrizationHitChain::init_hit(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
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

auto TFCSLateralShapeParametrizationHitChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  FastCaloSim::MSG::Level old_level = level();
  const bool debug = msgLvl(FastCaloSim::MSG::DEBUG);
  const bool verbose = msgLvl(FastCaloSim::MSG::VERBOSE);

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
  const double Elayer = simulstate.E(cs);
  if (FastCaloSim::exactly_equal(Elayer, 0.0)) {
    FCS_MSG_VERBOSE("Elayer=0, nothing to do");
    return FCSSuccess;
  }
  const double Ehit = get_E_hit(simulstate, truth, extrapol);
  if (Ehit * Elayer <= 0) {
    FCS_MSG_ERROR(
        "TFCSLateralShapeParametrizationHitChain::simulate(): Ehit "
        "and Elayer have different sign! Ehit="
        << Ehit << " Elayer=" << Elayer);
    return FCSFatal;
  }

  // Call get_number_of_hits() only once inside get_E_hit(...),
  // as it could contain a random number
  // TODO: should read simulate value!
  int nhit_signed = static_cast<int>(std::lrint(Elayer / Ehit));
  if (nhit_signed < 1)
    nhit_signed = 1;
  unsigned int nhit = static_cast<unsigned int>(nhit_signed);

  double sumEhit = 0;

  if (debug) {
    PropagateMSGLevel(old_level);
    FCS_MSG_DEBUG("E(" << cs << ")=" << simulstate.E(cs) << " #hits~" << nhit);
  }
  {
    auto hitloopstart = m_chain.begin() + get_nr_of_init();
    hit.set_idx(0);
    int ifail = 0;
    int itotalfail = 0;
    int retry_warning = 1;
    int retry = 0;
    bool failed = false;
    do {
      hit.reset();
      hit.set_E(Ehit);
      failed = false;
      if (debug)
        if (hit.idx() == 2)
          if (!verbose) {
            // Switch debug output back to INFO to avoid huge logs, but keep
            // full log in verbose
            PropagateMSGLevel(FastCaloSim::MSG::INFO);
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
        hit.set_idx(hit.idx() + 1);

        // TODO: Different signedness:
        if (((hit.idx() == 20 * nhit) || (hit.idx() == 100 * nhit))
            && hit.idx() >= 100)
        {
          FCS_MSG_DEBUG(
              "TFCSLateralShapeParametrizationHitChain::simulate(): Iterated "
              << hit.idx() << " times, expected " << nhit
              << " times. Deposited E(" << cs << ")=" << sumEhit
              << " expected E=" << Elayer);
        }
        if (hit.idx() >= 1000 * nhit && hit.idx() >= 1000) {
          FCS_MSG_DEBUG(
              "TFCSLateralShapeParametrizationHitChain::simulate():"
              " Aborting hit chain, iterated "
              << hit.idx() << " times, expected " << nhit
              << " times. Deposited E(" << cs << ")=" << sumEhit
              << " expected E=" << Elayer << ", caused by:");
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
    } while (
        !check_all_hits_simulated(hit, simulstate, truth, extrapol, !failed));
  }

  if (debug)
    PropagateMSGLevel(old_level);
  return FCSSuccess;
}

auto TFCSLateralShapeParametrizationHitChain::check_all_hits_simulated(
    TFCSLateralShapeParametrizationHitBase::Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol,
    bool success) const -> bool
{
  (void)truth;  // unused parameter
  (void)extrapol;  // unused parameter

  if (!success) {
    // The hit simulation failed => Nothing changed => Still not done
    return false;
  }

  // Called for every hit: read and update the sum with a single lookup
  auto [stored_sum, was_stored] =
      simulstate.auxInfoEntry("FCSHitChainEnergySum"_FCShash);
  double sumEhit = 0;
  if (was_stored && hit.idx() > 1) {
    sumEhit = static_cast<double>(stored_sum);
  }

  sumEhit += hit.E();

  stored_sum.set<double>(sumEhit);

  const double Elayer = simulstate.E(calosample());
  return (std::abs(sumEhit) >= std::abs(Elayer));
}

void TFCSLateralShapeParametrizationHitChain::Print(Option_t* option) const
{
  TFCSLateralShapeParametrization::Print(option);
  const std::string opt = option != nullptr ? option : "";
  const auto [shortprint, longprint, optprint] = print_options(option);

  if (m_number_of_hits_simul) {
    if (longprint)
      FCS_MSG_INFO(optprint << "#:Number of hits simulation:");
    m_number_of_hits_simul->Print((opt + "#:").c_str());
  }
  if (longprint && get_nr_of_init() > 0)
    FCS_MSG_INFO(optprint << "> Simulation init chain:");
  auto hitloopstart = m_chain.begin() + get_nr_of_init();
  for (auto hititr = m_chain.begin(); hititr != hitloopstart; ++hititr) {
    TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;
    hitsim->Print((opt + "> ").c_str());
  }
  if (longprint)
    FCS_MSG_INFO(optprint << "- Simulation chain:");
  char count = 'A';
  for (auto hititr = hitloopstart; hititr != m_chain.end(); ++hititr) {
    TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;
    hitsim->Print((opt + count + " ").c_str());
    count++;
  }
}

#pragma GCC diagnostic pop
