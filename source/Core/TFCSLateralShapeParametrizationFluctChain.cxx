#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFluctChain.h"

#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TMath.h"

//=============================================
//======= TFCSLateralShapeParametrizationFluctChain =========
//=============================================

TFCSLateralShapeParametrizationFluctChain::
    TFCSLateralShapeParametrizationFluctChain(const char* name,
                                              const char* title,
                                              float RMS)
    : TFCSLateralShapeParametrizationHitChain(name, title)
    , m_RMS(RMS)
{
}

TFCSLateralShapeParametrizationFluctChain::
    TFCSLateralShapeParametrizationFluctChain(
        TFCSLateralShapeParametrizationHitBase* hitsim)
    : TFCSLateralShapeParametrizationHitChain(hitsim)
{
}

float TFCSLateralShapeParametrizationFluctChain::get_E_hit(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  const float sigma2 = get_sigma2_fluctuation(simulstate, truth, extrapol);
  const int sample = calosample();
  if (sigma2 <= 0 || sample < 0)
    return -1.;
  const float maxWeight = getMaxWeight();

  if (maxWeight > 0)
    return simulstate.E(sample) * sigma2 / maxWeight;
  else
    return simulstate.E(sample) * sigma2;
}

FCSReturnCode TFCSLateralShapeParametrizationFluctChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const
{
  FCS_MSG::Level old_level = level();
  const bool debug = msgLvl(FCS_MSG::DEBUG);

  // Execute the first get_nr_of_init() simulate calls only once. Used for
  // example to initialize the center position
  TFCSLateralShapeParametrizationHitBase::Hit hit;
  if (init_hit(hit, simulstate, truth, extrapol) != FCSSuccess) {
    FCS_MSG_ERROR("init_hit() failed");
    return FCSFatal;
  }

  // Initialize hit energy only now, as init loop above might change the layer
  // energy
  const float Elayer = simulstate.E(calosample());
  if (Elayer == 0) {
    FCS_MSG_VERBOSE("Elayer=0, nothing to do");
    return FCSSuccess;
  }

  // Call get_sigma2_fluctuation only once, as it could contain a random number
  float sigma2 = get_sigma2_fluctuation(simulstate, truth, extrapol);
  if (sigma2 >= s_max_sigma2_fluctuation) {
    FCS_MSG_ERROR(
        "TFCSLateralShapeParametrizationFluctChain::simulate(): "
        "fluctuation of hits could not be calculated");
    return FCSFatal;
  }

  // Limit to relative precision of 10^-4=0.01%. ATLAS calorimeter is ~0.1% at
  // best
  if (sigma2 < 1e-8)
    sigma2 = 1e-8;

  // Make a good guess of the needed hit energy, assuming all hits would have
  // the same energy
  const float Eavghit = get_E_hit(simulstate, truth, extrapol);
  const float absEavghit_tenth = std::abs(Eavghit / 10);
  float sumEhit = 0;
  float error2_sumEhit = 0;
  float error2 = 2 * s_max_sigma2_fluctuation;

  if (debug) {
    PropagateMSGLevel(old_level);
    FCS_MSG_DEBUG("E(" << calosample() << ")=" << Elayer
                       << " sigma2=" << sigma2);
  }

  auto hitloopstart = m_chain.begin() + get_nr_of_init();
  int ihit = 0;
  int ifail = 0;
  int itotalfail = 0;
  int retry_warning = 1;
  int retry = 0;
  do {
    hit.reset();
    // hit.E()=Eavghit;
    do {
      hit.set_E(CLHEP::RandGaussZiggurat::shoot(
          simulstate.randomEngine(), Eavghit, m_RMS * Eavghit));
    } while (std::abs(hit.E()) < absEavghit_tenth);
    bool failed = false;
    if (debug)
      if (ihit == 2) {
        // Switch debug output back to INFO to avoid huge logs
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
      ++ihit;
      const float Ehit = hit.E();
      sumEhit += Ehit;
      float sumEhit2 = sumEhit * sumEhit;
      error2_sumEhit += Ehit * Ehit;
      if (sumEhit2 > 0)
        error2 = error2_sumEhit / sumEhit2;
    } else {
      if (ifail >= retry) {
        FCS_MSG_ERROR(
            "TFCSLateralShapeParametrizationFluctChain::simulate(): "
            "simulate_hit call failed after "
            << ifail << "/" << retry << "retries, total fails=" << itotalfail);
        if (debug)
          PropagateMSGLevel(old_level);
        return FCSFatal;
      }
      if (ifail >= retry_warning) {
        FCS_MSG_WARNING(
            "TFCSLateralShapeParametrizationFluctChain::simulate():"
            " retry simulate_hit calls "
            << ifail << "/" << retry << ", total fails=" << itotalfail);
      }
    }
  } while (error2 > sigma2);

  if (debug) {
    PropagateMSGLevel(old_level);
    FCS_MSG_DEBUG("E(" << calosample() << ")=" << Elayer << " sumE=" << sumEhit
                       << "+-" << TMath::Sqrt(error2_sumEhit) << " ~ "
                       << TMath::Sqrt(error2_sumEhit) / sumEhit * 100
                       << "% rel error^2=" << error2 << " sigma^2=" << sigma2
                       << " ~ " << TMath::Sqrt(sigma2) * 100
                       << "% hits=" << ihit << " fail=" << itotalfail);
  }

  return FCSSuccess;
}

void TFCSLateralShapeParametrizationFluctChain::Print(Option_t* option) const
{
  TString opt(option);
  bool shortprint = opt.Index("short") >= 0;
  bool longprint =
      msgLvl(FCS_MSG::DEBUG) || (msgLvl(FCS_MSG::INFO) && !shortprint);
  TString optprint = opt;
  optprint.ReplaceAll("short", "");
  TFCSLateralShapeParametrizationHitChain::Print(option);

  if (longprint)
    FCS_MSG_INFO(optprint << "  hit energy fluctuation RMS=" << m_RMS);
}

#pragma GCC diagnostic pop
