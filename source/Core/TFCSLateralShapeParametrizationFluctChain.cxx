// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFluctChain.h"

#include "CLHEP/Random/RandGaussZiggurat.h"
#include "FastCaloSim/Core/GaussZiggurat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Definitions/FloatCompare.h"

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

auto TFCSLateralShapeParametrizationFluctChain::get_E_hit(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> double
{
  const double sigma2 = get_sigma2_fluctuation(simulstate, truth, extrapol);
  if (sigma2 <= 0 || !has_calosample()) {
    return -1.;
  }
  const std::size_t sample = calosample();
  const double maxWeight = getMaxWeight();

  if (maxWeight > 0) {
    return simulstate.E(sample) * sigma2 / maxWeight;
  }
  return simulstate.E(sample) * sigma2;
}

auto TFCSLateralShapeParametrizationFluctChain::simulate(
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol) const -> FCSReturnCode
{
  FastCaloSim::MSG::Level const old_level = level();
  const bool debug = msgLvl(FastCaloSim::MSG::DEBUG);

  // Execute the first get_nr_of_init() simulate calls only once. Used for
  // example to initialize the center position
  TFCSLateralShapeParametrizationHitBase::Hit hit;
  if (init_hit(hit, simulstate, truth, extrapol) != FCSSuccess) {
    FCS_MSG_ERROR("init_hit() failed");
    return FCSFatal;
  }

  // Initialize hit energy only now, as init loop above might change the layer
  // energy
  const double Elayer = simulstate.E(calosample());
  if (FastCaloSim::exactly_equal(Elayer, 0.0)) {
    FCS_MSG_VERBOSE("Elayer=0, nothing to do");
    return FCSSuccess;
  }

  // Call get_sigma2_fluctuation only once, as it could contain a random number
  double sigma2 = get_sigma2_fluctuation(simulstate, truth, extrapol);
  if (sigma2 >= s_max_sigma2_fluctuation) {
    FCS_MSG_ERROR(
        "TFCSLateralShapeParametrizationFluctChain::simulate(): "
        "fluctuation of hits could not be calculated");
    return FCSFatal;
  }

  // Limit to relative precision of 10^-4=0.01%. ATLAS calorimeter is ~0.1% at
  // best
  sigma2 = std::max(sigma2, 1e-8);

  // Make a good guess of the needed hit energy, assuming all hits would have
  // the same energy
  const double Eavghit = get_E_hit(simulstate, truth, extrapol);
  const double absEavghit_tenth = std::abs(Eavghit / 10);
  // The number of hits expected with that energy. As in
  // TFCSLateralShapeParametrizationHitChain, the loop is aborted after 1000
  // times as many: hits that are all rejected (deposit no energy) would
  // otherwise never reach the precision
  const double nhit_expected =
      std::max(std::round(std::abs(Elayer / Eavghit)), 1.0);
  const double max_hits = std::max(1000 * nhit_expected, 1000.0);
  double sumEhit = 0;
  double error2_sumEhit = 0;
  double error2 = 2 * s_max_sigma2_fluctuation;

  if (debug) {
    PropagateMSGLevel(old_level);
    FCS_MSG_DEBUG("E(" << calosample() << ")=" << Elayer
                       << " sigma2=" << sigma2);
  }

  auto hitloopstart = m_chain.begin() + get_nr_of_init();
  std::size_t ihit = 0;
  int ifail = 0;
  int itotalfail = 0;
  int retry_warning = 1;
  int retry = 0;
  do {
    hit.reset();
    do {
      hit.set_E(
          FastCaloSim::gauss_ziggurat(simulstate.randomEngine(),
                                      Eavghit,
                                      static_cast<double>(m_RMS) * Eavghit));
    } while (std::abs(hit.E()) < absEavghit_tenth);
    bool failed = false;
    if (debug) {
      if (ihit == 2) {
        // Switch debug output back to INFO to avoid huge logs
        PropagateMSGLevel(FastCaloSim::MSG::INFO);
      }
    }
    for (auto hititr = hitloopstart; hititr != m_chain.end(); ++hititr) {
      TFCSLateralShapeParametrizationHitBase* hitsim = *hititr;

      FCSReturnCode const status =
          hitsim->simulate_hit(hit, simulstate, truth, extrapol);

      if (status == FCSSuccess) {
        continue;
      }
      if (status == FCSFatal) {
        if (debug) {
          PropagateMSGLevel(old_level);
        }
        return FCSFatal;
      }
      failed = true;
      ++ifail;
      ++itotalfail;
      retry = status.max_retries();
      retry_warning = retry / 2;
      retry_warning = std::max(retry_warning, 1);
      break;
    }
    if (!failed) {
      ifail = 0;
      ++ihit;
      const double Ehit = hit.E();
      sumEhit += Ehit;
      double const sumEhit2 = sumEhit * sumEhit;
      error2_sumEhit += Ehit * Ehit;
      if (sumEhit2 > 0) {
        error2 = error2_sumEhit / sumEhit2;
      }
      if (static_cast<double>(ihit) >= max_hits && error2 > sigma2) {
        FCS_MSG_WARNING(
            "TFCSLateralShapeParametrizationFluctChain::simulate(): Aborting "
            "hit chain, iterated "
            << ihit << " times, expected " << nhit_expected
            << " times. Deposited E(" << calosample() << ")=" << sumEhit
            << " expected E=" << Elayer);
        break;
      }
    } else {
      if (ifail >= retry) {
        FCS_MSG_ERROR(
            "TFCSLateralShapeParametrizationFluctChain::simulate(): "
            "simulate_hit call failed after "
            << ifail << "/" << retry << "retries, total fails=" << itotalfail);
        if (debug) {
          PropagateMSGLevel(old_level);
        }
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
                       << "+-" << std::sqrt(error2_sumEhit) << " ~ "
                       << std::sqrt(error2_sumEhit) / sumEhit * 100
                       << "% rel error^2=" << error2 << " sigma^2=" << sigma2
                       << " ~ " << std::sqrt(sigma2) * 100 << "% hits=" << ihit
                       << " fail=" << itotalfail);
  }

  return FCSSuccess;
}

void TFCSLateralShapeParametrizationFluctChain::Print(Option_t* option) const
{
  const auto [shortprint, longprint, optprint] = print_options(option);
  TFCSLateralShapeParametrizationHitChain::Print(option);

  if (longprint) {
    FCS_MSG_INFO(optprint << "  hit energy fluctuation RMS=" << m_RMS);
  }
}
