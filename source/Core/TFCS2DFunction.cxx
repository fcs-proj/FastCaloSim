// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <cmath>
#include <iostream>

#include "FastCaloSim/Core/TFCS2DFunction.h"

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Definitions/FloatCompare.h"
#include "TH2.h"
#include "TRandom.h"

//=============================================
//======= TFCS2DFunction =========
//=============================================

void TFCS2DFunction::rnd_to_fct(float value[], const float rnd[]) const
{
  double valuex = 0;
  double valuey = 0;
  rnd_to_fct(
      valuex, valuey, static_cast<double>(rnd[0]), static_cast<double>(rnd[1]));
  value[0] = static_cast<float>(valuex);
  value[1] = static_cast<float>(valuey);
}

//================================================================================================================================

auto TFCS2DFunction::CheckAndIntegrate2DHistogram(
    const TH2* hist, std::vector<double>& integral_vec, int& first, int& last)
    -> double
{
  FastCaloSim::MLogging logger;
  Int_t nbinsx = hist->GetNbinsX();
  Int_t nbinsy = hist->GetNbinsY();
  Int_t nbins = nbinsx * nbinsy;

  double integral = 0.0;
  const double hint = hist->Integral();
  integral_vec.resize(static_cast<std::size_t>(nbins));

  for (int ix = 1; ix <= nbinsx; ix++) {
    for (int iy = 1; iy <= nbinsy; iy++) {
      const auto globalbin =
          static_cast<std::size_t>((ix - 1) * nbinsy + iy - 1);
      double binval = hist->GetBinContent(ix, iy);
      if (binval < 0) {
        // Can't work if a bin is negative, forcing bins to 0 in this case
        double fraction = binval / hint;
        if (std::abs(fraction) > 1e-5) {
          FCS_MSG_NOCLASS(logger,
                          "Warning: bin content is negative in histogram "
                              << hist->GetName() << " : " << hist->GetTitle()
                              << " binval=" << binval << " " << fraction * 100
                              << "% of integral=" << hist->Integral()
                              << ". Forcing bin to 0.");
        }
        binval = 0;
      }
      integral += binval;
      integral_vec[globalbin] = integral;
    }
  }

  const auto integral_at = [&integral_vec](int bin)
  { return integral_vec[static_cast<std::size_t>(bin)]; };
  for (first = 0; first < nbins; first++) {
    if (!FastCaloSim::exactly_equal(integral_at(first), 0.0)) {
      break;
    }
  }
  for (last = nbins - 1; last > 0; last--) {
    if (!FastCaloSim::exactly_equal(integral_at(last), integral)) {
      break;
    }
  }
  last++;

  if (integral <= 0) {
    FCS_MSG_NOCLASS(logger,
                    "Error: histogram "
                        << hist->GetName() << " : " << hist->GetTitle()
                        << " integral=" << integral << " is <=0");
  }

  return integral;
}

auto create_random_TH2(int nbinsx = 64, int nbinsy = 64) -> TH2*
{
  TH2* hist = new TH2F("test2D", "test2D", nbinsx, 0, 1, nbinsy, 0, 1);
  hist->Sumw2();
  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy <= nbinsy; ++iy) {
      hist->SetBinContent(ix,
                          iy,
                          (0.5 + gRandom->Rndm()) * (nbinsx + ix)
                              * (nbinsy * nbinsy / 2 + iy * iy));
      if (gRandom->Rndm() < 0.1) {
        hist->SetBinContent(ix, iy, 0);
      }
      hist->SetBinError(ix, iy, 0);
    }
  }
  return hist;
}
