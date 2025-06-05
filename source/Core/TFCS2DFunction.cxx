// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

#include <iostream>

#include "FastCaloSim/Core/TFCS2DFunction.h"

#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "TH2.h"
#include "TRandom.h"

//=============================================
//======= TFCS2DFunction =========
//=============================================
using namespace FastCaloSim::Core;

void TFCS2DFunction::rnd_to_fct(float value[], const float rnd[]) const
{
  rnd_to_fct(value[0], value[1], rnd[0], rnd[1]);
}

//================================================================================================================================

double TFCS2DFunction::CheckAndIntegrate2DHistogram(
    const TH2* hist, std::vector<double>& integral_vec, int& first, int& last)
{
  MLogging logger;
  Int_t nbinsx = hist->GetNbinsX();
  Int_t nbinsy = hist->GetNbinsY();
  Int_t nbins = nbinsx * nbinsy;

  float integral = 0.0;
  float hint = hist->Integral();
  integral_vec.resize(nbins);

  for (int ix = 1; ix <= nbinsx; ix++) {
    for (int iy = 1; iy <= nbinsy; iy++) {
      int globalbin = (ix - 1) * nbinsy + iy - 1;
      float binval = hist->GetBinContent(ix, iy);
      if (binval < 0) {
        // Can't work if a bin is negative, forcing bins to 0 in this case
        double fraction = binval / hint;
        if (TMath::Abs(fraction) > 1e-5) {
          MSG_NOCLASS(logger,
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

  for (first = 0; first < nbins; first++)
    if (integral_vec[first] != 0)
      break;
  for (last = nbins - 1; last > 0; last--)
    if (integral_vec[last] != integral)
      break;
  last++;

  if (integral <= 0) {
    MSG_NOCLASS(logger,
                "Error: histogram " << hist->GetName() << " : "
                                    << hist->GetTitle()
                                    << " integral=" << integral << " is <=0");
  }

  return integral;
}

TH2* create_random_TH2(int nbinsx = 64, int nbinsy = 64)
{
  TH2* hist = new TH2F("test2D", "test2D", nbinsx, 0, 1, nbinsy, 0, 1);
  hist->Sumw2();
  for (int ix = 1; ix <= nbinsx; ++ix) {
    for (int iy = 1; iy <= nbinsy; ++iy) {
      hist->SetBinContent(ix,
                          iy,
                          (0.5 + gRandom->Rndm()) * (nbinsx + ix)
                              * (nbinsy * nbinsy / 2 + iy * iy));
      if (gRandom->Rndm() < 0.1)
        hist->SetBinContent(ix, iy, 0);
      hist->SetBinError(ix, iy, 0);
    }
  }
  return hist;
}

#pragma GCC diagnostic pop
