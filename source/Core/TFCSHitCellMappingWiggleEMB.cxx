/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/Core/TFCSHitCellMappingWiggleEMB.h"

#include "CLHEP/Random/RandFlat.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "TMath.h"
#include "TVector2.h"

//=============================================
//======= TFCSHitCellMappingWiggleEMB =========
//=============================================

TFCSHitCellMappingWiggleEMB::TFCSHitCellMappingWiggleEMB(const char* name,
                                                         const char* title,
                                                         ICaloGeometry* geo)
    : TFCSHitCellMapping(name, title, geo)
{
  double wiggleLayer1[] = {
      0.0110626, 0.0242509, 0.0398173, 0.055761, 0.0736173, 0.0938847, 0.115154,
      0.13639,   0.157644,  0.178934,  0.200182, 0.221473,  0.242745,  0.264019,
      0.285264,  0.306527,  0.327811,  0.349119, 0.370387,  0.391668,  0.412922,
      0.434208,  0.45546,   0.476732,  0.498023, 0.51931,   0.540527,  0.561799,
      0.583079,  0.604358,  0.625614,  0.646864, 0.668112,  0.689351,  0.710629,
      0.731894,  0.75318,   0.774426,  0.795695, 0.81699,   0.838258,  0.859528,
      0.880783,  0.90202,   0.922515,  0.941276, 0.958477,  0.975062,  0.988922,
      1};
  double wiggleLayer2[] = {
      0.0127507, 0.0255775, 0.0395137, 0.0542644, 0.0695555, 0.0858206,
      0.102274,  0.119653,  0.137832,  0.156777,  0.176938,  0.197727,
      0.217576,  0.236615,  0.256605,  0.277766,  0.2995,    0.321951,
      0.344663,  0.367903,  0.392401,  0.417473,  0.443514,  0.470867,
      0.498296,  0.52573,   0.553114,  0.57921,   0.604326,  0.628822,
      0.652191,  0.674853,  0.697268,  0.718983,  0.739951,  0.759866,
      0.778877,  0.798762,  0.819559,  0.839789,  0.858923,  0.877327,
      0.894831,  0.911693,  0.92821,   0.94391,   0.959156,  0.973593,
      0.986752,  1};
  double wiggleLayer3[] = {
      0.0217932, 0.0438502, 0.0670992, 0.091085, 0.11651,  0.143038, 0.169524,
      0.196205,  0.222944,  0.249703,  0.276629, 0.303559, 0.33034,  0.356842,
      0.383579,  0.410385,  0.437272,  0.464214, 0.49118,  0.518202, 0.545454,
      0.572667,  0.600037,  0.627544,  0.655072, 0.6826,   0.709824, 0.733071,
      0.754764,  0.775672,  0.793834,  0.810904, 0.828219, 0.844119, 0.858339,
      0.871248,  0.882485,  0.894889,  0.907955, 0.920289, 0.931136, 0.941039,
      0.949844,  0.957641,  0.965787,  0.97392,  0.981706, 0.988892, 0.994527,
      1};

  for (int i = 0; i < 50; i++) {
    m_wiggleLayer1[i] = wiggleLayer1[i];
    m_wiggleLayer2[i] = wiggleLayer2[i];
    m_wiggleLayer3[i] = wiggleLayer3[i];
  }
}

double TFCSHitCellMappingWiggleEMB::doWiggle(double searchRand)
{
  int layer = calosample();

  double wiggle = 0.0;
  double cell_dphi = 0.0;

  // Define cell size in phi
  // TODO: Should come from geometry!
  if (layer == 1)
    cell_dphi = 0.0981748;
  if (layer == 2 || layer == 3)
    cell_dphi = 0.0245437;
  if (layer == 0 || layer > 3) {
    return 0.0;
  }

  // Now for layer dependent approach
  if (layer == 1) {
    int chosenBin = (Int_t)TMath::BinarySearch(50, m_wiggleLayer1, searchRand);
    double x_wigg = ((-0.98) + (chosenBin + 1) * 0.04) / 2;
    wiggle = x_wigg * cell_dphi / 4;
  }

  if (layer == 2) {
    int chosenBin = (Int_t)TMath::BinarySearch(50, m_wiggleLayer2, searchRand);
    double x_wigg = ((-0.98) + (chosenBin + 1) * 0.04) / 2;
    wiggle = x_wigg * cell_dphi;
  }

  if (layer == 3) {
    int chosenBin = (Int_t)TMath::BinarySearch(50, m_wiggleLayer3, searchRand);
    double x_wigg = ((-0.98) + (chosenBin + 1) * 0.04) / 2;
    wiggle = x_wigg * cell_dphi;
  }

  return wiggle;
}

FCSReturnCode TFCSHitCellMappingWiggleEMB::simulate_hit(
    Hit& hit,
    TFCSSimulationState& simulstate,
    const TFCSTruthState* truth,
    const TFCSExtrapolationState* extrapol)
{
  if (!simulstate.randomEngine()) {
    return FCSFatal;
  }

  int cs = calosample();

  double wiggle = 0.0;
  if (cs < 4 && cs > 0)
    wiggle = doWiggle(CLHEP::RandFlat::shoot(simulstate.randomEngine()));

  ATH_MSG_DEBUG("HIT: E=" << hit.E() << " cs=" << cs << " eta=" << hit.eta()
                          << " phi=" << hit.phi() << " wiggle=" << wiggle);

  double hit_phi_shifted = hit.phi() - wiggle;
  hit.phi() = TVector2::Phi_mpi_pi(hit_phi_shifted);

  return TFCSHitCellMapping::simulate_hit(hit, simulstate, truth, extrapol);
}
