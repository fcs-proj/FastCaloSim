#include <memory>

#include <gtest/gtest.h>

#include <CLHEP/Random/RanluxEngine.h>
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Geometry/CaloGeometryFromFile.h"
#include "TFile.h"

class BasicSimTests : public ::testing::Test
{
protected:
  static CaloGeometryFromFile* geo;
  static TFile* param_file;

  // Per-test-suite set-up.
  // Called before the first test in this test suite.
  // Note: this does not work if incurred via ctest
  static void SetUpTestSuite()
  {
    // Load geometry
    geo = new CaloGeometryFromFile();
    const std::string geo_tag = "ATLAS-R2-2016-01-00-01";
    const std::string geo_path = std::string(TEST_GEO_DIR) + geo_tag + ".root";
    const std::string hash_file_path =
        std::string(TEST_GEO_DIR) + "cellId_vs_cellHashId_map.txt";
    const std::string fcal1_geo_path =
        std::string(TEST_GEO_DIR) + "FCal1-electrodes.sorted.HV.09Nov2007.dat";
    const std::string fcal2_geo_path =
        std::string(TEST_GEO_DIR) + "FCal2-electrodes.sorted.HV.April2011.dat";
    const std::string fcal3_geo_path =
        std::string(TEST_GEO_DIR) + "FCal3-electrodes.sorted.HV.09Nov2007.dat";
    geo->LoadGeometryFromFile(geo_path, geo_tag, hash_file_path);
    geo->LoadFCalGeometryFromFiles(
        {fcal1_geo_path, fcal2_geo_path, fcal3_geo_path});

    // Load Parmaetrization file
    const std::string input_file_path = std::string(TEST_PARAM_DIR)
        + "TFCSParam_pid22_Mom16384_131072_eta_15_30.root";
    // param_file = std::unique_ptr<TFile> {TFile::Open(input_file_path.c_str(),
    // "READ")};
    param_file = new TFile(input_file_path.c_str(), "READ");
  }

  // Per-test-suite tear-down.
  // Called after the last test in this test suite.
  // Note: this does not work if incurred via ctest
  static void TearDownTestSuite()
  {
    param_file->Close();
    delete param_file;
    delete geo;
  }
};

// Define static members to make them accessible in the tests
CaloGeometryFromFile* BasicSimTests::geo = nullptr;
TFile* BasicSimTests::param_file = nullptr;

TEST_F(BasicSimTests, ReadParamFile)
{
  Long64_t bytes_read = param_file->GetBytesRead();
  EXPECT_EQ(bytes_read, 9213);
}

TEST_F(BasicSimTests, CheckParamObject)
{
  const std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      param_file->Get(paramsObject.c_str()));
  EXPECT_NE(param, nullptr);
}

TEST_F(BasicSimTests, GetGeoCellBySamplingEtaPhi)
{
  const CaloDetDescrElement* cell = geo->getDDE(2, 0.24, 0.24);
  long int cell_id = cell->identify();
  EXPECT_EQ(cell_id, 3296654795753914368);
  EXPECT_EQ(cell->getSampling(), 2);
  EXPECT_NEAR(cell->eta(), 0.238859, 1e-4);
  EXPECT_NEAR(cell->phi(), 0.228048, 1e-4);
}

TEST_F(BasicSimTests, GetGeoCellByIdentifier)
{
  const CaloDetDescrElement* cell = geo->getDDE(3260641881524011008);
  long int cell_id = cell->identify();
  EXPECT_EQ(cell_id, 3260641881524011008);
  EXPECT_EQ(cell->getSampling(), 2);
  EXPECT_NEAR(cell->eta(), -0.411563, 1e-4);
  EXPECT_NEAR(cell->phi(), 1.63532, 1e-4);
}

TEST_F(BasicSimTests, DoDummySimulation)
{
  std::string paramsObject {"SelPDGID"};
  TFCSParametrizationBase* param = static_cast<TFCSParametrizationBase*>(
      param_file->Get(paramsObject.c_str()));

  // Set up the geometry
  param->set_geometry(geo);
  // Set logging level
  param->setLevel(MSG::Level::VERBOSE);

  // Set up a random engine
  std::unique_ptr<CLHEP::RanluxEngine> rnd_engine =
      std::make_unique<CLHEP::RanluxEngine>();
  rnd_engine->setSeed(42);

  // Set up simulation state and assign random engine
  TFCSSimulationState simul_state;
  simul_state.setRandomEngine(rnd_engine.get());

  // Set up particle to simulate
  TFCSTruthState truth_state;
  truth_state.SetPtEtaPhiM(65536, 0.225, 1.8, 0);
  truth_state.set_pdgid(22);

  // Set up a dummy extrapolation state
  TFCSExtrapolationState extrapol_state;
  extrapol_state.set_IDCaloBoundary_eta(truth_state.Eta());
  for (int i = 0; i < 24; ++i) {
    extrapol_state.set_eta(
        i, TFCSExtrapolationState::SUBPOS_ENT, truth_state.Eta());
    extrapol_state.set_eta(
        i, TFCSExtrapolationState::SUBPOS_EXT, truth_state.Eta());
    extrapol_state.set_eta(
        i, TFCSExtrapolationState::SUBPOS_MID, truth_state.Eta());
    extrapol_state.set_phi(
        i, TFCSExtrapolationState::SUBPOS_ENT, truth_state.Phi());
    extrapol_state.set_phi(
        i, TFCSExtrapolationState::SUBPOS_EXT, truth_state.Phi());
    extrapol_state.set_phi(
        i, TFCSExtrapolationState::SUBPOS_MID, truth_state.Phi());
    extrapol_state.set_r(i, TFCSExtrapolationState::SUBPOS_ENT, 1500 + i * 10);
    extrapol_state.set_r(i, TFCSExtrapolationState::SUBPOS_EXT, 1510 + i * 10);
    extrapol_state.set_r(i, TFCSExtrapolationState::SUBPOS_MID, 1505 + i * 10);
    extrapol_state.set_z(i, TFCSExtrapolationState::SUBPOS_ENT, 3500 + i * 10);
    extrapol_state.set_z(i, TFCSExtrapolationState::SUBPOS_EXT, 3510 + i * 10);
    extrapol_state.set_z(i, TFCSExtrapolationState::SUBPOS_MID, 3505 + i * 10);
  }

  // Perform the simulation call
  param->simulate(simul_state, &truth_state, &extrapol_state);

  // Check that simulation was successful
  EXPECT_NEAR(simul_state.E(), 66943.7, 1e-1);
  EXPECT_NEAR(simul_state.E(0), 1096.21, 1e-1);
  EXPECT_NEAR(simul_state.E(1), 19513.9, 1e-1);
  EXPECT_NEAR(simul_state.E(2), 46197.6, 1e-1);
  EXPECT_NEAR(simul_state.E(3), 85.8894, 1e-1);
  EXPECT_NEAR(simul_state.E(12), 50.1765, 1e-1);
}