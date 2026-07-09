// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef ISF_FASTCALOSIMEVENT_TFCSEnergyAndHitGANV2_h
#define ISF_FASTCALOSIMEVENT_TFCSEnergyAndHitGANV2_h

#include <mutex>
#include <string>

#include "FastCaloSim/Core/TFCSGANEtaSlice.h"
#include "FastCaloSim/Core/TFCSGANXMLParameters.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

class CaloGeo;

// forward declare lwtnn dependencies
namespace lwt
{
class LightweightGraph;
}

class TFCSEnergyAndHitGANV2 : public TFCSParametrizationBinnedChain
{
public:
  TFCSEnergyAndHitGANV2(const char* name = nullptr,
                        const char* title = nullptr,
                        CaloGeo* geo = nullptr);

  void set_geometry(CaloGeo* geo) override
  {
    m_geo = geo;
    TFCSParametrizationBase::set_geometry(geo);
  };

  ~TFCSEnergyAndHitGANV2() override;

  auto is_match_Ekin_bin(int /*Ekin_bin*/) const -> bool override
  {
    return true;
  };
  auto is_match_calosample(int calosample) const -> bool override;
  auto is_match_all_Ekin_bin() const -> bool override { return true; };
  auto is_match_all_calosample() const -> bool override { return false; };

  /// Status bit for chain persistency
  enum FCSGANfreemem
  {
    kGANfreemem =
        BIT(17)  ///< Set this bit in the TObject bit if the memory for
                 ///< m_input should be freed after reading in athena
  };

  auto GANfreemem() const -> bool { return TestBit(kGANfreemem); };
  void set_GANfreemem() { SetBit(kGANfreemem); };
  void reset_GANfreemem() { ResetBit(kGANfreemem); };

  /// Status bit for energy initialization
  enum FCSEnergyInitializationStatusBits
  {
    kOnlyScaleEnergy =
        BIT(18)  ///< Set this bit in the TObject bit field the simulated energy
                 ///< should only be scaled by the GAN
  };

  auto OnlyScaleEnergy() const -> bool { return TestBit(kOnlyScaleEnergy); };
  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); };
  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); };

  /// use the layer to be done as binning of the GAN chain
  auto get_bin(TFCSSimulationState& simulstate,
               const TFCSTruthState*,
               const TFCSExtrapolationState*) const -> int override
  {
    return simulstate.getAuxInfo<int>("GANlayer"_FCShash);
  };
  auto get_variable_text(TFCSSimulationState& simulstate,
                         const TFCSTruthState*,
                         const TFCSExtrapolationState*) const
      -> const std::string override;

  auto get_nr_of_init(unsigned int bin) const -> unsigned int;
  void set_nr_of_init(unsigned int bin, unsigned int ninit);

  auto get_Binning() const -> const TFCSGANXMLParameters::Binning&
  {
    return m_param.GetBinning();
  };
  auto get_ExtrapolationWeights() const
      -> const TFCSGANEtaSlice::ExtrapolatorWeights&
  {
    return m_slice->GetExtrapolatorWeights();
  };

  auto initializeNetwork(int const& pid,
                         int const& etaMin,
                         const std::string& FastCaloGANInputFolderName) -> bool;

  auto fillEnergy(TFCSSimulationState& simulstate,
                  const TFCSTruthState* truth,
                  const TFCSExtrapolationState* extrapol) const -> bool;
  auto simulate(TFCSSimulationState& simulstate,
                const TFCSTruthState* truth,
                const TFCSExtrapolationState* extrapol) const
      -> FCSReturnCode override;

  void Print(Option_t* option = "") const override;

protected:
  CaloGeo* m_geo;  //! do not persistify

  void SetRegionAndSliceFromXML(int pid,
                                int etaMax,
                                std::string FastCaloGANInputFolderName);

private:
  static auto GetBinsInFours(double const bins) -> int;
  auto GetAlphaBinsForRBin(const TAxis* x, int ix, int yBinNum) const -> int;

  std::vector<int> m_bin_ninit;

  // Persistify configuration in string m_input. A custom Streamer(...) builds
  // m_graph on the fly when reading from file. Inside Athena, if GANfreemem()
  // is true, the content of m_input is deleted after reading in order to free
  // memory

  TFCSGANEtaSlice* m_slice = nullptr;
  TFCSGANXMLParameters m_param;
  mutable std::mutex m_mutex;  //! Do not persistify

  ClassDefOverride(TFCSEnergyAndHitGANV2, 2)  // TFCSEnergyAndHitGANV2
};

#endif
