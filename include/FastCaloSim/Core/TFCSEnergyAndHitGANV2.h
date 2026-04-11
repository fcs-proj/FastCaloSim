// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#ifndef TFCSEnergyAndHitGANV2_h
#define TFCSEnergyAndHitGANV2_h

#include <string>

#include "FastCaloSim/Core/TFCSGANEtaSlice.h"
#include "FastCaloSim/Core/TFCSGANXMLParameters.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

namespace FastCaloSim::Geometry
{
class CaloGeo;
}
namespace FastCaloSim::Core
{
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

  virtual void set_geometry(CaloGeo* geo) override { m_geo = geo; };

  virtual ~TFCSEnergyAndHitGANV2();

  virtual bool is_match_Ekin_bin(int /*Ekin_bin*/) const override
  {
    return true;
  };
  virtual bool is_match_calosample(int calosample) const override;
  virtual bool is_match_all_Ekin_bin() const override { return true; };
  virtual bool is_match_all_calosample() const override { return false; };

  /// Status bit for chain persistency
  enum FCSGANfreemem
  {
    kGANfreemem =
        BIT(17)  ///< Set this bit in the TObject bit if the memory for
                 ///< m_input should be freed after reading in athena
  };

  bool GANfreemem() const { return TestBit(kGANfreemem); };
  void set_GANfreemem() { SetBit(kGANfreemem); };
  void reset_GANfreemem() { ResetBit(kGANfreemem); };

  /// Status bit for energy initialization
  enum FCSEnergyInitializationStatusBits
  {
    kOnlyScaleEnergy =
        BIT(18)  ///< Set this bit in the TObject bit field the simulated energy
                 ///< should only be scaled by the GAN
  };

  bool OnlyScaleEnergy() const { return TestBit(kOnlyScaleEnergy); };
  void set_OnlyScaleEnergy() { SetBit(kOnlyScaleEnergy); };
  void reset_OnlyScaleEnergy() { ResetBit(kOnlyScaleEnergy); };

  /// use the layer to be done as binning of the GAN chain
  virtual int get_bin(TFCSSimulationState& simulstate,
                      const TFCSTruthState*,
                      const TFCSExtrapolationState*) const override
  {
    return simulstate.getAuxInfo<int>("GANlayer"_FCShash);
  };
  virtual const std::string get_variable_text(
      TFCSSimulationState& simulstate,
      const TFCSTruthState*,
      const TFCSExtrapolationState*) const override;

  unsigned int get_nr_of_init(unsigned int bin) const;
  void set_nr_of_init(unsigned int bin, unsigned int ninit);

  const TFCSGANXMLParameters::Binning& get_Binning() const
  {
    return m_param.GetBinning();
  };
  const TFCSGANEtaSlice::ExtrapolatorWeights get_ExtrapolationWeights() const
  {
    return m_slice->GetExtrapolatorWeights();
  };

  bool initializeNetwork(int const& pid,
                         int const& etaMin,
                         const std::string& FastCaloGANInputFolderName);

  bool fillEnergy(TFCSSimulationState& simulstate,
                  const TFCSTruthState* truth,
                  const TFCSExtrapolationState* extrapol) const;
  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  virtual void Print(Option_t* option = "") const override;

protected:
  CaloGeo* m_geo;  //! do not persistify

  void SetRegionAndSliceFromXML(int pid,
                                int etaMax,
                                std::string FastCaloGANInputFolderName);

private:
  static int GetBinsInFours(double const& bins);
  int GetAlphaBinsForRBin(const TAxis* x, int ix, int yBinNum) const;

  std::vector<int> m_bin_ninit;

  // Persistify configuration in string m_input. A custom Streamer(...) builds
  // m_graph on the fly when reading from file. Inside Athena, if GANfreemem()
  // is true, the content of m_input is deleted after reading in order to free
  // memory

  TFCSGANEtaSlice* m_slice = nullptr;
  TFCSGANXMLParameters m_param;

  ClassDefOverride(TFCSEnergyAndHitGANV2, 2)  // TFCSEnergyAndHitGANV2
};

}  // namespace FastCaloSim::Core

#endif
