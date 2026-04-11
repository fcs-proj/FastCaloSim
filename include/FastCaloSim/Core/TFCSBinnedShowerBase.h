// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#ifndef TFCSBinnedShowerBase_h
#define TFCSBinnedShowerBase_h

// local includes
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

// External includes
#include <fstream>
#include <tuple>
#include <vector>

#include <RtypesCore.h>
#include <TMath.h>

class ICaloGeometry;

class TFCSBinnedShowerBase : public TFCSLateralShapeParametrizationHitBase
{
public:
  TFCSBinnedShowerBase(const char* name = nullptr, const char* title = nullptr);

  virtual ~TFCSBinnedShowerBase();

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

  // Fill layer energies
  virtual FCSReturnCode simulate(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) const override;

  // Do hit simulation
  virtual FCSReturnCode simulate_hit(
      Hit& hit,
      TFCSSimulationState& simulstate,
      const TFCSTruthState* truth,
      const TFCSExtrapolationState* extrapol) override;

  ICaloGeometry* get_geometry() { return m_geo; };

  virtual void set_geometry(ICaloGeometry* geo) override
  {
    m_geo = geo;
    TFCSParametrizationBase::set_geometry(geo);
  };

  virtual int get_number_of_hits(
      TFCSSimulationState& simulstate,
      const TFCSTruthState* /*truth*/,
      const TFCSExtrapolationState* /*extrapol*/) const override
  {
    long unsigned int nhits = get_n_hits(simulstate, calosample());
    if (nhits == 0)
      return -1;
    return static_cast<int>(nhits);
  }

protected:
  ICaloGeometry* m_geo;  //! do not persistify

  // Called at the beginning of the simulation to store and or generate the
  // needed shower data for the current event
  virtual void get_event(TFCSSimulationState& simulstate,
                         float eta_center,
                         float phi_center,
                         float e_init,
                         long unsigned int reference_layer_index) const = 0;

  // Returns the number of bins that are used in the given layer
  virtual long unsigned int get_n_hits(TFCSSimulationState& simulstate,
                                       long unsigned int layer_index) const = 0;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  virtual float get_layer_energy(TFCSSimulationState& simulstate,
                                 long unsigned int layer_index) const = 0;

  // Returns the position and energy of the corresponding hit in the given
  // event, layer and bin
  virtual std::tuple<float, float, float> get_hit_position_and_energy(
      TFCSSimulationState& simulstate,
      long unsigned int layer_index,
      long unsigned int hit_index) const = 0;

  // Delete all pointers that were created in get_event()
  virtual void delete_event(TFCSSimulationState& simulstate) const = 0;

private:
  const long unsigned int m_n_layers =
      CaloCell_ID_FCS::CaloSample_FCS::MaxSample;
  ClassDefOverride(TFCSBinnedShowerBase, 1)  // TFCSBinnedShowerBase
};

#endif
