// Copyright (c) 2025 CERN for the benefit of the FastCaloSim project
#ifndef FastCaloSim_LegacyCompatibility_h
#define FastCaloSim_LegacyCompatibility_h

#include <TBufferFile.h>
#include <TInterpreter.h>
#include <TROOT.h>

namespace FastCaloSim
{

class LegacyCompatibility
{
public:
  /// Call this function before opening any legacy ROOT files
  /// It sets up all the necessary typedefs for backward compatibility
  static void setupLegacyClassMapping()
  {
    static bool initialized = false;
    if (initialized)
      return;

    // Create all the necessary typedefs in ROOT's interpreter
    gInterpreter->Declare(R"(
            typedef FastCaloSim::Core::DoubleArray DoubleArray;
            typedef FastCaloSim::Core::IntArray IntArray;
            typedef FastCaloSim::Core::MLogging MLogging;
            typedef FastCaloSim::Core::TFCS1DFunction TFCS1DFunction;
            typedef FastCaloSim::Core::TFCS1DFunctionHistogram TFCS1DFunctionHistogram;
            typedef FastCaloSim::Core::TFCS1DFunctionInt16Histogram TFCS1DFunctionInt16Histogram;
            typedef FastCaloSim::Core::TFCS1DFunctionInt32Histogram TFCS1DFunctionInt32Histogram;
            typedef FastCaloSim::Core::TFCS1DFunctionRegression TFCS1DFunctionRegression;
            typedef FastCaloSim::Core::TFCS1DFunctionRegressionTF TFCS1DFunctionRegressionTF;
            typedef FastCaloSim::Core::TFCS1DFunctionSpline TFCS1DFunctionSpline;
            typedef FastCaloSim::Core::TFCS2DFunction TFCS2DFunction;
            typedef FastCaloSim::Core::TFCS2DFunctionHistogram TFCS2DFunctionHistogram;
            typedef FastCaloSim::Core::TFCS2DFunctionLateralShapeParametrization TFCS2DFunctionLateralShapeParametrization;
            typedef FastCaloSim::Core::TFCSCenterPositionCalculation TFCSCenterPositionCalculation;
            typedef FastCaloSim::Core::TFCSEnergyAndHitGANV2 TFCSEnergyAndHitGANV2;
            typedef FastCaloSim::Core::TFCSEnergyBinParametrization TFCSEnergyBinParametrization;
            typedef FastCaloSim::Core::TFCSEnergyInterpolationHistogram TFCSEnergyInterpolationHistogram;
            typedef FastCaloSim::Core::TFCSEnergyInterpolationLinear TFCSEnergyInterpolationLinear;
            typedef FastCaloSim::Core::TFCSEnergyInterpolationPiecewiseLinear TFCSEnergyInterpolationPiecewiseLinear;
            typedef FastCaloSim::Core::TFCSEnergyInterpolationSpline TFCSEnergyInterpolationSpline;
            typedef FastCaloSim::Core::TFCSEnergyParametrization TFCSEnergyParametrization;
            typedef FastCaloSim::Core::TFCSEnergyRenormalization TFCSEnergyRenormalization;
            typedef FastCaloSim::Core::TFCSExtrapolationState TFCSExtrapolationState;
            typedef FastCaloSim::Core::TFCSFlatLateralShapeParametrization TFCSFlatLateralShapeParametrization;
            typedef FastCaloSim::Core::TFCSFunction TFCSFunction;
            typedef FastCaloSim::Core::TFCSGANEtaSlice TFCSGANEtaSlice;
            typedef FastCaloSim::Core::TFCSGANLWTNNHandler TFCSGANLWTNNHandler;
            typedef FastCaloSim::Core::TFCSGANXMLParameters TFCSGANXMLParameters;
            typedef FastCaloSim::Core::TFCSHistoLateralShapeGausLogWeight TFCSHistoLateralShapeGausLogWeight;
            typedef FastCaloSim::Core::TFCSHistoLateralShapeGausLogWeightHitAndMiss TFCSHistoLateralShapeGausLogWeightHitAndMiss;
            typedef FastCaloSim::Core::TFCSHistoLateralShapeParametrization TFCSHistoLateralShapeParametrization;
            typedef FastCaloSim::Core::TFCSHistoLateralShapeParametrizationFCal TFCSHistoLateralShapeParametrizationFCal;
            typedef FastCaloSim::Core::TFCSHistoLateralShapeWeight TFCSHistoLateralShapeWeight;
            typedef FastCaloSim::Core::TFCSHistoLateralShapeWeightHitAndMiss TFCSHistoLateralShapeWeightHitAndMiss;
            typedef FastCaloSim::Core::TFCSHitCellMapping TFCSHitCellMapping;
            typedef FastCaloSim::Core::TFCSHitCellMappingFCal TFCSHitCellMappingFCal;
            typedef FastCaloSim::Core::TFCSHitCellMappingWiggle TFCSHitCellMappingWiggle;
            typedef FastCaloSim::Core::TFCSHitCellMappingWiggleEMB TFCSHitCellMappingWiggleEMB;
            typedef FastCaloSim::Core::TFCSInitWithEkin TFCSInitWithEkin;
            typedef FastCaloSim::Core::TFCSInvisibleParametrization TFCSInvisibleParametrization;
            typedef FastCaloSim::Core::TFCSLateralShapeParametrization TFCSLateralShapeParametrization;
            typedef FastCaloSim::Core::TFCSLateralShapeParametrizationFluctChain TFCSLateralShapeParametrizationFluctChain;
            typedef FastCaloSim::Core::TFCSLateralShapeParametrizationHitBase TFCSLateralShapeParametrizationHitBase;
            typedef FastCaloSim::Core::TFCSLateralShapeParametrizationHitChain TFCSLateralShapeParametrizationHitChain;
            typedef FastCaloSim::Core::TFCSLateralShapeParametrizationHitNumberFromE TFCSLateralShapeParametrizationHitNumberFromE;
            typedef FastCaloSim::Core::TFCSLateralShapeTuning TFCSLateralShapeTuning;
            typedef FastCaloSim::Core::TFCSONNXHandler TFCSONNXHandler;
            typedef FastCaloSim::Core::TFCSPCAEnergyParametrization TFCSPCAEnergyParametrization;
            typedef FastCaloSim::Core::TFCSParametrization TFCSParametrization;
            typedef FastCaloSim::Core::TFCSParametrizationAbsEtaSelectChain TFCSParametrizationAbsEtaSelectChain;
            typedef FastCaloSim::Core::TFCSParametrizationBase TFCSParametrizationBase;
            typedef FastCaloSim::Core::TFCSParametrizationBinnedChain TFCSParametrizationBinnedChain;
            typedef FastCaloSim::Core::TFCSParametrizationChain TFCSParametrizationChain;
            typedef FastCaloSim::Core::TFCSParametrizationEbinChain TFCSParametrizationEbinChain;
            typedef FastCaloSim::Core::TFCSParametrizationEkinSelectChain TFCSParametrizationEkinSelectChain;
            typedef FastCaloSim::Core::TFCSParametrizationEtaSelectChain TFCSParametrizationEtaSelectChain;
            typedef FastCaloSim::Core::TFCSParametrizationFloatSelectChain TFCSParametrizationFloatSelectChain;
            typedef FastCaloSim::Core::TFCSParametrizationPDGIDSelectChain TFCSParametrizationPDGIDSelectChain;
            typedef FastCaloSim::Core::TFCSParametrizationPlaceholder TFCSParametrizationPlaceholder;
            typedef FastCaloSim::Core::TFCSPredictExtrapWeights TFCSPredictExtrapWeights;
            typedef FastCaloSim::Core::TFCSSimpleLWTNNHandler TFCSSimpleLWTNNHandler;
            typedef FastCaloSim::Core::TFCSSimulationState TFCSSimulationState;
            typedef FastCaloSim::Core::TFCSTruthState TFCSTruthState;
            typedef FastCaloSim::Core::VNetworkBase VNetworkBase;
            typedef FastCaloSim::Core::VNetworkLWTNN VNetworkLWTNN;

            namespace ISF_FCS {
                typedef FastCaloSim::Core::MLogging MLogging;
            }

        )");

    initialized = true;
  }
};

}  // namespace FastCaloSim

#endif  // FastCaloSim_LegacyCompatibility_h
