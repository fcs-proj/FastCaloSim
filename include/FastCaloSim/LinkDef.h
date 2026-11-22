// clang-format off
// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/DoubleArray.h"
#include "FastCaloSim/Core/IntArray.h"
#include "FastCaloSim/Core/TFCSFunction.h"
#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "FastCaloSim/Core/TFCS1DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionInt32Histogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionRegression.h"
#include "FastCaloSim/Core/TFCS1DFunctionRegressionTF.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateHelpers.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateHistogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateInterpolationHistogram.h"
#include "FastCaloSim/Core/TFCS2DFunction.h"
#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"

#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSParametrizationPlaceholder.h"
#include "FastCaloSim/Core/TFCSParametrization.h"
#include "FastCaloSim/Core/TFCSInvisibleParametrization.h"
#include "FastCaloSim/Core/TFCSInitWithEkin.h"
#include "FastCaloSim/Core/TFCSEnergyInterpolationHistogram.h"
#include "FastCaloSim/Core/TFCSEnergyInterpolationLinear.h"
#include "FastCaloSim/Core/TFCSEnergyInterpolationPiecewiseLinear.h"
#include "FastCaloSim/Core/TFCSEnergyInterpolationSpline.h"
#include "FastCaloSim/Core/TFCSParametrizationChain.h"
#include "FastCaloSim/Core/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/Core/TFCSParametrizationFloatSelectChain.h"
#include "FastCaloSim/Core/TFCSParametrizationPDGIDSelectChain.h"
#include "FastCaloSim/Core/TFCSParametrizationEbinChain.h"
#include "FastCaloSim/Core/TFCSParametrizationEkinSelectChain.h"
#include "FastCaloSim/Core/TFCSParametrizationEtaSelectChain.h"
#include "FastCaloSim/Core/TFCSParametrizationAbsEtaSelectChain.h"

#include "FastCaloSim/Core/TFCSEnergyParametrization.h"
#include "FastCaloSim/Core/TFCSPCAEnergyParametrization.h"
#include "FastCaloSim/Core/TFCSEnergyBinParametrization.h"
#include "FastCaloSim/Core/TFCSPhiModulationCorrection.h"
#include "FastCaloSim/Core/TFCSBinnedShowerBase.h"
#include "FastCaloSim/Core/TFCSBinnedShower.h"
#include "FastCaloSim/Core/TFCSBinnedShowerONNX.h"
#include "FastCaloSim/Core/TFCSMLCalorimeterSimulator.h"

#include "FastCaloSim/Core/TFCSGANXMLParameters.h"
#include "FastCaloSim/Core/TFCSGANLWTNNHandler.h"
#include "FastCaloSim/Core/TFCSGANEtaSlice.h"
#include "FastCaloSim/Core/TFCSEnergyAndHitGANV2.h"
#include "FastCaloSim/Core/TFCSPredictExtrapWeights.h"
#include "FastCaloSim/Core/VNetworkBase.h"
#include "FastCaloSim/Core/VNetworkLWTNN.h"
#include "FastCaloSim/Core/TFCSSimpleLWTNNHandler.h"
#include "FastCaloSim/Core/TFCSONNXHandler.h"
#include "FastCaloSim/Core/TFCSNetworkFactory.h"

#include "FastCaloSim/Core/TFCSLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeGausLogWeight.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeGausLogWeightHitAndMiss.h"

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitChain.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFixedHitChain.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationFluctChain.h"
#include "FastCaloSim/Core/TFCSCenterPositionCalculation.h"
#include "FastCaloSim/Core/TFCSFlatLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeParametrizationFCal.h"
#include "FastCaloSim/Core/TFCS2DFunctionLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeWeight.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeWeightHitAndMiss.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitNumberFromE.h"
#include "FastCaloSim/Core/TFCSLateralShapeTuning.h"
#include "FastCaloSim/Core/TFCSHitCellMapping.h"
#include "FastCaloSim/Core/TFCSHitCellMappingFCal.h"
#include "FastCaloSim/Core/TFCSHitCellMappingWiggle.h"
#include "FastCaloSim/Core/TFCSHitCellMappingWiggleEMB.h"
#include "FastCaloSim/Core/TFCSEnergyRenormalization.h"
#include "FastCaloSim/Core/TFCSVoxelHistoLateralCovarianceFluctuations.h"

#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

#ifdef __CINT__
#pragma link C++ class DoubleArray + ;
#pragma link C++ class IntArray + ;
#pragma link C++ class TFCSFunction + ;
#pragma link C++ class TFCS1DFunction + ;
#pragma link C++ class TFCS1DFunctionHistogram + ;
#pragma link C++ class TFCS1DFunctionInt32Histogram + ;
// Build the transient search table after reading the contents
#pragma read sourceClass="TFCS1DFunctionInt32Histogram" targetClass="TFCS1DFunctionInt32Histogram" version="[1-]" source="" target="m_guide" code="{ m_guide.build(newObj->get_HistoContents()); }"
#pragma link C++ class TFCS1DFunctionRegression + ;
#pragma link C++ class TFCS1DFunctionRegressionTF + ;

/// Linkdefs needed for common classes for template based histograms
#pragma link C++ class TFCS1DFunction_Numeric < uint8_t, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < uint16_t, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < uint32_t, float> + ;

#pragma link C++ class TFCS1DFunction_Array < uint8_t> - ;
#pragma link C++ class TFCS1DFunction_Array < uint16_t> - ;
#pragma link C++ class TFCS1DFunction_Array < uint32_t> - ;

#pragma link C++ class TFCS1DFunction_HistogramContent < uint8_t, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < uint16_t, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < uint32_t, float> + ;

#pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float, uint8_t, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float, uint32_t, float> + ;

#pragma link C++ class TFCS1DFunction_HistogramInt8BinEdges + ;
#pragma link C++ class TFCS1DFunction_HistogramInt32BinEdges + ;

/// Linkdefs needed for 1D template based histograms
#pragma link C++ class TFCS1DFunctionTemplateHistogram < TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram < TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram < TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> + ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram < TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> + ;

#pragma link C++ class TFCS1DFunctionInt8Int8Histogram + ;
#pragma link C++ class TFCS1DFunctionInt8Int16Histogram + ;
#pragma link C++ class TFCS1DFunctionInt8Int32Histogram + ;
#pragma link C++ class TFCS1DFunctionInt32Int32Histogram + ;

#pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram < TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> + ;
#pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram < TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> + ;

#pragma link C++ class TFCS1DFunctionInt8Int8InterpolationHistogram + ;
#pragma link C++ class TFCS1DFunctionInt8Int16InterpolationHistogram + ;

/// End Linkdefs needed for template based histograms

#pragma link C++ class TFCS2DFunction + ;
#pragma link C++ class TFCS2DFunctionHistogram + ;
// Build the transient search table after reading the contents
#pragma read sourceClass="TFCS2DFunctionHistogram" targetClass="TFCS2DFunctionHistogram" version="[1-]" source="" target="m_guide" code="{ m_guide.build(newObj->get_HistoContents()); }"
#pragma link C++ class FastCaloSim::MLogging + ;
// MLogging was called ISF_FCS::MLogging when the library still lived in
// Athena. Old parametrization files reference that name in their streamer
// info; this rule maps it onto the renamed class so they read cleanly.
// The version selector must include 0: MLogging is persistified with class
// version 0 (ClassDef(MLogging, 0)), and rootcling defaults a missing
// version attribute to "[1-]", which would never match the legacy streamer
// info.
#pragma read sourceClass="ISF_FCS::MLogging" version="[0-]" targetClass="FastCaloSim::MLogging";
#pragma link C++ class TFCSParametrizationBase + ;
#pragma link C++ class TFCSParametrizationPlaceholder + ;
#pragma link C++ class TFCSParametrization + ;
#pragma link C++ class TFCSInvisibleParametrization + ;
#pragma link C++ class TFCSInitWithEkin + ;
#pragma link C++ class TFCSEnergyInterpolationHistogram + ;
#pragma link C++ class TFCSEnergyInterpolationLinear + ;
#pragma link C++ class TFCSEnergyInterpolationPiecewiseLinear - ;
#pragma link C++ class TFCSEnergyInterpolationSpline + ;
#pragma link C++ class TFCSParametrizationChain - ;
#pragma link C++ class TFCSParametrizationBinnedChain + ;
#pragma link C++ class TFCSParametrizationFloatSelectChain + ;
#pragma link C++ class TFCSParametrizationPDGIDSelectChain + ;
#pragma link C++ class TFCSParametrizationEbinChain + ;
#pragma link C++ class TFCSParametrizationEkinSelectChain + ;
#pragma link C++ class TFCSParametrizationEtaSelectChain + ;
#pragma link C++ class TFCSParametrizationAbsEtaSelectChain + ;

#pragma link C++ class TFCSEnergyParametrization + ;
#pragma link C++ class TFCSPCAEnergyParametrization - ;
#pragma link C++ class TFCSEnergyBinParametrization + ;
#pragma link C++ class TFCSPhiModulationCorrection + ;
#pragma link C++ class TFCSBinnedShowerBase + ;
#pragma link C++ class TFCSBinnedShower - ;
#pragma link C++ class TFCSBinnedShowerONNX - ;
#pragma link C++ class TFCSMLCalorimeterSimulator + ;
#pragma link C++ class TFCSGANXMLParameters + ;
#pragma link C++ class TFCSGANLWTNNHandler - ;
#pragma link C++ class TFCSGANEtaSlice + ;
#pragma link C++ class TFCSEnergyAndHitGANV2 + ;
#pragma link C++ class TFCSPredictExtrapWeights - ;

#pragma link C++ class VNetworkBase + ;
#pragma link C++ class VNetworkLWTNN + ;
#pragma link C++ class TFCSSimpleLWTNNHandler - ;
#pragma link C++ class TFCSONNXHandler - ;

#pragma link C++ class TFCSLateralShapeParametrization + ;
#pragma link C++ class TFCSLateralShapeParametrizationHitBase + ;
#pragma link C++ class TFCSLateralShapeParametrizationHitChain + ;
#pragma link C++ class TFCSLateralShapeParametrizationFixedHitChain + ;
#pragma link C++ class TFCSLateralShapeParametrizationFluctChain + ;
#pragma link C++ class TFCSCenterPositionCalculation + ;
#pragma link C++ class TFCSFlatLateralShapeParametrization + ;
#pragma link C++ class TFCSHistoLateralShapeParametrization + ;
#pragma link C++ class TFCSHistoLateralShapeParametrizationFCal + ;
#pragma link C++ class TFCS2DFunctionLateralShapeParametrization + ;
#pragma link C++ class TFCSHistoLateralShapeWeight + ;
#pragma link C++ class TFCSHistoLateralShapeGausLogWeight + ;
#pragma link C++ class TFCSHistoLateralShapeWeightHitAndMiss + ;
#pragma link C++ class TFCSHistoLateralShapeGausLogWeightHitAndMiss + ;
#pragma link C++ class TFCSLateralShapeParametrizationHitNumberFromE + ;
#pragma link C++ class TFCSLateralShapeTuning + ;
#pragma link C++ class TFCSHitCellMapping + ;
#pragma link C++ class TFCSHitCellMappingFCal + ;
#pragma link C++ class TFCSHitCellMappingWiggle + ;
#pragma link C++ class TFCSHitCellMappingWiggleEMB + ;
#pragma link C++ class TFCSEnergyRenormalization + ;
#pragma link C++ class TFCSVoxelHistoLateralCovarianceFluctuations + ;

#pragma link C++ class TFCSTruthState + ;
#pragma link C++ class TFCSExtrapolationState + ;
#pragma link C++ struct TFCSExtrapolationState::LayerPosition + ;
#pragma link C++ class std::vector<TFCSExtrapolationState::LayerPosition> + ;
#pragma link C++ class TFCSSimulationState + ;

#pragma link C++ class std::vector<TFCSTruthState> + ;
#pragma link C++ class std::vector<TFCSExtrapolationState> + ;
#pragma link C++ class std::vector<TFCSSimulationState> + ;

#endif
// clang-format on
