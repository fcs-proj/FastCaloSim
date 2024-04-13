/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#include "FastCaloSim/DoubleArray.h"
#include "FastCaloSim/IntArray.h"
#include "FastCaloSim/TFCSFunction.h"
#include "FastCaloSim/TFCS1DFunction.h"
#include "FastCaloSim/TFCS1DFunctionHistogram.h"
#include "FastCaloSim/TFCS1DFunctionInt16Histogram.h"
#include "FastCaloSim/TFCS1DFunctionInt32Histogram.h"
#include "FastCaloSim/TFCS1DFunctionRegression.h"
#include "FastCaloSim/TFCS1DFunctionRegressionTF.h"
#include "FastCaloSim/TFCS1DFunctionSpline.h"
#include "FastCaloSim/TFCS1DFunctionTemplateHelpers.h"
#include "FastCaloSim/TFCS1DFunctionTemplateHistogram.h"
#include "FastCaloSim/TFCS1DFunctionTemplateInterpolationHistogram.h"
#include "FastCaloSim/TFCS2DFunction.h"
#include "FastCaloSim/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/TFCS2DFunctionTemplateHistogram.h"
#include "FastCaloSim/TFCS2DFunctionTemplateInterpolationHistogram.h"
#include "FastCaloSim/TFCS2DFunctionTemplateInterpolationExpHistogram.h"

#include "FastCaloSim/MLogging.h"
#include "FastCaloSim/TFCSParametrizationBase.h"
#include "FastCaloSim/TFCSParametrizationPlaceholder.h"
#include "FastCaloSim/TFCSParametrization.h"
#include "FastCaloSim/TFCSInvisibleParametrization.h"
#include "FastCaloSim/TFCSInitWithEkin.h"
#include "FastCaloSim/TFCSEnergyInterpolationHistogram.h"
#include "FastCaloSim/TFCSEnergyInterpolationLinear.h"
#include "FastCaloSim/TFCSEnergyInterpolationPiecewiseLinear.h"
#include "FastCaloSim/TFCSEnergyInterpolationSpline.h"
#include "FastCaloSim/TFCSParametrizationChain.h"
#include "FastCaloSim/TFCSParametrizationBinnedChain.h"
#include "FastCaloSim/TFCSParametrizationFloatSelectChain.h"
#include "FastCaloSim/TFCSParametrizationPDGIDSelectChain.h"
#include "FastCaloSim/TFCSParametrizationEbinChain.h"
#include "FastCaloSim/TFCSParametrizationEkinSelectChain.h"
#include "FastCaloSim/TFCSParametrizationEtaSelectChain.h"
#include "FastCaloSim/TFCSParametrizationAbsEtaSelectChain.h"

#include "FastCaloSim/TFCSEnergyParametrization.h"
#include "FastCaloSim/TFCSPCAEnergyParametrization.h"
#include "FastCaloSim/TFCSEnergyBinParametrization.h"
#include "FastCaloSim/TFCSGANXMLParameters.h"
#include "FastCaloSim/TFCSGANLWTNNHandler.h"
#include "FastCaloSim/TFCSGANEtaSlice.h"
#include "FastCaloSim/TFCSEnergyAndHitGANV2.h"
#include "FastCaloSim/TFCSPredictExtrapWeights.h"
#include "FastCaloSim/VNetworkBase.h"
#include "FastCaloSim/VNetworkLWTNN.h"
#include "FastCaloSim/TFCSSimpleLWTNNHandler.h"
#include "FastCaloSim/TFCSONNXHandler.h"
#include "FastCaloSim/TFCSNetworkFactory.h"

#include "FastCaloSim/TFCSLateralShapeParametrization.h"
#include "FastCaloSim/TFCSHistoLateralShapeGausLogWeight.h"
#include "FastCaloSim/TFCSHistoLateralShapeGausLogWeightHitAndMiss.h"

#include "FastCaloSim/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/TFCSLateralShapeParametrizationHitChain.h"
#include "FastCaloSim/TFCSLateralShapeParametrizationFluctChain.h"
#include "FastCaloSim/TFCSCenterPositionCalculation.h"
#include "FastCaloSim/TFCSFlatLateralShapeParametrization.h"
#include "FastCaloSim/TFCSHistoLateralShapeParametrization.h"
#include "FastCaloSim/TFCSHistoLateralShapeParametrizationFCal.h"
#include "FastCaloSim/TFCS2DFunctionLateralShapeParametrization.h"
#include "FastCaloSim/TFCSHistoLateralShapeWeight.h"
#include "FastCaloSim/TFCSHistoLateralShapeWeightHitAndMiss.h"
#include "FastCaloSim/TFCSLateralShapeParametrizationHitNumberFromE.h"
#include "FastCaloSim/TFCSLateralShapeTuning.h"
#include "FastCaloSim/TFCSHitCellMapping.h"
#include "FastCaloSim/TFCSHitCellMappingFCal.h"
#include "FastCaloSim/TFCSHitCellMappingWiggle.h"
#include "FastCaloSim/TFCSHitCellMappingWiggleEMB.h"
#include "FastCaloSim/TFCSEnergyRenormalization.h"
#include "FastCaloSim/TFCSVoxelHistoLateralCovarianceFluctuations.h"

#include "FastCaloSim/TFCSTruthState.h"
#include "FastCaloSim/TFCSExtrapolationState.h"
#include "FastCaloSim/TFCSSimulationState.h"


#ifdef __CLING__
// Standard preamble: turn off creation of dictionaries for "everything":
// we then turn it on only for the types we are interested in.
#  pragma link off all globals;
#  pragma link off all classes;
#  pragma link off all functions;
// Turn on creation of dictionaries for nested classes
#  pragma link C++ nestedclasses;

// Custom classes
#pragma link C++ class DoubleArray + ;
#pragma link C++ class IntArray + ;
#pragma link C++ class TFCSFunction + ;
#pragma link C++ class TFCS1DFunction + ;
#pragma link C++ class TFCS1DFunctionHistogram + ;
#pragma link C++ class TFCS1DFunctionInt16Histogram + ;
#pragma link C++ class TFCS1DFunctionInt32Histogram + ;
#pragma link C++ class TFCS1DFunctionRegression + ;
#pragma link C++ class TFCS1DFunctionRegressionTF + ;
#pragma link C++ class TFCS1DFunctionSpline + ;

/// Linkdefs needed for common classes for template based histograms
#pragma link C++ class TFCS1DFunction_Numeric < uint8_t, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < uint16_t, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < uint32_t, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < float, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < double, float> + ;
#pragma link C++ class TFCS1DFunction_Numeric < double, double> + ;

#pragma link C++ class TFCS1DFunction_Array < float> - ;
#pragma link C++ class TFCS1DFunction_Array < double> - ;
#pragma link C++ class TFCS1DFunction_Array < uint8_t> - ;
#pragma link C++ class TFCS1DFunction_Array < uint16_t> - ;
#pragma link C++ class TFCS1DFunction_Array < uint32_t> - ;

#pragma link C++ class TFCS1DFunction_HistogramContent < float, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < double, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < double, double> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < uint8_t, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < uint16_t, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramContent < uint32_t, float> + ;

#pragma link C++ class TFCS1DFunction_HistogramBinEdges < float, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramBinEdges < double, float> + ;
#pragma link C++ class TFCS1DFunction_HistogramBinEdges < double, double> + ;

#pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float,        \
    uint8_t, float> +                                                          \
    ;
#pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float,        \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunction_HistogramCompactBinEdges < float,        \
    uint32_t, float> +                                                         \
    ;

#pragma link C++ class TFCS1DFunction_HistogramInt8BinEdges + ;
#pragma link C++ class TFCS1DFunction_HistogramInt16BinEdges + ;
#pragma link C++ class TFCS1DFunction_HistogramInt32BinEdges + ;
#pragma link C++ class TFCS1DFunction_HistogramFloatBinEdges + ;
#pragma link C++ class TFCS1DFunction_HistogramDoubleBinEdges + ;

/// Linkdefs needed for 1D template based histograms
#pragma link C++ class TFCS1DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint8_t, float> +                                                          \
    ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint32_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint32_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    uint32_t, float> +                                                         \
    ;

#pragma link C++ class TFCS1DFunctionInt8Int8Histogram + ;
#pragma link C++ class TFCS1DFunctionInt8Int16Histogram + ;
#pragma link C++ class TFCS1DFunctionInt8Int32Histogram + ;
#pragma link C++ class TFCS1DFunctionInt16Int16Histogram + ;
#pragma link C++ class TFCS1DFunctionInt16Int32Histogram + ;
#pragma link C++ class TFCS1DFunctionInt32Int32Histogram + ;

#pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint8_t, float> +                                                          \
    ;
#pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class TFCS1DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint32_t, float> +                                                         \
    ;

#pragma link C++ class TFCS1DFunctionInt8Int8InterpolationHistogram + ;
#pragma link C++ class TFCS1DFunctionInt8Int16InterpolationHistogram + ;
#pragma link C++ class TFCS1DFunctionInt16Int16InterpolationHistogram + ;
#pragma link C++ class TFCS1DFunctionInt16Int32InterpolationHistogram + ;

/// Linkdefs needed for 2D template based histograms
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateHistogram <                       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;

#pragma link C++ class TFCS2DFunctionInt8Int8Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int8Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int8Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int16Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int16Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int16Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int32Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int32Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int32Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int8Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int8Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int8Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int16Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int16Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int16Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int32Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int32Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int32Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int8Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int8Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int8Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int16Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int16Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int16Int32Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int32Int8Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int32Int16Histogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int32Int32Histogram + ;

#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationHistogram <          \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;

#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt8BinEdges,                                      \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt16BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    TFCS1DFunction_HistogramInt32BinEdges,                                     \
    TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;

#pragma link C++ class TFCS2DFunctionInt8Int8Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int8Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int8Int32InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int16Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int16Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int16Int32InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int32Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int32Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt8Int32Int32InterpolationHistogram + ;

#pragma link C++ class TFCS2DFunctionInt16Int8Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int8Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int8Int32InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int16Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int16Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int16Int32InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int32Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int32Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt16Int32Int32InterpolationHistogram + ;

#pragma link C++ class TFCS2DFunctionInt32Int8Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int8Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int8Int32InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int16Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int16Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int16Int32InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int32Int8InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int32Int16InterpolationHistogram + ;
#pragma link C++ class TFCS2DFunctionInt32Int32Int32InterpolationHistogram + ;

/// End Linkdefs needed for template based histograms

#pragma link C++ class TFCS2DFunction + ;
#pragma link C++ class TFCS2DFunctionHistogram + ;
#pragma link C++ class ISF_FCS::MLogging + ;
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
#pragma link C++ class TFCSGANXMLParameters + ;
#pragma link C++ class TFCSGANLWTNNHandler - ;
#pragma link C++ class TFCSGANEtaSlice + ;
#pragma link C++ class TFCSEnergyAndHitGANV2 + ;
#pragma link C++ class TFCSPredictExtrapWeights - ;

#pragma link C++ class VNetworkBase + ;
#pragma link C++ class VNetworkLWTNN + ;
#pragma link C++ class TFCSSimpleLWTNNHandler - ;
#pragma link C++ class TFCSONNXHandler - ;
#pragma link C++ class TFCSNetworkFactory + ;

#pragma link C++ class TFCSLateralShapeParametrization + ;
#pragma link C++ class TFCSLateralShapeParametrizationHitBase + ;
#pragma link C++ class TFCSLateralShapeParametrizationHitChain + ;
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
#pragma link C++ class TFCSSimulationState + ;


#endif  // __CLING__