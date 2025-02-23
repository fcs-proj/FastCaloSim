// clang-format off
// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#include "FastCaloSim/Core/DoubleArray.h"
#include "FastCaloSim/Core/IntArray.h"
#include "FastCaloSim/Core/TFCSFunction.h"
#include "FastCaloSim/Core/TFCS1DFunction.h"
#include "FastCaloSim/Core/TFCS1DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionInt16Histogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionInt32Histogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionRegression.h"
#include "FastCaloSim/Core/TFCS1DFunctionRegressionTF.h"
#include "FastCaloSim/Core/TFCS1DFunctionSpline.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateHelpers.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateHistogram.h"
#include "FastCaloSim/Core/TFCS1DFunctionTemplateInterpolationHistogram.h"
#include "FastCaloSim/Core/TFCS2DFunction.h"
#include "FastCaloSim/Core/TFCS2DFunctionHistogram.h"
#include "FastCaloSim/Core/TFCS2DFunctionTemplateHistogram.h"
#include "FastCaloSim/Core/TFCS2DFunctionTemplateInterpolationHistogram.h"
#include "FastCaloSim/Core/TFCS2DFunctionTemplateInterpolationExpHistogram.h"

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
#ifndef __FastCaloSimNoLWTNN__
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
#endif

#include "FastCaloSim/Core/TFCSLateralShapeParametrization.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeGausLogWeight.h"
#include "FastCaloSim/Core/TFCSHistoLateralShapeGausLogWeightHitAndMiss.h"

#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitBase.h"
#include "FastCaloSim/Core/TFCSLateralShapeParametrizationHitChain.h"
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

#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSSimulationState.h"

#ifdef __CINT__
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
#ifndef __FastCaloSimNoLWTNN__
#pragma link C++ class TFCSGANXMLParameters + ;
#pragma link C++ class TFCSGANLWTNNHandler - ;
#pragma link C++ class TFCSGANEtaSlice + ;
#pragma link C++ class TFCSEnergyAndHitGANV2 + ;
#pragma link C++ class TFCSPredictExtrapWeights - ;

#pragma link C++ class VNetworkBase + ;
#pragma link C++ class VNetworkLWTNN + ;
#pragma link C++ class TFCSSimpleLWTNNHandler - ;
#pragma link C++ class TFCSONNXHandler - ;
#endif

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

#pragma link C++ class TFCSTruthState + ;
#pragma link C++ class TFCSExtrapolationState + ;
#pragma link C++ class TFCSSimulationState + ;

#pragma link C++ class std::vector<TFCSTruthState> + ;
#pragma link C++ class std::vector<TFCSExtrapolationState> + ;
#pragma link C++ class std::vector<TFCSSimulationState> + ;

#endif
// clang-format on
