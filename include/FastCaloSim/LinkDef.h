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
#pragma link C++ class FastCaloSim::Core::DoubleArray + ;
#pragma link C++ class FastCaloSim::Core::IntArray + ;
#pragma link C++ class FastCaloSim::Core::TFCSFunction + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionRegression + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionRegressionTF + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionSpline + ;

/// Linkdefs needed for common classes for template based histograms
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Numeric < uint8_t, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Numeric < uint16_t, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Numeric < uint32_t, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Numeric < float, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Numeric < double, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Numeric < double, double> + ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Array < float> - ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Array < double> - ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Array < uint8_t> - ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Array < uint16_t> - ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_Array < uint32_t> - ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramContent < float, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramContent < double, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramContent < double, double> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramContent < uint8_t, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramContent < uint16_t, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramContent < uint32_t, float> + ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramBinEdges < float, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramBinEdges < double, float> + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramBinEdges < double, double> + ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramCompactBinEdges < float,        \
    uint8_t, float> +                                                          \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramCompactBinEdges < float,        \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramCompactBinEdges < float,        \
    uint32_t, float> +                                                         \
    ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramFloatBinEdges + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunction_HistogramDoubleBinEdges + ;

/// Linkdefs needed for 1D template based histograms
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint8_t, float> +                                                          \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint32_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint32_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    uint32_t, float> +                                                         \
    ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt8Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt8Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt8Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt16Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt16Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt32Int32Histogram + ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint8_t, float> +                                                          \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint16_t, float> +                                                         \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    uint32_t, float> +                                                         \
    ;

#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt8Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt8Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt16Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS1DFunctionInt16Int32InterpolationHistogram + ;

/// Linkdefs needed for 2D template based histograms
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateHistogram <                       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;

#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int32Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int8Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int16Histogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int32Histogram + ;

#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationHistogram <          \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;

#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges,                                      \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint8_t, float> +                    \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint16_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt8BinEdges, uint32_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt16BinEdges, uint32_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint8_t, float> +                   \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint16_t, float> +                  \
    ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionTemplateInterpolationExpHistogram <       \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges,                                     \
    FastCaloSim::Core::TFCS1DFunction_HistogramInt32BinEdges, uint32_t, float> +                  \
    ;

#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int8Int32InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int16Int32InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt8Int32Int32InterpolationHistogram + ;

#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int8Int32InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int16Int32InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt16Int32Int32InterpolationHistogram + ;

#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int8Int32InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int16Int32InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int8InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int16InterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionInt32Int32Int32InterpolationHistogram + ;

/// End Linkdefs needed for template based histograms

#pragma link C++ class FastCaloSim::Core::TFCS2DFunction + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionHistogram + ;
#pragma link C++ class FastCaloSim::Core::MLogging + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationBase + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationPlaceholder + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSInvisibleParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSInitWithEkin + ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyInterpolationHistogram + ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyInterpolationLinear + ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyInterpolationPiecewiseLinear - ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyInterpolationSpline + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationChain - ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationBinnedChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationFloatSelectChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationPDGIDSelectChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationEbinChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationEkinSelectChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationEtaSelectChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSParametrizationAbsEtaSelectChain + ;

#pragma link C++ class FastCaloSim::Core::TFCSEnergyParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSPCAEnergyParametrization - ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyBinParametrization + ;
#ifndef __FastCaloSimNoLWTNN__
#pragma link C++ class FastCaloSim::Core::TFCSGANXMLParameters + ;
#pragma link C++ class FastCaloSim::Core::TFCSGANLWTNNHandler - ;
#pragma link C++ class FastCaloSim::Core::TFCSGANEtaSlice + ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyAndHitGANV2 + ;
#pragma link C++ class FastCaloSim::Core::TFCSPredictExtrapWeights - ;

#pragma link C++ class FastCaloSim::Core::VNetworkBase + ;
#pragma link C++ class FastCaloSim::Core::VNetworkLWTNN + ;
#pragma link C++ class FastCaloSim::Core::TFCSSimpleLWTNNHandler - ;
#pragma link C++ class FastCaloSim::Core::TFCSONNXHandler - ;
#endif

#pragma link C++ class FastCaloSim::Core::TFCSLateralShapeParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSLateralShapeParametrizationHitBase + ;
#pragma link C++ class FastCaloSim::Core::TFCSLateralShapeParametrizationHitChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSLateralShapeParametrizationFluctChain + ;
#pragma link C++ class FastCaloSim::Core::TFCSCenterPositionCalculation + ;
#pragma link C++ class FastCaloSim::Core::TFCSFlatLateralShapeParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSHistoLateralShapeParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSHistoLateralShapeParametrizationFCal + ;
#pragma link C++ class FastCaloSim::Core::TFCS2DFunctionLateralShapeParametrization + ;
#pragma link C++ class FastCaloSim::Core::TFCSHistoLateralShapeWeight + ;
#pragma link C++ class FastCaloSim::Core::TFCSHistoLateralShapeGausLogWeight + ;
#pragma link C++ class FastCaloSim::Core::TFCSHistoLateralShapeWeightHitAndMiss + ;
#pragma link C++ class FastCaloSim::Core::TFCSHistoLateralShapeGausLogWeightHitAndMiss + ;
#pragma link C++ class FastCaloSim::Core::TFCSLateralShapeParametrizationHitNumberFromE + ;
#pragma link C++ class FastCaloSim::Core::TFCSLateralShapeTuning + ;
#pragma link C++ class FastCaloSim::Core::TFCSHitCellMapping + ;
#pragma link C++ class FastCaloSim::Core::TFCSHitCellMappingFCal + ;
#pragma link C++ class FastCaloSim::Core::TFCSHitCellMappingWiggle + ;
#pragma link C++ class FastCaloSim::Core::TFCSHitCellMappingWiggleEMB + ;
#pragma link C++ class FastCaloSim::Core::TFCSEnergyRenormalization + ;

#pragma link C++ class FastCaloSim::Core::TFCSTruthState + ;
#pragma link C++ class FastCaloSim::Core::TFCSExtrapolationState + ;
#pragma link C++ class FastCaloSim::Core::TFCSSimulationState + ;

#pragma link C++ class std::vector<FastCaloSim::Core::TFCSTruthState> + ;
#pragma link C++ class std::vector<FastCaloSim::Core::TFCSExtrapolationState> + ;
#pragma link C++ class std::vector<FastCaloSim::Core::TFCSSimulationState> + ;

#endif
// clang-format on
