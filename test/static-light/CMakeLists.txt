# Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

# Verifies the opposite corner from test/static-embed: a consumer that uses
# only the lightweight logging utility (FastCaloSim/Core/MLogging.h) must be
# able to build the relevant part of the FastCaloSim archive into its own
# shared library WITHOUT FastCaloSim's heavy dependency closure (ROOT
# RDataFrame, HDF5, lwtnn, onnxruntime, libspatialindex, libxml2, ...).
#
# The archive is linked as a raw file (not the CMake target) so that none of
# the target's interface link dependencies reach the link line, and
# --no-undefined turns any reference that escapes the light path -- for
# example the logging class dragging in the monolithic ROOT dictionary
# translation unit -- into a hard link error.
add_library(fcs_light_consumer SHARED light_consumer.cxx)
target_include_directories(
    fcs_light_consumer PRIVATE
    "$<TARGET_PROPERTY:FastCaloSim::FastCaloSim,INTERFACE_INCLUDE_DIRECTORIES>"
)
target_link_libraries(
    fcs_light_consumer PRIVATE
    "$<TARGET_FILE:FastCaloSim::FastCaloSim>"
)
target_link_options(fcs_light_consumer PRIVATE "LINKER:--no-undefined")
add_dependencies(fcs_light_consumer FastCaloSim_FastCaloSim)
deactivate_checks(fcs_light_consumer)

# A thin driver that emits a log line through the absorbed logging class.
add_executable(static_light_check main.cxx)
target_link_libraries(static_light_check PRIVATE fcs_light_consumer)
deactivate_checks(static_light_check)

add_test(NAME StaticLightLogging COMMAND static_light_check)
