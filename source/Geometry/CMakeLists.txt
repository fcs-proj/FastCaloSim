# Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

file(GLOB FastCaloSimGeometry_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimGeometry_SOURCES}
)
