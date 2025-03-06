# Copyright (c) 2025 CERN for the benefit of the FastCaloSim project

file(GLOB FastCaloSimParam_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimParam_SOURCES}
)
