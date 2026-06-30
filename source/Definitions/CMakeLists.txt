# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

file(GLOB FastCaloSimDefinitions_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimDefinitions_SOURCES}
)
