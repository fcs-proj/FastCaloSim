# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

file(GLOB FastCaloSimExtrapolation_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimExtrapolation_SOURCES}
)
