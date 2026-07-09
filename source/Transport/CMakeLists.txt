# Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

file(GLOB FastCaloSimTransport_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimTransport_SOURCES}
)
