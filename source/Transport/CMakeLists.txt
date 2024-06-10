file(GLOB FastCaloSimTransport_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimTransport_SOURCES}
)