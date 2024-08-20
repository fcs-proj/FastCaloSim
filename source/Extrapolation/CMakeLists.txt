file(GLOB FastCaloSimExtrapolation_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimExtrapolation_SOURCES}
)
