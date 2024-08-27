file(GLOB FastCaloSimCore_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimCore_SOURCES}
)