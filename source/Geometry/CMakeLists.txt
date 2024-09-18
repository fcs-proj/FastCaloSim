file(GLOB FastCaloSimGeometry_SOURCES "*.cxx")
target_sources(
    FastCaloSim_FastCaloSim
    PRIVATE
    ${FastCaloSimGeometry_SOURCES}
)
