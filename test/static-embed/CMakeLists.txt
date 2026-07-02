# Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

# Verifies that FastCaloSim can be absorbed, whole-archive, into a single shared
# library and that its ROOT dictionary still self-registers at runtime. This is
# the configuration an embedding consumer uses when it links FastCaloSim as a
# static, position-independent archive into one of its own shared libraries.
#
# The whole-archive link feature ($<LINK_LIBRARY:WHOLE_ARCHIVE,...>) needs
# CMake >= 3.24.
if(CMAKE_VERSION VERSION_LESS 3.24)
  message(FATAL_ERROR
    "FastCaloSim_BUILD_EMBED_TEST requires CMake >= 3.24 for whole-archive linking")
endif()

# Absorb FastCaloSim into a single shared library. Whole-archive keeps the ROOT
# dictionary translation unit (which self-registers the TFCS* classes) from being
# dropped by the linker as unreferenced when FastCaloSim is a static archive.
# Building this shared library against a static archive also exercises that the
# archive is position-independent.
add_library(fcs_embed SHARED probe.cxx)
target_link_libraries(
    fcs_embed PRIVATE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,FastCaloSim::FastCaloSim>"
)
deactivate_checks(fcs_embed)

# A thin driver that asks whether the dictionaries are present in the absorbed
# library. Kept free of any FastCaloSim/ROOT headers so it only depends on the
# absorbing shared library.
add_executable(static_embed_check main.cxx)
target_link_libraries(static_embed_check PRIVATE fcs_embed)
deactivate_checks(static_embed_check)

add_test(NAME StaticEmbedDictionary COMMAND static_embed_check)
