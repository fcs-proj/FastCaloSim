# Locate the pinned libxml-2.0 API via pkg-config and expose it as a target.
find_package(PkgConfig REQUIRED)
pkg_check_modules(LibXml2 REQUIRED IMPORTED_TARGET libxml-2.0)
