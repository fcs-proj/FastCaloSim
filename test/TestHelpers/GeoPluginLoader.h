// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>

#include <dlfcn.h>

template<class Plugin>
class GeoPluginLoader
{
public:
  // Constructor:
  GeoPluginLoader() = default;

  // Destructor:
  ~GeoPluginLoader() = default;

  GeoPluginLoader(const GeoPluginLoader&) = delete;
  auto operator=(const GeoPluginLoader&) -> GeoPluginLoader& = delete;
  GeoPluginLoader(GeoPluginLoader&&) = delete;
  auto operator=(GeoPluginLoader&&) -> GeoPluginLoader& = delete;

  // Load the plugin library pString and create its plugin; throws
  // std::runtime_error if the library, its create function or the plugin is
  // missing
  auto load(const std::string& pString) const -> Plugin*;
};

template<class Plugin>
auto GeoPluginLoader<Plugin>::load(const std::string& pString) const -> Plugin*
{
  std::string bNameString = std::filesystem::path(pString)
                                .filename()
                                .string();  // Strip off the directory
  bNameString = bNameString.substr(3);  // Strip off leading "lib"
  bNameString =
      bNameString.substr(0, bNameString.find('.'));  // Strip off extensions

  std::string const createFunctionName = std::string("create") + bNameString;

  //
  // Loads the library:
  //
  void* handle = dlopen(pString.c_str(), RTLD_NOW);
  if (handle == nullptr) {
    throw std::runtime_error("Cannot load plugin library: "
                             + std::string(dlerror()));
  }

  //
  // Gets the function
  //
  void* f = dlsym(handle, createFunctionName.c_str());
  if (f == nullptr) {
    throw std::runtime_error("Cannot find " + createFunctionName + " in "
                             + pString + ": " + std::string(dlerror()));
  }

  using CreationMethod = void* (*)();
  auto F = reinterpret_cast<CreationMethod>(f);

  //
  //
  //
  auto* factory = static_cast<Plugin*>(F());
  if (factory == nullptr) {
    throw std::runtime_error(createFunctionName + " in " + pString
                             + " returned no plugin");
  }
  return factory;
}
