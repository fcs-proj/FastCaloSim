// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <filesystem>
#include <iostream>
#include <string>

#include <dlfcn.h>

template<class Plugin>
class GeoPluginLoader
{
public:
  // Constructor:
  GeoPluginLoader();

  // Destructor:
  ~GeoPluginLoader();

  // load G4Solid plugin
  Plugin* load(const std::string& path) const;

private:
  GeoPluginLoader(const GeoPluginLoader&) = delete;
  GeoPluginLoader& operator=(const GeoPluginLoader&) = delete;
};

template<class Plugin>
GeoPluginLoader<Plugin>::GeoPluginLoader() = default;

template<class Plugin>
GeoPluginLoader<Plugin>::~GeoPluginLoader() = default;

template<class Plugin>
Plugin* GeoPluginLoader<Plugin>::load(const std::string& pString) const
{
  std::string bNameString = std::filesystem::path(pString)
                                .filename()
                                .string();  // Strip off the directory
  bNameString = bNameString.substr(3);  // Strip off leading "lib"
  bNameString =
      bNameString.substr(0, bNameString.find("."));  // Strip off extensions

  std::string createFunctionName = std::string("create") + bNameString;

  //
  // Loads the library:
  //
  void* handle = dlopen(pString.c_str(), RTLD_NOW);
  if (!handle) {
    std::cerr << dlerror() << '\n';
    return nullptr;
  }

  //
  // Gets the function
  //
  void* f = dlsym(handle, createFunctionName.c_str());
  if (!f) {
    std::cerr << dlerror() << '\n';
    return nullptr;
  }

  using CreationMethod = void* (*)();
  auto F = reinterpret_cast<CreationMethod>(f);

  //
  //
  //
  Plugin* factory = static_cast<Plugin*>(F());
  return factory;
}
