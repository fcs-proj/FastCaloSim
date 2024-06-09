/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#pragma once

#include <iostream>
#include <string>

#include <dlfcn.h>
#include <libgen.h>

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
GeoPluginLoader<Plugin>::GeoPluginLoader()
{
}

template<class Plugin>
GeoPluginLoader<Plugin>::~GeoPluginLoader()
{
}

template<class Plugin>
Plugin* GeoPluginLoader<Plugin>::load(const std::string& pString) const
{
  std::string bNameString =
      basename((char*)pString.c_str());  // Strip off the directory
  bNameString = bNameString.substr(3);  // Strip off leading "lib"
  bNameString =
      bNameString.substr(0, bNameString.find("."));  // Strip off extensions

  std::string createFunctionName = std::string("create") + bNameString;

  //
  // Loads the library:
  //
  void* handle = dlopen(pString.c_str(), RTLD_NOW);
  if (!handle) {
    std::cerr << dlerror() << std::endl;
    return nullptr;
  }

  //
  // Gets the function
  //
  void* f = dlsym(handle, createFunctionName.c_str());
  if (!f) {
    std::cerr << dlerror() << std::endl;
    return nullptr;
  }

  typedef void* (*CreationMethod)();
  CreationMethod F = (CreationMethod)f;

  //
  //
  //
  Plugin* factory = (Plugin*)F();
  return factory;
}