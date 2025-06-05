/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/
#ifndef G4ATLASTOOLS_THREADLOCALHOLDER_H
#define G4ATLASTOOLS_THREADLOCALHOLDER_H

/// @file ThreadLocalHolder.h
/// @brief Defines template helper classes for thread-local storage.
///
/// The helpers here provide an implementation of thread-local storage which
/// allows for non-static dynamic memory management as well as element
/// iteration. This extra functionality makes these containers much slower than
/// the builtin thread_local and should only be used when needed.
///
/// @todo Should we optimize these implementations for non-MT environments? One
/// could alternatively just optimize the clients where necessary. I could
/// replace these classes with a simple pointer wrapper if ATHENAHIVE and
/// G4MULTITHREADED are not defined, but the map iteration functionality would
/// need an awkward implementation.
///

// System includes
#include <thread>
#include <utility>

// Other includes
#include "tbb/concurrent_unordered_map.h"

namespace FastCaloSim::Transport::thread_utils
{

/// @class ThreadLocalHolder
/// @brief A thread-local storage wrapper.
///
/// This container is implemented as a wrapper for a concurrent map
/// keyed by std thread ID. It is thus fully thread-safe (in theory).
///
/// This container is non-owning. It doesn't clean up memory.
/// For a corresponding container which owns its objects, see
/// ThreadLocalOwner.
///
/// @author Steve Farrell <Steven.Farrell@cern.ch>
///
template<class T>
class ThreadLocalHolder
{
public:
  using MapKey_t = std::thread::id;
  using MapVal_t = T*;
  using MapHash_t = std::hash<MapKey_t>;
  using ThreadMap_t =
      tbb::concurrent_unordered_map<MapKey_t, MapVal_t, MapHash_t>;

  /// Get the object of the current thread.
  T* get()
  {
    auto itr = m_threadMap.find(std::this_thread::get_id());
    if (itr == m_threadMap.end())
      return nullptr;
    return itr->second;
  }

  /// Assign the object of the current thread.
  void set(T* obj)
  {
    const auto tid = std::this_thread::get_id();
    m_threadMap.insert(std::make_pair(tid, obj));
  }

  /// Constant access for iteration, etc.
  const ThreadMap_t& getMap() const { return m_threadMap; }

protected:
  /// The wrapped thread-local storage container
  ThreadMap_t m_threadMap;

};  // class ThreadLocalHolder

/// @class ThreadLocalOwner
/// @brief A ThreadLocalHolder which owns its objects.
///
/// This class merely adds deletion of the map elements in the destructor.
///
/// @author Steve Farrell <Steven.Farrell@cern.ch>
///
template<class T>
class ThreadLocalOwner : public ThreadLocalHolder<T>
{
public:
  /// Destructor will clean up the storage
  ~ThreadLocalOwner()
  {
    for (auto& mapPair : this->m_threadMap)
      delete mapPair.second;
  }
};  // class ThreadLocalOwner

}  // namespace FastCaloSim::Transport::thread_utils

#endif
