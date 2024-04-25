/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FSmap_H
#define ISF_FSmap_H

#include <map>

#include <math.h>

template<class _Key, class _Tp>
class FSmap : public std::map<_Key, _Tp>
{
public:
  typedef _Key key_type;
  typedef typename std::map<_Key, _Tp>::iterator iterator;
  typedef typename std::map<_Key, _Tp>::const_iterator const_iterator;

  const_iterator find_closest(const key_type& k) const
  {
    if (std::map<_Key, _Tp>::size() == 0)
      return std::map<_Key, _Tp>::end();
    const_iterator i = std::map<_Key, _Tp>::lower_bound(k);
    if (i == std::map<_Key, _Tp>::end()) {
      --i;
      return i;
    }
    if (i == std::map<_Key, _Tp>::begin())
      return i;

    const_iterator il = i;
    --il;

    if (fabs(il->first - k) > fabs(i->first - k))
      return i;
    else
      return il;
  };
  iterator find_closest(const key_type& k)
  {
    if (std::map<_Key, _Tp>::size() == 0)
      return std::map<_Key, _Tp>::end();
    iterator i = std::map<_Key, _Tp>::lower_bound(k);
    if (i == std::map<_Key, _Tp>::end()) {
      --i;
      return i;
    }
    if (i == std::map<_Key, _Tp>::begin())
      return i;

    iterator il = i;
    --il;

    if (fabs(il->first - k) > fabs(i->first - k))
      return i;
    else
      return il;
  };
};

#endif
