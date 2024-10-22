// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

// ***************************************************************************
// Liquid Argon FCAL detector description package
// -----------------------------------------
//
//
// 10-Sep-2000 S.Simion   Handling of the FCAL read-out identifiers
//    Jan-2001 R.Sobie    Modify for persistency
//    Feb-2002 R.Sobie    Use same FCAL geometry files as simulation
// ***************************************************************************

#ifndef FCAL_CHANNELMAP_H
#define FCAL_CHANNELMAP_H

#include <map>
#include <string>
#include <vector>

/** This class contains the tube and tile maps for the FCAL <br>
 * A tile is of a set of FCAL tubes
 */
class FCAL_ChannelMap
{
public:
  //
  //	 Typedefs:
  //
  typedef unsigned int tileName_t;
  typedef unsigned int tubeID_t;

  //------------ TubePosition:  X,Y Position of an FCAL Tube---- //
  //                                                             //
  class TubePosition
  {  //
  public:  //
    TubePosition();  //
    TubePosition(tileName_t name,
                 float x,
                 float y,
                 const std::string& hvFT);  //
    tileName_t get_tileName() const;  //
    float x() const;  //
    float y() const;  //
    const std::string& getHVft() const;  // Gabe
  private:  //
    tileName_t m_tileName;  //
    float m_x;  //
    float m_y;  //
    std::string m_hvFT;  // Gabe
  };  //
  //-------------------------------------------------------------//

  /**  Constructors: */
  FCAL_ChannelMap(int itemp);

  typedef std::map<tubeID_t, TubePosition> tubeMap_t;
  typedef tubeMap_t::size_type tubemap_sizetype;
  typedef tubeMap_t::value_type tubemap_valuetype;
  typedef tubeMap_t::const_iterator tubemap_const_iterator;

  /** tubeMap access functions */
  tubemap_const_iterator tubemap_begin(int isam) const;
  tubemap_const_iterator tubemap_end(int isam) const;
  tubemap_sizetype tubemap_size(int isam) const;

  // Get the tube by copy number:
  tubemap_const_iterator getTubeByCopyNumber(int isam, int copyNo) const;

  /**  ---- For the new LArFCAL_ID Identifier
   */

  bool getTileID(int isam, float x, float y, int& eta, int& phi) const;

  /** For reconstruction, decoding of tile identifiers */
  float x(int isam, int eta, int phi) const;
  float y(int isam, int eta, int phi) const;

  void tileSize(int sam, int eta, int phi, float& dx, float& dy) const;

  void tileSize(int isam, int ntubes, float& dx, float& dy) const;

  /** print functions */
  void print_tubemap(int isam) const;

  /** set and get for the inversion flags
   */

  bool invert_x() const;
  bool invert_xy() const;
  void set_invert_x(bool);
  void set_invert_xy(bool);

  // Fill this.  The information comes from the database.
  void add_tube(const std::string& tileName,
                int mod,
                int id,
                int i,
                int j,
                double xCm,
                double yCm);  // original
  void add_tube(const std::string& tileName,
                int mod,
                int id,
                int i,
                int j,
                double xCm,
                double yCm,
                const std::string& hvFT);  // 29-03-07 include HV

  // Finish the job. Create the tile map.
  void finish();

  class TilePosition
  {
  public:
    TilePosition();
    TilePosition(float x, float y, int ntubes);
    float x() const;
    float y() const;
    unsigned int ntubes() const;

  private:
    float m_x;
    float m_y;
    unsigned int m_ntubes;
  };

  /** TileMap */
  typedef std::map<tileName_t, TilePosition> tileMap_t;
  typedef tileMap_t::size_type tileMap_sizetype;
  typedef tileMap_t::value_type tileMap_valuetype;
  typedef tileMap_t::const_iterator tileMap_const_iterator;

  // Methods to iterate over the tile map:
  tileMap_const_iterator begin(int isam) const
  {
    return m_tileMap[isam - 1].begin();
  }
  tileMap_const_iterator end(int isam) const
  {
    return m_tileMap[isam - 1].end();
  };

private:
  /** Geometrical parameters here, in CLHEP::cm please to be compatible with G3
   */
  static const double m_tubeSpacing[];
  double m_tubeDx[3];
  double m_tubeDy[3];
  double m_tileDx[3];
  double m_tileDy[3];
  mutable bool m_invert_x;  // Some geometry need x inverted
  mutable bool m_invert_xy;  // Some geometry need xy crossed

  tileMap_t m_tileMap[3];
  void create_tileMap(int isam);

  /** TubeMap */
  tubeMap_t m_tubeMap[3];
  std::vector<tubemap_const_iterator> m_tubeIndex[3];
};

// Tube Position functions

inline FCAL_ChannelMap::TubePosition::TubePosition()
    : m_tileName(0)
    , m_x(0)
    , m_y(0)
    , m_hvFT("")
{
}

inline FCAL_ChannelMap::TubePosition::TubePosition(tileName_t name,
                                                   float x,
                                                   float y,
                                                   const std::string& hvFT)
    : m_tileName(name)
    , m_x(x)
    , m_y(y)
    , m_hvFT(hvFT)
{
}

inline FCAL_ChannelMap::tileName_t FCAL_ChannelMap::TubePosition::get_tileName()
    const
{
  return m_tileName;
}

inline float FCAL_ChannelMap::TubePosition::x() const
{
  return m_x;
}

inline float FCAL_ChannelMap::TubePosition::y() const
{
  return m_y;
}

inline const std::string& FCAL_ChannelMap::TubePosition::getHVft() const
{
  return m_hvFT;
}

// Tile Position functions

inline FCAL_ChannelMap::TilePosition::TilePosition()
    : m_x(0)
    , m_y(0)
    , m_ntubes(0)
{
}

inline FCAL_ChannelMap::TilePosition::TilePosition(float x, float y, int ntubes)
    : m_x(x)
    , m_y(y)
    , m_ntubes(ntubes)
{
}

inline float FCAL_ChannelMap::TilePosition::x() const
{
  return m_x;
}

inline float FCAL_ChannelMap::TilePosition::y() const
{
  return m_y;
}

inline unsigned int FCAL_ChannelMap::TilePosition::ntubes() const
{
  return m_ntubes;
}

inline bool FCAL_ChannelMap::invert_x() const
{
  return m_invert_x;
}

inline bool FCAL_ChannelMap::invert_xy() const
{
  return m_invert_xy;
}

inline void FCAL_ChannelMap::set_invert_x(bool flag)
{
  m_invert_x = flag;
  return;
}

inline void FCAL_ChannelMap::set_invert_xy(bool flag)
{
  m_invert_xy = flag;
  return;
}

// #include "CLIDSvc/CLASS_DEF.h"
// CLASS_DEF(FCAL_ChannelMap, 74242524, 1)

#endif  // LARDETDESCR_FCAL_CHANNELMAP_H
