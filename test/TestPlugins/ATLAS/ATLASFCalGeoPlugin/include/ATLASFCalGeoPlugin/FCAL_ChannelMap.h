// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

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
  using tileName_t = unsigned int;
  using tubeID_t = unsigned int;

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
    auto get_tileName() const -> tileName_t;  //
    auto x() const -> float;  //
    auto y() const -> float;  //
    auto getHVft() const -> const std::string&;  // Gabe
  private:  //
    tileName_t m_tileName;  //
    float m_x;  //
    float m_y;  //
    std::string m_hvFT;  // Gabe
  };  //
  //-------------------------------------------------------------//

  /**  Constructors: */
  FCAL_ChannelMap(int itemp);

  using tubeMap_t = std::map<tubeID_t, TubePosition>;
  using tubemap_sizetype = tubeMap_t::size_type;
  using tubemap_valuetype = tubeMap_t::value_type;
  using tubemap_const_iterator = tubeMap_t::const_iterator;

  /** tubeMap access functions */
  auto tubemap_begin(int isam) const -> tubemap_const_iterator;
  auto tubemap_end(int isam) const -> tubemap_const_iterator;
  auto tubemap_size(int isam) const -> tubemap_sizetype;

  // Get the tube by copy number:
  auto getTubeByCopyNumber(int isam, int copyNo) const
      -> tubemap_const_iterator;

  /**  ---- For the new LArFCAL_ID Identifier
   */

  auto getTileID(int isam, float x, float y, int& eta, int& phi) const -> bool;

  /** For reconstruction, decoding of tile identifiers */
  auto x(int isam, int eta, int phi) const -> float;
  auto y(int isam, int eta, int phi) const -> float;

  void tileSize(int sam, int eta, int phi, float& dx, float& dy) const;

  void tileSize(int isam, int ntubes, float& dx, float& dy) const;

  /** print functions */
  void print_tubemap(int isam) const;

  /** set and get for the inversion flags
   */

  auto invert_x() const -> bool;
  auto invert_xy() const -> bool;
  void set_invert_x(bool /*flag*/);
  void set_invert_xy(bool /*flag*/);

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
    TilePosition(float x, float y, unsigned int ntubes);
    auto x() const -> float;
    auto y() const -> float;
    auto ntubes() const -> unsigned int;

  private:
    float m_x;
    float m_y;
    unsigned int m_ntubes;
  };

  /** TileMap */
  using tileMap_t = std::map<tileName_t, TilePosition>;
  using tileMap_sizetype = tileMap_t::size_type;
  using tileMap_valuetype = tileMap_t::value_type;
  using tileMap_const_iterator = tileMap_t::const_iterator;

  // Methods to iterate over the tile map:
  auto begin(int isam) const -> tileMap_const_iterator
  {
    return m_tileMap[isam - 1].begin();
  }
  auto end(int isam) const -> tileMap_const_iterator
  {
    return m_tileMap[isam - 1].end();
  }

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

inline auto FCAL_ChannelMap::TubePosition::get_tileName() const
    -> FCAL_ChannelMap::tileName_t
{
  return m_tileName;
}

inline auto FCAL_ChannelMap::TubePosition::x() const -> float
{
  return m_x;
}

inline auto FCAL_ChannelMap::TubePosition::y() const -> float
{
  return m_y;
}

inline auto FCAL_ChannelMap::TubePosition::getHVft() const -> const std::string&
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

inline FCAL_ChannelMap::TilePosition::TilePosition(float x,
                                                   float y,
                                                   unsigned int ntubes)
    : m_x(x)
    , m_y(y)
    , m_ntubes(ntubes)
{
}

inline auto FCAL_ChannelMap::TilePosition::x() const -> float
{
  return m_x;
}

inline auto FCAL_ChannelMap::TilePosition::y() const -> float
{
  return m_y;
}

inline auto FCAL_ChannelMap::TilePosition::ntubes() const -> unsigned int
{
  return m_ntubes;
}

inline auto FCAL_ChannelMap::invert_x() const -> bool
{
  return m_invert_x;
}

inline auto FCAL_ChannelMap::invert_xy() const -> bool
{
  return m_invert_xy;
}

inline void FCAL_ChannelMap::set_invert_x(bool flag)
{
  m_invert_x = flag;
}

inline void FCAL_ChannelMap::set_invert_xy(bool flag)
{
  m_invert_xy = flag;
}

// #include "CLIDSvc/CLASS_DEF.h"
// CLASS_DEF(FCAL_ChannelMap, 74242524, 1)

#endif  // LARDETDESCR_FCAL_CHANNELMAP_H
