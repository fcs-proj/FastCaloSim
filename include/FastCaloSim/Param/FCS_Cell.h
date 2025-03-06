/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef FCS_Cell
#define FCS_Cell

#include <vector>

#include <Rtypes.h>
#include <TLorentzVector.h>

struct FCS_cell
{
  Long64_t cell_identifier;
  int sampling;
  float energy;
  float center_x;
  float center_y;
  float center_z;  // to be updated later
  auto operator<(const FCS_cell& rhs) const -> bool
  {
    return energy > rhs.energy;
  };
};

struct FCS_hit  // this is the FCS detailed hit
{
  Long64_t identifier;  // hit in the same tile cell can have two identifiers
                        // (for two PMTs)
  Long64_t cell_identifier;
  int sampling;  // calorimeter layer
  float hit_energy;  // energy is already scaled for the sampling fraction
  float hit_time;
  float hit_x;
  float hit_y;
  float hit_z;
  auto operator<(const FCS_hit& rhs) const -> bool
  {
    return hit_energy > rhs.hit_energy;
  };
  // float  hit_sampfrac;
};

struct FCS_g4hit  // this is the standard G4Hit
{
  Long64_t identifier;
  Long64_t cell_identifier;
  int sampling;
  float hit_energy;
  float hit_time;
  // float  hit_sampfrac;
  auto operator<(const FCS_g4hit& rhs) const -> bool
  {
    return hit_energy > rhs.hit_energy;
  };
};

struct FCS_matchedcell  // this is the matched structure for a single cell
{
  FCS_cell cell;
  std::vector<FCS_g4hit> g4hit;
  std::vector<FCS_hit> hit;
  inline void clear()
  {
    g4hit.clear();
    hit.clear();
  };
  inline auto scalingfactor() -> float
  {
    float hitsum = 0.;
    for (auto& i : hit) {
      hitsum += i.hit_energy;
    };
    return cell.energy / hitsum;
  };  // doesn't check for 0!
  auto operator<(const FCS_matchedcell& rhs) const -> bool
  {
    return cell.energy > rhs.cell.energy;
  };
  inline void sorthit() { std::sort(hit.begin(), hit.end()); };
  inline void sortg4hit() { std::sort(g4hit.begin(), g4hit.end()); };
  inline void sort()
  {
    sorthit();
    sortg4hit();
  };
  inline void time_trim(float timing_cut)
  {
    hit.erase(std::remove_if(hit.begin(),
                             hit.end(),
                             [&timing_cut](const FCS_hit& rhs)
                             { return rhs.hit_time > timing_cut; }),
              hit.end());
    g4hit.erase(std::remove_if(g4hit.begin(),
                               g4hit.end(),
                               [&timing_cut](const FCS_g4hit& rhs)
                               { return rhs.hit_time > timing_cut; }),
                g4hit.end());
  };
};

struct FCS_matchedcellvector  // this is the matched structure for the whole
                              // event (or single layer) - vector of
                              // FCS_matchedcell
{
  // Note that struct can have methods
  // Note the overloaded operator(s) to access the underlying vector
  std::vector<FCS_matchedcell> m_vector;
  inline auto GetLayer(int layer) -> std::vector<FCS_matchedcell>
  {
    std::vector<FCS_matchedcell> ret;
    for (auto& i : m_vector) {
      if (i.cell.sampling == layer) {
        ret.push_back(i);
      }
    };
    return ret;
  };
  inline auto operator[](unsigned int place) -> FCS_matchedcell
  {
    return m_vector[place];
  };
  inline auto size() const -> unsigned int { return m_vector.size(); };
  inline void push_back(FCS_matchedcell cell) { m_vector.push_back(cell); };
  inline void sort_cells() { std::sort(m_vector.begin(), m_vector.end()); };
  inline void sort()
  {
    std::sort(m_vector.begin(), m_vector.end());
    for (auto& i : m_vector) {
      i.sort();
    };
  };
  inline void time_trim(float timing_cut)
  {
    for (auto& i : m_vector) {
      i.time_trim(timing_cut);
    };
    m_vector.erase(std::remove_if(m_vector.begin(),
                                  m_vector.end(),
                                  [](const FCS_matchedcell& rhs)
                                  {
                                    return (rhs.hit.empty() && rhs.g4hit.empty()
                                            && fabs(rhs.cell.energy) < 1e-3);
                                  }),
                   m_vector.end());
  };
  inline auto scalingfactor() -> float
  {
    float cellsum = 0.;
    float hitsum = 0.;
    for (auto& i : m_vector) {
      cellsum += i.cell.energy;
      for (auto& j : i.hit) {
        hitsum += j.hit_energy;
      };
    };
    return cellsum / hitsum;
  };  // doesn't check for 0!
};

#endif
