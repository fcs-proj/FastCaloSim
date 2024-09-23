/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef FastCaloSimCaloExtrapolation_H
#define FastCaloSimCaloExtrapolation_H

#include <CLHEP/Vector/ThreeVector.h>
#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Core/TFCSExtrapolationState.h"

class TFCSTruthState;
class G4FieldTrack;
class CaloGeo;

// use CLHEP vector
using Vector3D = CLHEP::Hep3Vector;

struct CylinderIntersections
{
  Vector3D first;
  Vector3D second;
  unsigned int number;
};

class FASTCALOSIM_EXPORT FastCaloSimCaloExtrapolation : public ISF_FCS::MLogging
{
public:
  FastCaloSimCaloExtrapolation();
  ~FastCaloSimCaloExtrapolation() = default;

  void set_geometry(CaloGeo* geo) { m_geo = geo; };

  enum HITPOSITION
  {
    INSIDE,  // hit position is inside cylinder bounds
    OUTSIDE,  // hit position is outside cylinder bounds
    ON  // hit position is on cylinder bounds
  };

  virtual void extrapolate(TFCSExtrapolationState& result,
                           const TFCSTruthState* truth,
                           const std::vector<G4FieldTrack>& caloSteps) const;

private:
  /*Main extrapolation methods*/

  /// Finds best extrapolation extPos from the caloSteps for a cylinder defined
  /// by radius cylR and half-length cylZ as well as corresponding momentum
  /// direction
  bool extrapolateToCylinder(const std::vector<G4FieldTrack>& caloSteps,
                             float cylR,
                             float cylZ,
                             Vector3D& extPos,
                             Vector3D& momDir) const;
  /// Extrapolates to ID using three uniquely defined cylinder surfaces
  void extrapolateToID(TFCSExtrapolationState& result,
                       const std::vector<G4FieldTrack>& caloSteps,
                       const TFCSTruthState* truth) const;
  /// Extrapolates to all other layers of the calorimeter
  void extrapolateToLayers(TFCSExtrapolationState& result,
                           const std::vector<G4FieldTrack>& caloSteps,
                           const TFCSTruthState* truth) const;

  /*Extrapolator helper methods*/

  /// Finds Point of Closest Approach (PCA) on the cylinder defined by radius
  /// cylR and half-length cylZ of a line segment spanned by two hit positions
  /// to a cylinder
  void findPCA(float cylR,
               float cylZ,
               Vector3D& hitPos1,
               Vector3D& hitPos2,
               Vector3D& PCA) const;
  /// Computes the distance between a point and the line segment spanned by
  /// hitPos1 and hitPos2
  static double getPointLineSegmentDistance(Vector3D& point,
                                            Vector3D& hitPos1,
                                            Vector3D& hitPos2);
  /// Finds PCA iteratively given two bounds A and B on a line segment, used for
  /// (rare) cases with no easy analytical solutions
  void getIterativePCA(float cylR,
                       float cylZ,
                       Vector3D& BoundA,
                       Vector3D& BoundB,
                       Vector3D& PCA) const;
  /// Returns true if point lies on the line segment spanned by hitPos1 and
  /// hitPos2, otherwise returns false
  static bool isOnSegment(Vector3D& point,
                          Vector3D& hitPos1,
                          Vector3D& hitPos2);
  /// Computes intersection between the (infinite) line spanned by pointA and
  /// pointB with the positive (negative) endcap of a cylinder, returns true if
  /// intersection is found
  static bool cylinderEndcapIntersection(float cylR,
                                         float cylZ,
                                         bool positiveEndcap,
                                         Vector3D& pointA,
                                         Vector3D& pointB,
                                         Vector3D& intersection);
  /*!Extrapolates position on cylinder by finding intersections of subsequent
   hit positions, intersection is considered if we detect a travel through the
   surface with the line segment or we find a forward intersection (in the
   travel direction of the particle) which lies on the line segment, returns
   false if no such position is found*/
  bool extrapolateWithIntersection(const std::vector<G4FieldTrack>& caloSteps,
                                   float cylR,
                                   float cylZ,
                                   Vector3D& extPos,
                                   Vector3D& momDir) const;
  /// Extrapolates to the cylinder using the PCA to the polygon spanned by the
  /// individual line segments from the caloSteps
  bool extrapolateWithPCA(const std::vector<G4FieldTrack>& caloSteps,
                          float cylR,
                          float cylZ,
                          Vector3D& extPos,
                          Vector3D& momDir) const;
  /// Returns true if the line segment spanned by hitPos1 and hitPos2 crosses
  /// the cylinder surface, false otherwise
  static bool doesTravelThroughSurface(float cylR,
                                       float cylZ,
                                       Vector3D& hitPos1,
                                       Vector3D& hitPos2);
  /// Returns ID of more sensible intersection between line segment spanned by
  /// hitPos1 and hitPos2 and cylinder
  int whichIntersection(float cylR,
                        float cylZ,
                        Vector3D& hitPos1,
                        Vector3D& hitPos2,
                        Vector3D& intersectionA,
                        Vector3D intersectionB) const;
  /// Analytically computes 2D intersections between circle of radius circR and
  /// (infinite) line spanned by pointA nad pointB
  int circleLineIntersection2D(float circR,
                               Vector3D& pointA,
                               Vector3D& pointB,
                               Vector3D& intersectA,
                               Vector3D& intersectB) const;
  /// Analytically computes the intersection between the (infinite) line defined
  /// by pointA and pointB and the cylinder cover (without endcaps)
  int cylinderLineIntersection(float cylR,
                               float cylZ,
                               Vector3D& pointA,
                               Vector3D& pointB,
                               Vector3D& intersectA,
                               Vector3D& intersectB) const;
  /// Checks if position of hitPos is inside, outside or on the cylinder bounds
  static enum HITPOSITION whereOnCylinder(float cylR,
                                          float cylZ,
                                          Vector3D& hitPos);
  /// Projects position hitPos onto the cylinder surface and returns projected
  /// position
  static Vector3D projectOnCylinder(float cylR, float cylZ, Vector3D& hitPos);
  /// Analytically computes the intersection between the (infinite) line spanned
  /// by hitPos1 and hitPos2 with a cylinder
  CylinderIntersections getCylinderIntersections(float cylR,
                                                 float cylZ,
                                                 Vector3D& hitPos1,
                                                 Vector3D& hitPos2) const;

  auto get3DLinePCA(const Vector3D& posA,
                    const Vector3D& dirA,
                    const Vector3D& posB,
                    const Vector3D& dirB) const -> Vector3D;

  // uniquely defined ID-Calo surfaces
  std::vector<double> m_CaloBoundaryR {1148.0, 120.0, 41.0};
  std::vector<double> m_CaloBoundaryZ {3550.0, 4587.0, 4587.0};

  // FastCaloSim geometry helper
  CaloGeo* m_geo;

  /** calculates the distance between two point in 3D space */
  static inline auto distance(const Vector3D& p1, const Vector3D& p2) -> double
  {
    return (p2 - p1).mag();
  }

  /** checks if a vector is zero */
  static inline auto isZero(const Vector3D& v) -> bool
  {
    return v.x() == 0.0 && v.y() == 0.0 && v.z() == 0.0;
  }
};

#endif  // FastCaloSimCaloExtrapolation_H
