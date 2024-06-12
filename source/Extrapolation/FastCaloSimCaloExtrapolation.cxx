/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

/* Header include */
#include "FastCaloSim/Extrapolation/FastCaloSimCaloExtrapolation.h"

/* ISF includes */
#include "FastCaloSim/Core/MLogging.h"
#include "FastCaloSim/Core/TFCSExtrapolationState.h"
#include "FastCaloSim/Core/TFCSTruthState.h"
#include "FastCaloSim/Geometry/FastCaloSim_CaloCell_ID.h"
#include "FastCaloSim/Geometry/ICaloGeometry.h"

/* G4FieldTrack used to store transportation steps */
#include "G4FieldTrack.hh"

/* Preprocessor macro to use
   -- DEBUG   if CONDITION is True
   -- WARNING if CONDITION is False
 */
#define ATH_MSG_COND(MSG, CONDITION) \
  do { \
    if (CONDITION) { \
      ATH_MSG_DEBUG(MSG); \
    } else { \
      ATH_MSG_WARNING(MSG); \
    } \
  } while (0)

FastCaloSimCaloExtrapolation::FastCaloSimCaloExtrapolation()
    : m_geo(nullptr)
{
}

void FastCaloSimCaloExtrapolation::extrapolate(
    TFCSExtrapolationState& result,
    const TFCSTruthState* truth,
    const std::vector<G4FieldTrack>& caloSteps) const
{
  ATH_MSG_DEBUG("[extrapolate] Initializing extrapolation to ID-Calo boundary");
  extrapolateToID(result, caloSteps, truth);

  ATH_MSG_DEBUG(
      "[extrapolate] Initializing extrapolation to calorimeter layers");
  extrapolateToLayers(result, caloSteps, truth);

  ATH_MSG_DEBUG("[extrapolate] Extrapolation done");
}

void FastCaloSimCaloExtrapolation::extrapolateToID(
    TFCSExtrapolationState& result,
    const std::vector<G4FieldTrack>& caloSteps,
    const TFCSTruthState* truth) const
{
  ATH_MSG_DEBUG("Start extrapolateToID()");

  // pT threshold of truth particles over which extrapolation failures will be
  // printed as warnings
  const float transverseMomWarningLimit = 500;

  // initialize values
  result.set_IDCaloBoundary_eta(-999.);
  result.set_IDCaloBoundary_phi(-999.);
  result.set_IDCaloBoundary_r(0);
  result.set_IDCaloBoundary_z(0);
  result.set_IDCaloBoundary_AngleEta(-999.);
  result.set_IDCaloBoundary_Angle3D(-999.);

  // magnitude of extrapolated position
  double extPosDist = -1;

  for (unsigned int surfID = 0; surfID < 3; surfID++) {
    double R = m_CaloBoundaryR.at(surfID);
    double Z = m_CaloBoundaryZ.at(surfID);

    ATH_MSG_DEBUG("[ExtrapolateToID] Extrapolating to ID-Calo boundary with ID="
                  << surfID << " R=" << R << " Z=" << Z);

    // extrapolated position and momentum direction at IDCaloBoundary
    Vector3D extPos, momDir;

    // main extrapolation call
    if (!extrapolateToCylinder(caloSteps, R, Z, extPos, momDir))
      continue;

    double tolerance = 0.001;

    // test if z inside previous cylinder within some tolerance
    ATH_MSG_DEBUG(
        "[ExtrapolateToID] Testing condition 1: hit z=" << extPos.z());
    if (surfID > 0
        && std::abs(extPos.z()) < m_CaloBoundaryZ.at(surfID - 1) - tolerance)
      continue;
    ATH_MSG_DEBUG("[ExtrapolateToID] Passed condition 1.");

    // test if r inside next cylinder within some tolerance
    ATH_MSG_DEBUG(
        "[ExtrapolateToID] Testing condition 2: hit r=" << extPos.perp());
    if (surfID < m_CaloBoundaryR.size() - 1
        && extPos.perp() < m_CaloBoundaryR.at(surfID + 1) - tolerance)
      continue;
    ATH_MSG_DEBUG("[ExtrapolateToID] Passed condition 2.");

    ATH_MSG_DEBUG("[ExtrapolateToID] Testing condition 3: hit magnitude="
                  << extPos.mag());
    if (extPosDist >= 0 && extPos.mag() > extPosDist)
      continue;
    ATH_MSG_DEBUG("[ExtrapolateToID] Passed condition 3.");

    extPosDist = extPos.mag();

    result.set_IDCaloBoundary_eta(extPos.eta());
    result.set_IDCaloBoundary_phi(extPos.phi());
    result.set_IDCaloBoundary_r(extPos.perp());
    result.set_IDCaloBoundary_z(extPos.z());

    ATH_MSG_DEBUG("[ExtrapolateToID] Setting IDCaloBoundary to eta="
                  << extPos.eta() << " phi=" << extPos.phi()
                  << " r=" << extPos.perp() << " z=" << extPos.z());

    // compute angle between extrapolated position vector and momentum at
    // IDCaloBoundary can be used to correct shower shapes for particles which
    // do not originate from {0,0,0}
    double Angle3D = extPos.angle(momDir);
    double AngleEta = extPos.theta() - momDir.theta();
    result.set_IDCaloBoundary_AngleEta(AngleEta);
    result.set_IDCaloBoundary_Angle3D(Angle3D);

  }  // end of loop over surfaces

  if (result.IDCaloBoundary_eta() == -999)
    ATH_MSG_COND(
        "[ExtrapolateToID] Failed extrapolation to ID-Calo boundary. "
        "\n[ExtrapolateToID] Particle with truth vertex at ("
            << truth->vertex().X() << "," << truth->vertex().Y() << ","
            << truth->vertex().Z() << ")"
            << " with"
            << " PdgId=" << truth->pdgid() << " pT=" << truth->Pt()
            << " eta=" << truth->Eta() << " phi=" << truth->Phi()
            << " E=" << truth->E() << " Ekin_off=" << truth->Ekin_off(),
        truth->Pt() < transverseMomWarningLimit);

  ATH_MSG_DEBUG("[ExtrapolateToID] End extrapolateToID()");
}

void FastCaloSimCaloExtrapolation::extrapolateToLayers(
    TFCSExtrapolationState& result,
    const std::vector<G4FieldTrack>& caloSteps,
    const TFCSTruthState* truth) const
{
  ATH_MSG_DEBUG("[extrapolateToLayers] Start extrapolate");

  // pT threshold of truth particles over which extrapolation failures will be
  // printed as warnings
  const float transverseMomWarningLimit = 500;

  //////////////////////////////////////
  // Start calo extrapolation
  //////////////////////////////////////

  // only continue if inside the calo
  if (std::abs(result.IDCaloBoundary_eta()) < 6) {
    // now try to extrapolate to all calo layers that contain energy
    for (int sample = CaloCell_ID_FCS::FirstSample;
         sample < CaloCell_ID_FCS::MaxSample;
         ++sample)
    {
      for (int subpos = SUBPOS_MID; subpos <= SUBPOS_EXT; ++subpos) {
        float cylR, cylZ;
        if (m_geo->isCaloBarrel(sample)) {
          cylR = std::abs(
              m_geo->rpos(sample, result.IDCaloBoundary_eta(), subpos));
          // EMB0 - EMB3 use z position of EME1 front end surface for
          // extrapolation else extrapolate to cylinder with symmetrized maximum
          // Z bounds set eta to a dummy value of 1000 and -1000 to force
          // detector side
          if (sample < 4)
            cylZ = result.IDCaloBoundary_eta() > 0
                ? std::abs(m_geo->zpos(5, 1000, 1))
                : std::abs(m_geo->zpos(5, -1000, 1));
          else
            cylZ = 0.5
                * (std::abs(m_geo->zpos(sample, 1000, subpos))
                   + std::abs(m_geo->zpos(sample, -1000, subpos)));
        } else {
          // if we are not at barrel surface, extrapolate to cylinder with
          // maximum R to reduce extrapolation length
          cylZ = std::abs(
              m_geo->zpos(sample, result.IDCaloBoundary_eta(), subpos));
          // calculate radius of cylinder we will extrapolate to
          double mineta, maxeta, eta;
          m_geo->minmaxeta(sample, result.IDCaloBoundary_eta(), mineta, maxeta);
          // get eta where we will look up the layer radius
          eta = result.IDCaloBoundary_eta() > 0 ? mineta : maxeta;
          // calculate azimuthal angle from pseudorapidity
          double theta = 2 * std::atan(std::exp(-eta));
          // calculate maximum R of last cell of layer from z and theta
          cylR = std::abs(
              cylZ * std::sqrt((1 / (std::cos(theta) * std::cos(theta))) - 1));
        }

        Vector3D extPos, momDir;
        if (extrapolateToCylinder(caloSteps, cylR, cylZ, extPos, momDir)) {
          // scale the extrapolation to fit the radius of the cylinder in the
          // case of barrel and scale extrapolation to fit z component in case
          // of endcap layer scale is only non-unitary in case we extrapolate to
          // the endcaps of the cylinder for barrel and in case we extrapolate
          // to cover for endcaps this will keep phi, eta intact and only scale
          // r and z to fit a sensible position on the cylinder
          double scale = 1;
          if (m_geo->isCaloBarrel(sample) && std::abs(extPos.perp()) > 1e-6)
            scale = cylR / extPos.perp();
          else if (!m_geo->isCaloBarrel(sample) && std::abs(extPos.z()) > 1e-6)
            scale = cylZ / std::abs(extPos.z());
          // scale extrapolated position accordingly
          extPos = scale * extPos;

          result.set_OK(sample, subpos, true);
          result.set_phi(sample, subpos, extPos.phi());
          result.set_z(sample, subpos, extPos.z());
          result.set_eta(sample, subpos, extPos.eta());
          result.set_r(sample, subpos, extPos.perp());
        } else {
          ATH_MSG_COND(
              " [extrapolateToLayers] Extrapolation to cylinder failed. Sample="
                  << sample << " subpos=" << subpos
                  << " eta=" << result.IDCaloBoundary_eta()
                  << " phi=" << result.IDCaloBoundary_phi()
                  << " r=" << result.IDCaloBoundary_r()
                  << " z=" << result.IDCaloBoundary_z(),
              truth->Pt() < transverseMomWarningLimit);
        }
      }  // for pos
    }  // for sample
  }  // inside calo

  else
    ATH_MSG_COND(
        "[extrapolateToLayers] Ups. Not inside calo. "
        "result.IDCaloBoundary_eta()="
            << result.IDCaloBoundary_eta()
            << "\n[extrapolateToLayers] Particle with truth vertex at ("
            << truth->vertex().X() << "," << truth->vertex().Y() << ","
            << truth->vertex().Z() << ")"
            << " with"
            << " PdgId=" << truth->pdgid() << " pT=" << truth->Pt()
            << " eta=" << truth->Eta() << " phi=" << truth->Phi()
            << " E=" << truth->E() << " Ekin_off=" << truth->Ekin_off(),
        truth->Pt() < transverseMomWarningLimit);

  ATH_MSG_DEBUG("[extrapolateToLayers] End extrapolateToLayers()");
}

bool FastCaloSimCaloExtrapolation::extrapolateToCylinder(
    const std::vector<G4FieldTrack>& caloSteps,
    float cylR,
    float cylZ,
    Vector3D& extPos,
    Vector3D& momDir) const
{
  if (caloSteps.size() == 1) {
    Vector3D hitPos = caloSteps.at(0).GetPosition();
    ATH_MSG_DEBUG("[extrapolateWithPCA(R="
                  << cylR << ",Z=" << cylZ
                  << ")] Extrapolating single hit position to surface.");
    extPos = projectOnCylinder(cylR, cylZ, hitPos);
    momDir = caloSteps.at(0).GetMomentum();
    return true;
  }

  // if we do not find any good intersections, extrapolate to closest point on
  // surface
  bool foundHit =
      extrapolateWithIntersection(caloSteps, cylR, cylZ, extPos, momDir)
      ? true
      : extrapolateWithPCA(caloSteps, cylR, cylZ, extPos, momDir);

  if (foundHit) {
    ATH_MSG_DEBUG("[extrapolateToCylinder(R="
                  << cylR << ",Z=" << cylZ
                  << ")::END] Extrapolated to cylinder with R=" << cylR
                  << " and Z=" << cylZ << " at (" << extPos.x() << ","
                  << extPos.y() << "," << extPos.z() << ")");
  } else {
    // this is not expected to ever happen
    ATH_MSG_DEBUG("(R=" << cylR << ", Z=" << cylZ
                        << "::END) Extrapolation to cylinder surface failed!");
  }

  return foundHit;
}

bool FastCaloSimCaloExtrapolation::extrapolateWithIntersection(
    const std::vector<G4FieldTrack>& caloSteps,
    float cylR,
    float cylZ,
    Vector3D& extPos,
    Vector3D& momDir) const
{
  ATH_MSG_DEBUG("[extrapolateWithIntersection(R="
                << cylR << ",Z=" << cylZ
                << ")] Checking for cylinder intersections of line segments.");

  // counter for number of computed extrapolations, does not count cases of
  // rejected extrapolations due to close by hit positions
  unsigned int nExtrapolations = 0;
  for (size_t hitID = 1; hitID < caloSteps.size(); hitID++) {
    // initialize intersection result variables
    // get current and consecutive hit position and build hitLine
    Vector3D hitPos1 = caloSteps.at(hitID - 1).GetPosition();
    Vector3D hitPos2 = caloSteps.at(hitID).GetPosition();
    Vector3D hitDir = hitPos2 - hitPos1;

    ATH_MSG_DEBUG(
        "[extrapolateWithIntersection(R="
        << cylR << ",Z=" << cylZ << ")] Considering line segment between ("
        << hitPos1.x() << "," << hitPos1.y() << "," << hitPos1.z() << ") and ("
        << hitPos2.x() << "," << hitPos2.y() << "," << hitPos2.z() << ")");
    // get position of the hit positions on the cylinder
    HITPOSITION cylPosHit1 = whereOnCylinder(cylR, cylZ, hitPos1);
    HITPOSITION cylPosHit2 = whereOnCylinder(cylR, cylZ, hitPos2);

    // check if one of the hit positions already lays on the cylinder surface
    if (cylPosHit1 == ON || cylPosHit2 == ON) {
      extPos = cylPosHit1 == ON ? hitPos1 : hitPos2;
      momDir = cylPosHit1 == ON ? caloSteps.at(hitID - 1).GetMomentum()
                                : caloSteps.at(hitID).GetMomentum();
      ATH_MSG_DEBUG("[extrapolateWithIntersection(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Hit position already on cylinder surface.");
      return true;
    }

    // do not try to extrapolate with intersections if the hit position are very
    // close together
    if (hitDir.mag() < 0.01)
      continue;

    // get intersections through cylinder
    CylinderIntersections intersections =
        getCylinderIntersections(cylR, cylZ, hitPos1, hitPos2);
    nExtrapolations++;

    Vector3D selectedIntersection(0, 0, 0);

    // select the best intersection
    if (intersections.number == 1)
      selectedIntersection = intersections.first;
    else if (intersections.number > 1)
      selectedIntersection = whichIntersection(cylR,
                                               cylZ,
                                               hitPos1,
                                               hitPos2,
                                               intersections.first,
                                               intersections.second)
              == 0
          ? intersections.first
          : intersections.second;

    if (intersections.number > 0) {
      bool isForwardExtrapolation =
          (selectedIntersection.x() - hitPos1.x()) / (hitPos2.x() - hitPos1.x())
          >= 0;
      bool travelThroughSurface =
          doesTravelThroughSurface(cylR, cylZ, hitPos1, hitPos2);

      // do not allow for backward extrapolation except in the case of first two
      // (distinguishable) hit positions outside cylinder and in the case we
      // detect a travel though the surface
      if (nExtrapolations > 1 && !isForwardExtrapolation
          && !travelThroughSurface)
        continue;

      // check if the intersection between infinite line and cylinder lays on
      // segment spanned by hit positions
      bool intersectionOnSegment =
          isOnSegment(selectedIntersection, hitPos1, hitPos2);
      // check if both hit positions lay outside of the cylinder
      bool hitPosOutside = cylPosHit1 == OUTSIDE && cylPosHit2 == OUTSIDE;

      // we found our extrapolated hit position in case that either
      // we detect that the line segment crosses the surface of the cylinder
      // the intersection between the infinite lines and the cylinder lays on
      // the line segment both hit positions are outside of the cylinder and
      // there is a backwards extrapolation for the first two hit positions if
      // this is not the case for any of the hit position pairs we will use the
      // last two hit position for the linear extrapolation if these do not have
      // any intersection, then we will pass back to extrapolateWithPCA
      if (travelThroughSurface || intersectionOnSegment
          || (hitPosOutside && !isForwardExtrapolation && nExtrapolations == 1)
          || caloSteps.size() - 1 == hitID)
      {
        // take momentum direction of hit position closest to cylinder surface
        // alternatively one could also take the extrapolated direction normDir
        // = hitPos2 - hitPos1
        double distHitPos1 =
            (hitPos1 - projectOnCylinder(cylR, cylZ, hitPos1)).mag();
        double distHitPos2 =
            (hitPos2 - projectOnCylinder(cylR, cylZ, hitPos2)).mag();
        momDir = distHitPos1 < distHitPos2
            ? caloSteps.at(hitID - 1).GetMomentum()
            : caloSteps.at(hitID).GetMomentum();
        extPos = selectedIntersection;
        return true;
      }
      ATH_MSG_DEBUG("[extrapolateWithIntersection(R="
                    << cylR << ",Z=" << cylZ << ")] Extrapolated position at ("
                    << selectedIntersection.x() << ","
                    << selectedIntersection.y() << ","
                    << selectedIntersection.z() << ")");
    }
  }  // end of loop over hit positions

  return false;
}

bool FastCaloSimCaloExtrapolation::extrapolateWithPCA(
    const std::vector<G4FieldTrack>& caloSteps,
    float cylR,
    float cylZ,
    Vector3D& extPos,
    Vector3D& momDir) const
{
  bool foundHit = false;
  ATH_MSG_DEBUG("[extrapolateWithPCA(R="
                << cylR << ",Z=" << cylZ
                << ")] No forward intersections with cylinder surface. "
                   "Extrapolating to closest point on surface.");

  // here we also need to consider distances from line segments to the cylinder
  double minDistToSurface = 100000;
  for (size_t hitID = 1; hitID < caloSteps.size(); hitID++) {
    Vector3D hitPos1 = caloSteps.at(hitID - 1).GetPosition();
    Vector3D hitPos2 = caloSteps.at(hitID).GetPosition();

    ATH_MSG_DEBUG(
        "[extrapolateWithPCA(R="
        << cylR << ",Z=" << cylZ << ")] Considering line segment between ("
        << hitPos1.x() << "," << hitPos1.y() << "," << hitPos1.z() << ") and ("
        << hitPos2.x() << "," << hitPos2.y() << "," << hitPos2.z() << ")");

    Vector3D PCA;
    // find the point of closest approach (PCA) to the cylinder on the line
    // segment
    findPCA(cylR, cylZ, hitPos1, hitPos2, PCA);
    // compute distance between PCA and cylinder
    Vector3D cylinderSurfacePCA = projectOnCylinder(cylR, cylZ, PCA);
    double tmpMinDistToSurface = (PCA - cylinderSurfacePCA).mag();

    ATH_MSG_DEBUG(
        "[extrapolateWithPCA(R="
        << cylR << ",Z=" << cylZ << ")] Extrapolated line segment to ("
        << cylinderSurfacePCA.x() << "," << cylinderSurfacePCA.y() << ","
        << cylinderSurfacePCA.z() << ") with distance " << tmpMinDistToSurface);

    if (tmpMinDistToSurface < minDistToSurface) {
      foundHit = true;
      extPos = cylinderSurfacePCA;
      // take momentum direction of hit position closest to cylinder surface
      // alternatively one could also take the extrapolated direction normDir =
      // hitPos2 - hitPos1
      double distHitPos1 =
          (hitPos1 - projectOnCylinder(cylR, cylZ, hitPos1)).mag();
      double distHitPos2 =
          (hitPos2 - projectOnCylinder(cylR, cylZ, hitPos2)).mag();
      momDir = distHitPos1 < distHitPos2 ? caloSteps.at(hitID - 1).GetMomentum()
                                         : caloSteps.at(hitID).GetMomentum();

      minDistToSurface = tmpMinDistToSurface;
    }
  }  // end over loop of hit positions

  return foundHit;
}

void FastCaloSimCaloExtrapolation::findPCA(float cylR,
                                           float cylZ,
                                           Vector3D& hitPos1,
                                           Vector3D& hitPos2,
                                           Vector3D& PCA) const
{
  // in the following we will try to find the closest point-of-approach (PCA) to
  // the cylinder on the line segment hit direction
  Vector3D hitDir = hitPos2 - hitPos1;

  // project both hit positions onto the cylinder
  Vector3D projCylinderHitPos1 = projectOnCylinder(cylR, cylZ, hitPos1);
  Vector3D projCylinderHitPos2 = projectOnCylinder(cylR, cylZ, hitPos2);
  // direction of line spanned by the two projected points on the cylinder
  // surface
  Vector3D cylinderProjDir = projCylinderHitPos2 - projCylinderHitPos1;

  // CASE A: projections on the cylinder are close enough, take one of the hit
  // positions as PCA
  if (cylinderProjDir.mag() < 0.0001) {
    PCA = hitPos1;
    return;
  };

  // CASE B: we are outside the Z bounds of the cylinder
  if ((hitPos1.z() > cylZ || hitPos1.z() < -cylZ)
      || (hitPos2.z() > cylZ || hitPos2.z() < -cylZ))
  {
    // calculate PCA to point on endcap
    Vector3D cylZEndcap(0, 0, cylZ);
    bool isParallelToEndcap = std::abs(hitPos1.z() - hitPos2.z()) < 0.00001;

    // Check if parallel to endcap plane
    if (isParallelToEndcap) {
      // if both inside there are infinite solutions take one in the middle
      Vector3D intersectA(0, 0, 0);
      Vector3D intersectB(0, 0, 0);

      int nIntersections = circleLineIntersection2D(
          cylR, hitPos1, hitPos2, intersectA, intersectB);

      if (nIntersections == 2) {
        bool IntAOnSegment = isOnSegment(intersectA, hitPos1, hitPos2);
        bool IntBOnSegment = isOnSegment(intersectB, hitPos1, hitPos2);

        if (IntAOnSegment && IntBOnSegment)
          PCA = intersectA + 0.5 * (intersectB - intersectA);
        else if (IntAOnSegment)
          PCA = hitPos1.perp() <= cylR
              ? intersectA + 0.5 * (hitPos1 - intersectA)
              : intersectA + 0.5 * (hitPos2 - intersectA);
        else if (IntBOnSegment)
          PCA = hitPos1.perp() <= cylR
              ? intersectB + 0.5 * (hitPos1 - intersectB)
              : intersectB + 0.5 * (hitPos2 - intersectB);
        // intersections are not on line segment, i.e. line segment is within
        // extended cylinder
        else
          PCA = hitPos1 + 0.5 * hitDir;

      } else if (!isZero(intersectA) || !isZero(intersectB)) {
        // this can only happen if the extended line is tangential to the
        // cylinder if intersection lays on segment PCA will be intersection, if
        // not it will be the corresponding end points
        Vector3D intersect = isZero(intersectA) ? intersectB : intersectA;
        Vector3D hitPos =
            (hitPos1 - intersect).mag() < (hitPos2 - intersect).mag() ? hitPos1
                                                                      : hitPos2;
        bool IntOnSegment = isOnSegment(intersectA, hitPos1, hitPos2);
        PCA = IntOnSegment ? intersect : hitPos;

      } else {
        // line segment is outside extended cylinder
        // PCA corresponds to closest distance to center {0, 0, cylZ}
        Vector3D infLinePCA = hitPos1
            + ((cylZEndcap - hitPos1).dot(hitDir) / hitDir.dot(hitDir))
                * (hitDir);

        if (isOnSegment(infLinePCA, hitPos1, hitPos2))
          PCA = infLinePCA;
        else
          PCA = (hitPos1 - infLinePCA).mag() < (hitPos2 - infLinePCA).mag()
              ? hitPos1
              : hitPos2;
      }
    }

    else
    {
      // figure out all other cases iteratively beginning with BoundA and BoundB
      Vector3D BoundA, BoundB;
      // this is point on line closest to {0, 0, cylZ}, always on segment
      double t = ((cylZEndcap - hitPos1).dot(hitDir) / hitDir.dot(hitDir));
      BoundA = t <= 0 ? hitPos1 : (t >= 1 ? hitPos2 : hitPos1 + t * hitDir);

      // calculate intersection point of line segment and endcap plane and
      // project intersection onto cylinder check if t is between 0 and 1, if
      // not, take hitpos as starting bound
      t = (cylZ - hitPos1.z()) / hitDir.z();
      BoundB = t <= 0 ? hitPos1 : (t >= 1 ? hitPos2 : hitPos1 + t * hitDir);
      // looks for the PCA iteratively in cases there is no easy analytical
      // solution
      getIterativePCA(cylR, cylZ, BoundA, BoundB, PCA);
    }

    return;
  }

  // compute point of closest approach to z axis
  // this is analogous to finding the PCA of two 3D lines
  Vector3D PCACylBounds = get3DLinePCA(
      hitPos1, hitDir.unit(), Vector3D(0, 0, 0), Vector3D(0, 0, 1));

  // take PCA on line in case it lays on segment, otherwise take closest hit
  // position to surface

  double distSurfHit1 = (projCylinderHitPos1 - hitPos1).mag();
  double distSurfHit2 = (projCylinderHitPos2 - hitPos2).mag();

  PCA = isOnSegment(PCACylBounds, hitPos1, hitPos2)
      ? PCACylBounds
      : (distSurfHit1 < distSurfHit2 ? hitPos1 : hitPos2);
}

void FastCaloSimCaloExtrapolation::getIterativePCA(float cylR,
                                                   float cylZ,
                                                   Vector3D& BoundA,
                                                   Vector3D& BoundB,
                                                   Vector3D& PCA) const
{
  ATH_MSG_DEBUG("[getIterativePCA] Finding PCA iteratively.");

  Vector3D boundDir = BoundB - BoundA;
  double distBounds = boundDir.mag();

  // this sets the precision of the iterative finding procedure
  const double stepSize = 0.01;

  // if bounds are close enough together, there is nothing to do
  // should make sure that nHalfDivisions >= 2
  if (distBounds <= 4 * stepSize) {
    PCA = BoundA + 0.5 * (BoundB - BoundA);
    return;
  }

  Vector3D tmpBoundA, tmpBoundB, tmpOnCylinderBoundA, tmpOnCylinderBoundB;
  Vector3D resBoundA, resBoundB, resOnCylinderBoundA, resOnCylinderBoundB;

  // initial positions on cylinder and distance to line segment
  Vector3D OnCylinderBoundA = projectOnCylinder(cylR, cylZ, BoundA);
  Vector3D OnCylinderBoundB = projectOnCylinder(cylR, cylZ, BoundB);

  double minDistA = (BoundA - OnCylinderBoundA).mag();
  double minDistB = (BoundB - OnCylinderBoundB).mag();

  // initialize result bounds with closest input bounds as fall back option
  if (minDistA < minDistB) {
    resBoundA = BoundA;
    resBoundB = BoundA;
  } else {
    resBoundA = BoundB;
    resBoundB = BoundB;
  }
  double tmpMinDistA, tmpMinDistB;
  unsigned int nHalfDivisions = (distBounds / stepSize) / 2;
  for (unsigned int step = 0; step < nHalfDivisions; step++) {
    // temporary bounds on line segment
    tmpBoundA = BoundA + (step + 1) * stepSize * (boundDir / distBounds);
    tmpBoundB = BoundB - (step + 1) * stepSize * (boundDir / distBounds);

    // temporary projected bounds on cylinder
    tmpOnCylinderBoundA = projectOnCylinder(cylR, cylZ, tmpBoundA);
    tmpOnCylinderBoundB = projectOnCylinder(cylR, cylZ, tmpBoundB);

    // temporary minimum distance between bound on segment and bound on cylinder
    tmpMinDistA = (tmpBoundA - tmpOnCylinderBoundA).mag();
    tmpMinDistB = (tmpBoundB - tmpOnCylinderBoundB).mag();

    if (minDistA >= tmpMinDistA) {
      minDistA = tmpMinDistA;
    } else {
      double t = (step * stepSize) / distBounds;
      resBoundA = BoundA + t * boundDir;
      resBoundB = tmpBoundA;
      break;
    }

    if (minDistB >= tmpMinDistB) {
      minDistB = tmpMinDistB;
    } else {
      double t = (step * stepSize) / distBounds;
      resBoundB = BoundB - t * boundDir;
      resBoundA = tmpBoundB;
      break;
    }
  }

  // return middle of best bounds
  PCA = resBoundA + 0.5 * (resBoundB - resBoundA);
}

auto FastCaloSimCaloExtrapolation::circleLineIntersection2D(
    float circR,
    Vector3D& pointA,
    Vector3D& pointB,
    Vector3D& intersectA,
    Vector3D& intersectB) const -> int
{
  // find intersections intA and intB with line spanned by pointA and pointB
  // returns number of intersections
  // assumes circle lays in xy plane

  double dx, dy, A, B, C, det, t;

  dx = pointB.x() - pointA.x();
  dy = pointB.y() - pointA.y();

  A = dx * dx + dy * dy;
  B = 2 * (dx * pointA.x() + dy * pointA.y());
  C = pointA.x() * pointA.x() + pointA.y() * pointA.y() - circR * circR;

  det = B * B - 4 * A * C;

  if (A <= 0.0000001 || det < 0) {
    ATH_MSG_DEBUG("[circleLineIntersection2D] No intersections.");
    return 0;
  } else if (std::abs(det) < 0.00001) {
    // one solution, tangential case.
    t = -B / (2 * A);
    intersectA = {pointA.x() + t * dx, pointA.y() + t * dy, pointA.z()};
    ATH_MSG_DEBUG("[circleLineIntersection2D] One intersection at ("
                  << intersectA.x() << "," << intersectA.y() << ","
                  << intersectA.z() << ").");
    return 1;
  } else {
    // two solutions
    t = (-B + std::sqrt(det)) / (2 * A);
    intersectA = {pointA.x() + t * dx, pointA.y() + t * dy, pointA.z()};
    t = (-B - std::sqrt(det)) / (2 * A);
    intersectB = {pointA.x() + t * dx, pointA.y() + t * dy, pointB.z()};
    ATH_MSG_DEBUG("[circleLineIntersection2D] Two intersections at ("
                  << intersectA.x() << "," << intersectA.y() << ","
                  << intersectA.z() << ") and at (" << intersectB.x() << ","
                  << intersectB.y() << "," << intersectB.z() << ").");
    return 2;
  }
}

Vector3D FastCaloSimCaloExtrapolation::projectOnCylinder(float cylR,
                                                         float cylZ,
                                                         Vector3D& hitPos)
{
  Vector3D closestPointOnCylinder;
  Vector3D cylAxis(0, 0, cylZ);

  // positive side
  if (hitPos.z() >= cylZ) {
    // project hit position on x-y plane at positive side
    Vector3D projHitPos(hitPos.x(), hitPos.y(), cylZ);

    // if r of hit position outside cylinder, closest hit is always on edge
    if (hitPos.perp() > cylR)
      closestPointOnCylinder = cylAxis + cylR * (projHitPos - cylAxis).unit();
    else
      closestPointOnCylinder =
          cylAxis + hitPos.perp() * (projHitPos - cylAxis).unit();

  }
  // negative side
  else if (hitPos.z() <= -cylZ)
  {
    // project hit position on x-y plane at negative side
    Vector3D projHitPos(hitPos.x(), hitPos.y(), -cylZ);

    if (hitPos.perp() > cylR)
      closestPointOnCylinder = -cylAxis + cylR * (projHitPos + cylAxis).unit();
    else
      closestPointOnCylinder =
          -cylAxis + hitPos.perp() * (projHitPos + cylAxis).unit();

  } else {
    Vector3D hitPosZ(0, 0, hitPos.z());
    closestPointOnCylinder = hitPosZ + cylR * (hitPos - hitPosZ).unit();
  }

  return closestPointOnCylinder;
}

CylinderIntersections FastCaloSimCaloExtrapolation::getCylinderIntersections(
    float cylR, float cylZ, Vector3D& hitPos1, Vector3D& hitPos2) const
{
  // calculates intersection of infinite line with cylinder --> can have 0 or 2
  // intersections
  CylinderIntersections intersections;

  // look for intersections with the cover of the cylinder
  unsigned int nCoverIntersections = cylinderLineIntersection(
      cylR, cylZ, hitPos1, hitPos2, intersections.first, intersections.second);
  if (nCoverIntersections == 2) {
    ATH_MSG_DEBUG(
        "[getCylinderIntersections(R="
        << cylR << ",Z=" << cylZ
        << ")] Found two cylinder intersections through cylinder cover.");
    intersections.number = 2;
    return intersections;
  } else if (nCoverIntersections == 1) {
    Vector3D positiveEndcapIntersection, negativeEndcapIntersection;
    bool IsPositiveEndcapIntersection = cylinderEndcapIntersection(
        cylR, cylZ, true, hitPos1, hitPos2, positiveEndcapIntersection);
    bool IsNegativeEndcapIntersection = cylinderEndcapIntersection(
        cylR, cylZ, false, hitPos1, hitPos2, negativeEndcapIntersection);

    if (IsPositiveEndcapIntersection && IsNegativeEndcapIntersection) {
      // if we have a cover intersection we only expect one additional endcap
      // intersection both endcap intersections can be valid in case the
      // intersection is at the edge of the cylinder cover and endcap in that
      // case take the endcap intersection which is further away from the
      // cylinder cover intersection to prevent taking the same intersection
      // twice
      ATH_MSG_DEBUG(
          "[getCylinderIntersections(R="
          << cylR << ",Z=" << cylZ
          << ")] Found intersection through cylinder cover and both endcaps. "
             "Intersection seems to be at edge of cover and endcap.");
      intersections.second =
          (positiveEndcapIntersection - intersections.first).mag()
              > (negativeEndcapIntersection - intersections.first).mag()
          ? positiveEndcapIntersection
          : negativeEndcapIntersection;
      intersections.number = 2;
    } else if (IsPositiveEndcapIntersection) {
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found intersection through cylinder cover and "
                       "positive endcap.");
      intersections.second = positiveEndcapIntersection;
      intersections.number = 2;
    } else if (IsNegativeEndcapIntersection) {
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found intersection through cylinder cover and "
                       "negative endcap.");
      intersections.second = negativeEndcapIntersection;
      intersections.number = 2;
    } else {
      // line is tangential to cylinder cover
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found single intersection through cylinder cover.");
      intersections.number = 1;
    }

  } else {
    // no cylinder cover intersections
    Vector3D positiveEndcapIntersection, negativeEndcapIntersection;
    bool IsPositiveEndcapIntersection = cylinderEndcapIntersection(
        cylR, cylZ, true, hitPos1, hitPos2, positiveEndcapIntersection);
    bool IsNegativeEndcapIntersection = cylinderEndcapIntersection(
        cylR, cylZ, false, hitPos1, hitPos2, negativeEndcapIntersection);

    if (IsPositiveEndcapIntersection && IsNegativeEndcapIntersection) {
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found intersections through both endcaps.");
      intersections.first = positiveEndcapIntersection;
      intersections.second = negativeEndcapIntersection;
      intersections.number = 2;
    } else if (IsPositiveEndcapIntersection) {
      // don't expect this to ever happen
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found single intersection through positive endcap. "
                       "This should not happen.");
      intersections.first = positiveEndcapIntersection;
      intersections.number = 1;
    } else if (IsNegativeEndcapIntersection) {
      // don't expect this to ever happen
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found single intersection through negative endcap. "
                       "This should not happen.");
      intersections.first = negativeEndcapIntersection;
      intersections.number = 1;
    } else {
      ATH_MSG_DEBUG("[getCylinderIntersections(R="
                    << cylR << ",Z=" << cylZ
                    << ")] Found no cylinder intersections.");
      // no intersections at all
      intersections.number = 0;
    }
  }

  return intersections;
}

auto FastCaloSimCaloExtrapolation::get3DLinePCA(const Vector3D& posA,
                                                const Vector3D& dirA,
                                                const Vector3D& posB,
                                                const Vector3D& dirB) const
    -> Vector3D
{
  {
    // line A
    const Vector3D& ma = posA;
    const Vector3D& ea = dirA.unit();
    // line B
    const Vector3D& mb = posB;
    const Vector3D& eb = dirB.unit();
    // solve for the closest approach
    Vector3D mab(mb - ma);
    double eaTeb = ea.dot(eb);
    double denom = 1 - eaTeb * eaTeb;
    if (std::abs(denom) > 10e-7) {
      double lambda0 = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
      // evaluate in terms of direction
      // evaluate validaty in terms of bounds
      Vector3D result = (ma + lambda0 * ea);
      // return the result
      return result;
    }
    return posA;
  }
}

// calculates the intersection between the line defined by pointA and pointB and
// the cylinder cover defined by cylR and cylZ
int FastCaloSimCaloExtrapolation::cylinderLineIntersection(
    float cylR,
    float cylZ,
    Vector3D& pointA,
    Vector3D& pointB,
    Vector3D& intersectA,
    Vector3D& intersectB) const
{
  // projections of points spanning the line onto the xy plane
  Vector3D projPointA(pointA.x(), pointA.y(), 0);
  Vector3D projPointB(pointB.x(), pointB.y(), 0);
  Vector3D projDiff = projPointA - projPointB;

  // calculate distance from (0,0,0) to line spanned by projPointA and
  // projPointB
  double projDiffNorm2 = projDiff.dot(projDiff);
  double t = projPointA.dot(projDiff) / projDiffNorm2;
  double d2 = projPointA.dot(projPointA) - t * t * projDiffNorm2;

  if (d2 < 0) {
    ATH_MSG_COND("[cylinderLineIntersection] Got negative distance (d2="
                     << d2 << "). Forcing to 0.",
                 d2 > -0.001);
    d2 = 0;
  }

  // if distance larger than cylinder radius then there are no intersection and
  // we are done
  if (d2 > cylR * cylR)
    return 0;

  double k = std::sqrt((cylR * cylR - d2) / projDiffNorm2);

  intersectA = pointA + (t + k) * (pointB - pointA);
  intersectB = pointA + (t - k) * (pointB - pointA);

  // check if intersection is outside z bounds
  bool IntAisValid = (intersectA.z() <= cylZ && intersectA.z() >= -cylZ);
  bool IntBisValid = (intersectB.z() <= cylZ && intersectB.z() >= -cylZ);

  if (IntAisValid && IntBisValid)
    return 2;
  else if (IntAisValid)
    return 1;
  else if (IntBisValid) {
    intersectA = intersectB;
    return 1;
  }

  return 0;
}

bool FastCaloSimCaloExtrapolation::cylinderEndcapIntersection(
    float cylR,
    float cylZ,
    bool positiveEndcap,
    Vector3D& pointA,
    Vector3D& pointB,
    Vector3D& intersection)
{
  // normal and point on endcap defines the plane
  Vector3D pointOnEndcap;
  Vector3D normal(0, 0, 1);
  positiveEndcap ? pointOnEndcap = {0, 0, cylZ} : pointOnEndcap = {0, 0, -cylZ};
  Vector3D hitDir = (pointB - pointA);

  double denom = normal.dot(hitDir);
  if (std::abs(denom) > 1e-6) {
    double t = normal.dot(pointOnEndcap - pointB) / denom;
    // compute intersection regardless of direction (t>0 or t<0)
    intersection = pointB + t * hitDir;
    Vector3D v = intersection - pointOnEndcap;

    // check if intersection is within cylR bounds
    return std::sqrt(v.dot(v)) <= cylR;
  }

  return false;
}

int FastCaloSimCaloExtrapolation::whichIntersection(
    float cylR,
    float cylZ,
    Vector3D& hitPos1,
    Vector3D& hitPos2,
    Vector3D& intersectionA,
    Vector3D intersectionB) const
{
  // check if the hit positions are outside or inside the cylinder surface
  HITPOSITION cylPosHit1 = whereOnCylinder(cylR, cylZ, hitPos1);
  HITPOSITION cylPosHit2 = whereOnCylinder(cylR, cylZ, hitPos2);

  if ((cylPosHit1 == INSIDE) ^ (cylPosHit2 == INSIDE)) {
    /* CASE A: one hit position inside and one outside of the cylinder (travel
    through surface), one intersection is on cylinder, take intersection closest
    to line segment */
    ATH_MSG_DEBUG("[whichIntersection] Travel through surface.");
    return getPointLineSegmentDistance(intersectionA, hitPos1, hitPos2)
        > getPointLineSegmentDistance(intersectionB, hitPos1, hitPos2);
  } else if (cylPosHit1 == INSIDE && cylPosHit2 == INSIDE) {
    /* CASE B: both hit position inside, take intersection which points towards
     * travel direction of particle */
    Vector3D dirA = intersectionA - hitPos2;
    Vector3D dirB = intersectionB - hitPos2;
    Vector3D hitDir = hitPos2 - hitPos1;
    ATH_MSG_DEBUG("[whichIntersection] Both hit positions inside.");
    return dirA.dot(hitDir) < dirB.dot(hitDir);
  } else {
    // /* CASE C: both hit position outside and the intersections lay on the
    // segment, take intersection closest to second hit position */
    // /* CASE D: both hit positions are outside and the intersections are not
    // on the line segment, take intersection closest to one of the hit
    // positions */
    double distHitPosIntersectA = (hitPos2 - intersectionA).mag();
    double distHitPosIntersectB = (hitPos2 - intersectionB).mag();
    ATH_MSG_DEBUG("[whichIntersection] Both hit positions outside.");
    return distHitPosIntersectA > distHitPosIntersectB;
  }
}

double FastCaloSimCaloExtrapolation::getPointLineSegmentDistance(
    Vector3D& point, Vector3D& hitPos1, Vector3D& hitPos2)
{
  Vector3D hitDir = hitPos2 - hitPos1;
  Vector3D w = point - hitPos1;

  double c1 = w.dot(hitDir);
  if (c1 <= 0)
    return distance(point, hitPos1);
  double c2 = hitDir.dot(hitDir);
  if (c2 <= c1)
    return distance(point, hitPos2);
  double t = c1 / c2;
  Vector3D vec = hitPos1 + t * hitDir;
  return distance(point, vec);
}

enum FastCaloSimCaloExtrapolation::HITPOSITION
FastCaloSimCaloExtrapolation::whereOnCylinder(float cylR,
                                              float cylZ,
                                              Vector3D& hitPos)
{
  // set a 1mm tolerance within which the hit position is considered to be on
  // the cylinder surface setting this higher can lead to extrapolation failures
  // around truth particle eta ~4
  float tolerance = 1;

  bool isOnEndcap = hitPos.perp() <= cylR + tolerance
      && (hitPos.z() > 0 ? std::abs(hitPos.z() - cylZ) < tolerance
                         : std::abs(hitPos.z() + cylZ) < tolerance);
  bool isOnCover = std::abs(hitPos.perp() - cylR) < tolerance
      && hitPos.z() < cylZ && hitPos.z() > -cylZ;

  // check if hit position is on endcap or cover of cylinder
  if (isOnEndcap || isOnCover)
    return HITPOSITION::ON;

  // check if hit position is inside cover
  if (hitPos.z() < cylZ && hitPos.z() > -cylZ && hitPos.perp() < cylR)
    return HITPOSITION::INSIDE;

  return HITPOSITION::OUTSIDE;
}

bool FastCaloSimCaloExtrapolation::doesTravelThroughSurface(float cylR,
                                                            float cylZ,
                                                            Vector3D& hitPos1,
                                                            Vector3D& hitPos2)
{
  // travel through surface in case one hit position is outside and the other
  // outside of cylinder surface
  return (whereOnCylinder(cylR, cylZ, hitPos1) == INSIDE)
      ^ (whereOnCylinder(cylR, cylZ, hitPos2) == INSIDE);
}

bool FastCaloSimCaloExtrapolation::isOnSegment(Vector3D& point,
                                               Vector3D& hitPos1,
                                               Vector3D& hitPos2)
{
  return getPointLineSegmentDistance(point, hitPos1, hitPos2) < 0.001;
}