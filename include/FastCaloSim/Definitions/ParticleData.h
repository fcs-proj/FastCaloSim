// Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <vector>
/** Implementation of classification functions according to PDG2022.
 *  https://pdg.lbl.gov/2023/reviews/rpp2022-rev-monte-carlo-numbering.pdf
 */

class DecodedPID : public std::pair<int, std::vector<int>>
{
public:
  inline DecodedPID(const int& p)
  {
    this->first = p;
    this->second.reserve(10);
    int ap = std::abs(p);
    for (; ap; ap /= 10)
      this->second.push_back(ap % 10);
    std::reverse(this->second.begin(), this->second.end());
  }
  inline DecodedPID shift(const size_t n) const
  {
    return DecodedPID(this->first % int(std::pow(10, ndigits() - n)));
  }
  inline const int& operator()(const size_t n) const
  {
    return this->second.at(n);
  }
  inline const int& last() const { return this->second.back(); }
  inline const int& pid() const { return this->first; }
  inline int max_digit(const int m, const int n) const
  {
    return *std::max_element(second.rbegin() + m, second.rbegin() + n);
  }
  inline int min_digit(const int m, const int n) const
  {
    return *std::min_element(second.rbegin() + m, second.rbegin() + n);
  }
  inline size_t ndigits() const { return this->second.size(); }
};

namespace ParticleData
{

static const int TABLESIZE = 100;
static const std::array<int, TABLESIZE> triple_charge = {
    +0, -1, +2, -1, +2, -1, +2, -1, +2, +0, +0, -3, +0, -3, +0, -3, +0,
    -3, +0, +0, +0, +0, +0, +0, +3, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +3, +0, +0, +3, +0, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0,
    +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0};

static const int UQUARK = 1;
static const int DQUARK = 2;
static const int SQUARK = 3;
static const int CQUARK = 4;
static const int BQUARK = 5;
static const int TQUARK = 6;

static const int ELECTRON = 11;
static const int POSITRON = -ELECTRON;
static const int NU_E = 12;
static const int MUON = 13;
static const int NU_MU = 14;
static const int TAU = 15;
static const int NU_TAU = 16;

static const int GLUON = 21;
static const int PHOTON = 22;
static const int Z0BOSON = 23;
static const int WPLUSBOSON = 24;
static const int HIGGSBOSON = 25;
static const int GRAVITON = 39;
static const int LEPTOQUARK = 42;

/// PDG Ids for Mavtop madgraph UFO model found under DarkX. The
/// mavtop is a vector-like top partner with coupling to a dark photon.
/// Theory paper: https://arxiv.org/abs/1904.05893
/// Pheno paper: https://arxiv.org/pdf/2112.08425
static const int DARKPHOTON = 60000;
static const int MAVTOP = 60001;

static const int PIPLUS = 211;
static const int PIMINUS = -PIPLUS;
static const int PI0 = 111;
static const int K0L = 130;

static const int K0S = 310;
static const int K0 = 311;
static const int KPLUS = 321;
static const int DPLUS = 411;
static const int DSTAR = 413;
static const int D0 = 421;
static const int DSPLUS = 431;
static const int JPSI = 443;
static const int B0 = 511;
static const int BCPLUS = 541;
static const int PROTON = 2212;
static const int LAMBDA0 = 3122;
static const int LAMBDACPLUS = 4122;
static const int LAMBDAB0 = 5122;
static const int PSI2S = 20443;
/// PDG rule 8:
/// The pomeron and odderon trajectories and a generic reggeon trajectory
/// of states in QCD areassigned codes 990, 9990, and 110 respectively
static const int POMERON = 990;
static const int ODDERON = 9990;
static const int REGGEON = 110;

/// PDG rule 10:
/// Codes 81–100 are reserved for generator-specific pseudoparticles and
/// concepts. Codes 901–930, 1901–1930, 2901–2930, and 3901–3930 are for
/// additional components of Standard Modelparton distribution functions, where
/// the latter three ranges are intended to distinguish left/right/ longitudinal
/// components. Codes 998 and 999 are reserved for GEANT tracking pur-poses.
static const int GEANTINOPLUS = 998;
static const int GEANTINO0 = 999;

/// PDG rule 2:
/// Quarks and leptons are numbered consecutively starting from 1 and 11
/// respectively; to dothis they are first ordered by family and within
/// families by weak isospin.
/// APID: the fourth generation quarks are quarks.
template<class T>
inline bool isQuark(const T& p)
{
  return isQuark(p->pdg_id());
}
template<>
inline bool isQuark(const int& p)
{
  return p != 0 && (std::abs(p) <= 8 || std::abs(p) == MAVTOP);
}
template<>
inline bool isQuark(const DecodedPID& p)
{
  return isQuark(p.pid());
}

template<class T>
inline bool isStrange(const T& p)
{
  return isStrange(p->pdg_id());
}
template<>
inline bool isStrange(const int& p)
{
  return std::abs(p) == 3;
}

template<class T>
inline bool isCharm(const T& p)
{
  return isCharm(p->pdg_id());
}
template<>
inline bool isCharm(const int& p)
{
  return std::abs(p) == 4;
}

template<class T>
inline bool isBottom(const T& p)
{
  return isBottom(p->pdg_id());
}
template<>
inline bool isBottom(const int& p)
{
  return std::abs(p) == 5;
}

template<class T>
inline bool isTop(const T& p)
{
  return isTop(p->pdg_id());
}
template<>
inline bool isTop(const int& p)
{
  return std::abs(p) == 6;
}

/// APID: the fourth generation leptons are leptons.
template<class T>
inline bool isLepton(const T& p)
{
  return isLepton(p->pdg_id());
}
template<>
inline bool isLepton(const int& p)
{
  auto sp = std::abs(p);
  return sp >= 11 && sp <= 18;
}
template<>
inline bool isLepton(const DecodedPID& p)
{
  return isLepton(p.pid());
}

template<class T>
inline bool isSMLepton(const T& p)
{
  return isSMLepton(p->pdg_id());
}
template<>
inline bool isSMLepton(const int& p)
{
  auto sp = std::abs(p);
  return sp >= 11 && sp <= 16;
}
template<>
inline bool isSMLepton(const DecodedPID& p)
{
  return isSMLepton(p.pid());
}

/// APID: the fourth generation leptons are leptons.
template<class T>
inline bool isChLepton(const T& p)
{
  return isChLepton(p->pdg_id());
}
template<>
inline bool isChLepton(const int& p)
{
  auto sp = std::abs(p);
  return sp >= 11 && sp <= 18 && sp % 2 == 1;
}

template<class T>
inline bool isElectron(const T& p)
{
  return isElectron(p->pdg_id());
}
template<>
inline bool isElectron(const int& p)
{
  return std::abs(p) == ELECTRON;
}

template<class T>
inline bool isMuon(const T& p)
{
  return isMuon(p->pdg_id());
}
template<>
inline bool isMuon(const int& p)
{
  return std::abs(p) == MUON;
}

template<class T>
inline bool isTau(const T& p)
{
  return isTau(p->pdg_id());
}
template<>
inline bool isTau(const int& p)
{
  return std::abs(p) == TAU;
}

/// APID: the fourth generation neutrinos are neutrinos.
template<class T>
inline bool isNeutrino(const T& p)
{
  return isNeutrino(p->pdg_id());
}
template<>
inline bool isNeutrino(const int& p)
{
  auto sp = std::abs(p);
  return sp == 12 || sp == 14 || sp == 16 || sp == 18;
}

template<class T>
inline bool isSMNeutrino(const T& p)
{
  return isSMNeutrino(p->pdg_id());
}
template<>
inline bool isSMNeutrino(const int& p)
{
  auto sp = std::abs(p);
  return sp == 12 || sp == 14 || sp == 16;
}

template<class T>
inline bool isGluon(const T& p)
{
  return isGluon(p->pdg_id());
}
template<>
inline bool isGluon(const int& p)
{
  return p == GLUON;
}

template<class T>
inline bool isPhoton(const T& p)
{
  return isPhoton(p->pdg_id());
}
template<>
inline bool isPhoton(const int& p)
{
  return p == PHOTON;
}

template<class T>
inline bool isZ(const T& p)
{
  return isZ(p->pdg_id());
}
template<>
inline bool isZ(const int& p)
{
  return p == Z0BOSON;
}

template<class T>
inline bool isW(const T& p)
{
  return isW(p->pdg_id());
}
template<>
inline bool isW(const int& p)
{
  return std::abs(p) == WPLUSBOSON;
}

/// PDG rule 11j:
/// The nature of Dark Matter (DM) is not known, and therefore a definitive
/// classificationis too early. Candidates within specific scenarios are
/// classified therein, such as 1000022 for the lightest neutralino.
/// Generic fundamental states can be given temporary codes in the range 51 -
/// 60, with 51, 52 and 53 reserved for spin 0, 1/2 and 1 ones (this could also
/// be an axion state). Generic mediators of s-channel DM pair creation of
/// annihilationcan be given codes 54 and 55 for spin 0 or 1 ones. Separate
/// antiparticles, with negativecodes, may or may not exist. More elaborate new
/// scenarios should be constructed with n= 5 and nr = 9. APID: Only the 51-60
/// range is considered DM. The antiparticles are assumed to be existing.
template<class T>
inline bool isDM(const T& p)
{
  return isDM(p->pdg_id());
}
template<>
inline bool isDM(const int& p)
{
  auto sp = std::abs(p);
  return (sp >= 51 && sp <= 60) || sp == DARKPHOTON;
}

/// PDG rule 8:
/// The pomeron and odderon trajectories and a generic reggeon trajectory
/// of states in QCD areassigned codes 990, 9990, and 110 respectively
template<class T>
inline bool isTrajectory(const T& p)
{
  return isTrajectory(p->pdg_id());
}
template<>
inline bool isTrajectory(const int& p)
{
  return std::abs(p) == POMERON || std::abs(p) == ODDERON
      || std::abs(p) == REGGEON;
}

/// APID: HIGGS boson is only one particle.
template<class T>
inline bool isHiggs(const T& p)
{
  return isHiggs(p->pdg_id());
}
template<>
inline bool isHiggs(const int& p)
{
  return p == HIGGSBOSON;
}

template<class T>
inline bool isResonance(const T& p)
{
  return isZ(p) || isW(p) || isHiggs(p) || isTop(p);
}

template<class T>
inline bool isGraviton(const T& p)
{
  return isGraviton(p->pdg_id());
}
template<>
inline bool isGraviton(const int& p)
{
  return p == 39;
}

template<class T>
inline bool isLeptoQuark(const T& p)
{
  return isLeptoQuark(p->pdg_id());
}
template<>
inline bool isLeptoQuark(const int& p)
{
  return std::abs(p) == LEPTOQUARK;
}

template<class T>
inline bool isSUSY(const T& p)
{
  return isSUSY(p->pdg_id());
}
template<class T>
inline bool isTechnicolor(const T& p)
{
  return isTechnicolor(p->pdg_id());
}
template<class T>
inline bool isExcited(const T& p)
{
  return isExcited(p->pdg_id());
}
template<class T>
inline bool isKK(const T& p)
{
  return isKK(p->pdg_id());
}
template<class T>
inline bool isHiddenValley(const T& p)
{
  return isHiddenValley(p->pdg_id());
}
template<class T>
inline bool isDiquark(const T& p)
{
  return isDiquark(p->pdg_id());
}
template<class T>
inline bool isHadron(const T& p)
{
  return isHadron(p->pdg_id());
}
template<class T>
inline bool isMeson(const T& p)
{
  return isMeson(p->pdg_id());
}
template<class T>
inline bool isBaryon(const T& p)
{
  return isBaryon(p->pdg_id());
}
template<class T>
inline bool isTetraquark(const T& p)
{
  return isTetraquark(p->pdg_id());
}
template<class T>
inline bool isPentaquark(const T& p)
{
  return isPentaquark(p->pdg_id());
}
template<class T>
inline bool isNucleus(const T& p)
{
  return isNucleus(p->pdg_id());
}
template<class T>
inline bool isBSM(const T& p)
{
  return isBSM(p->pdg_id());
}
template<class T>
inline bool isValid(const T& p)
{
  return isValid(p->pdg_id());
}
template<>
inline bool isValid(const DecodedPID& p);
template<class T>
inline bool isTransportable(const T& p)
{
  return isTransportable(p->pdg_id());
}
template<class T>
inline bool isGenSpecific(const T& p)
{
  return isGenSpecific(p->pdg_id());
}
template<class T>
inline bool isGeantino(const T& p)
{
  return isGeantino(p->pdg_id());
}
template<class T>
inline bool isPythia8Specific(const T& p)
{
  return isPythia8Specific(p->pdg_id());
}
template<>
inline bool isPythia8Specific(const DecodedPID& p)
{
  return (p(0) == 9 && p(1) == 9 && p.ndigits() == 7);
}
/// Main Table
/// for MC internal use 81–100,901–930,998-999,1901–1930,2901–2930, and
/// 3901–3930
template<>
inline bool isGenSpecific(const int& p)
{
  if (p >= 81 && p <= 100)
    return true;
  if (p >= 901 && p <= 930)
    return true;
  if (p >= 998 && p <= 999)
    return true;
  if (p >= 1901 && p <= 1930)
    return true;
  if (p >= 2901 && p <= 2930)
    return true;
  if (p >= 3901 && p <= 3930)
    return true;
  return false;
}

template<>
inline bool isGeantino(const int& p)
{
  return (std::abs(p) == GEANTINO0 || std::abs(p) == GEANTINOPLUS);
}

/// PDG rule 11d
/// Fundamental supersymmetric particles are identified by adding a nonzeronto
/// the par-ticle number. T he superpartner of a boson or a left-handed fermion
/// hasn= 1whilethe superpartner of a right-handed fermion hasn= 2. When mixing
/// occurs, such asbetween the winos and charged Higgsinos to give charginos, or
/// between left and rightsfermions, the lighter physical state is given the
/// smaller basis state number.
template<>
inline bool isSUSY(const DecodedPID& p)
{
  return (p.ndigits() == 7 && (p(0) == 1 || p(0) == 2) && isValid(p.shift(2))
          && !isGenSpecific(p.shift(2).pid()));
}

/// PDG rule 11e
/// Technicolor states have n= 3, with technifermions treated like ordinary
/// fermions. States which are ordinary color singlets have nr= 0. Color octets
/// have nr= 1. If a state has non-trivial quantum numbers under the topcolor
/// groups SU(3)1×SU(3)2, the quantum numbers are specified by tech, ij, where i
/// and j are 1 or 2. nLis then 2i+j. The coloron V8, is a heavy gluon color
/// octet and thus is 3100021
template<>
inline bool isTechnicolor(const DecodedPID& p)
{
  return (p.ndigits() == 7 && p(0) == 3 && (p(1) == 0 || p(0) == 1)
          && isValid(p.shift(2)) && !isGenSpecific(p.shift(2).pid()));
}
template<>
inline bool isTechnicolor(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isTechnicolor(value_digits);
}

/// PDG rule 11f
/// Excited (composite) quarks and leptons are identified by setting n= 4 and
/// nr= 0
template<>
inline bool isExcited(const DecodedPID& p)
{
  return (p.ndigits() == 7 && (p(0) == 4 && p(1) == 0)
          && (isLepton(p.shift(2)) || isQuark(p.shift(2))));
}
template<>
inline bool isExcited(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isExcited(value_digits);
}

/// PDG rule 11h
/// A black hole in models with extra dimensions has code 5000040. Kaluza-Klein
/// excitations in models with extra dimensions have n= 5 or n= 6, to
/// distinguish excitations of left-or right-handed fermions or, in case of
/// mixing, the lighter or heavier state (cf. 11d). The non zero nr digit gives
/// the radial excitation number, in scenarios where the level spacing allows
/// these to be
///  distinguished. Should the model also contain supersymmetry, excited SUSY
///  states would be denoted by a nnr>0, with n= 1 or 2 as usual.
/// Should some colored states be long-lived enough that hadrons would form
/// around them, the coding strategy of 11g applies, with the initial two nnr
/// digits preserved in the combined code.
template<>
inline bool isKK(const DecodedPID& p)
{
  return (p.ndigits() == 7 && (p(0) == 5 || p(0) == 6));
}
template<>
inline bool isKK(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isExcited(value_digits);
}

/// PDG rule 11k
/// Hidden Valley particles have n= 4 and nr= 9, and trailing numbers in
/// agreement with their nearest-analog standard particles, as far as possible.
/// Thus 4900021 is the gauge boson gv of a confining gauge field, 490000 nqv
/// and 490001 nlv fundamental constituents charged or not under this, 4900022
/// is the γv of a non-confining field, and 4900 nqv1 nqv2 nJ a Hidden Valley
/// meson.
template<>
inline bool isHiddenValley(const DecodedPID& p)
{
  return (p.ndigits() == 7 && p(0) == 4 && p(1) == 9 && isValid(p.shift(2)));
}
template<>
inline bool isHiddenValley(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isHiddenValley(value_digits);
}

/// PDG rule 4
/// Diquarks have 4-digit numbers with nq1 >= nq2 and nq3 = 0
/// APID: the diquarks with fourth generation are not diquarks
template<>
inline bool isDiquark(const DecodedPID& p)
{
  if (p.ndigits() == 4 && p(0) >= p(1) && p(2) == 0 && p.last() % 2 == 1
      && p.max_digit(1, 3) <= 6)
    return true;
  return false;
}

/// Table 43.1
///  PDG rule 5a:
///  The numbers specifying the meson’s quark content conform to the convention
///  nq1= 0 and nq2 >= nq3. The special case K0L is the sole exception to this
///  rule. PDG rule 5C: The special numbers 310 and 130 are given to the K0S and
///  K0L respectively. APID: The special code K0 is used when a generator uses
///  K0S/K0L
template<>
inline bool isMeson(const DecodedPID& p)
{
  if (p.ndigits() < 3)
    return false;
  if (std::abs(p.pid()) == K0S)
    return true;
  if (std::abs(p.pid()) == K0L)
    return true;
  if (std::abs(p.pid()) == K0)
    return true;
  if (p.last() % 2 != 1)
    return false;
  if (p.max_digit(1, 3) >= 6)
    return false;
  if (p.max_digit(1, 3) == 0)
    return false;

  if (p.ndigits() == 3 && p(0) == p(1) && p.pid() < 0)
    return false;
  if (p.ndigits() == 5 && p(2) == p(3) && p.pid() < 0)
    return false;
  if (p.ndigits() == 7 && p(4) == p(5) && p.pid() < 0)
    return false;

  if (p.ndigits() == 3 && p(0) >= p(1) && p(1) != 0)
    return true;
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 1 && p(1) == 0)
    return true;
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 2 && p(1) == 0
      && p.last() > 1)
    return true;
  if (p.ndigits() == 5 && p(2) >= p(3) && p(3) != 0 && p(0) == 3 && p(1) == 0
      && p.last() > 1)
    return true;

  if (p.ndigits() == 6 && p(3) >= p(4) && p(4) != 0 && p.last() % 2 == 1)
    return true;

  if (p.ndigits() == 7 && p(4) >= p(5) && p(5) != 0)
    return true;

  return false;
}
/// Table 43.2
template<>
inline bool isBaryon(const DecodedPID& p)
{
  if (p.ndigits() < 4)
    return false;
  if (p.max_digit(1, 4) >= 6)
    return false;
  if (p.min_digit(1, 4) == 0)
    return false;
  if (p.ndigits() == 4
      && (p.last() == 2 || p.last() == 4 || p.last() == 6 || p.last() == 8))
    return true;

  if (p.ndigits() == 5 && p(0) == 1 && (p.last() == 2 || p.last() == 4))
    return true;
  if (p.ndigits() == 5 && p(0) == 3 && (p.last() == 2 || p.last() == 4))
    return true;

  if (p.ndigits() == 6) {
    if (p(0) == 1 && p(1) == 0 && p.last() == 2)
      return true;
    if (p(0) == 1 && p(1) == 1 && p.last() == 2)
      return true;
    if (p(0) == 1 && p(1) == 2 && p.last() == 4)
      return true;

    if (p(0) == 2 && p(1) == 0 && p.last() == 2)
      return true;
    if (p(0) == 2 && p(1) == 0 && p.last() == 4)
      return true;
    if (p(0) == 2 && p(1) == 1 && p.last() == 2)
      return true;

    if (p(0) == 1 && p(1) == 0 && p.last() == 4)
      return true;
    if (p(0) == 1 && p(1) == 0 && p.last() == 6)
      return true;
    if (p(0) == 2 && p(1) == 0 && p.last() == 6)
      return true;
    if (p(0) == 2 && p(1) == 0 && p.last() == 8)
      return true;
  }

  if (p.ndigits() == 5) {
    if (p(0) == 2 && p.last() == 2)
      return true;
    if (p(0) == 2 && p.last() == 4)
      return true;
    if (p(0) == 2 && p.last() == 6)
      return true;
    if (p(0) == 5 && p.last() == 2)
      return true;
    if (p(0) == 1 && p.last() == 6)
      return true;
    if (p(0) == 4 && p.last() == 2)
      return true;
  }
  return false;
}
/// PDG rule 15
/// The 9-digit penta-quark codes are±1nrnLnq1nq2nq3nq4nq5nJ, sorted such
/// thatnq1≥nq2≥nq3≥nq4. In the particle the first four are quarks and the fifth
/// an antiquark while t
/// heopposite holds in the antiparticle, which is given with a negative sign.
/// Thenr,nL, andnJnumbers have the same meaning as for ordinary hadrons.

template<>
inline bool isPentaquark(const DecodedPID& p)
{
  return (p.ndigits() == 9 && p(0) == 1 && p.max_digit(1, 6) <= 6
          && p.min_digit(1, 6) > 0
          && (p(3) >= p(4) && p(4) >= p(5) && p(5) >= p(6)));
}
/// PDG rule 14
/// The 9-digit tetra-quark codes are±1nrnLnq1nq20nq3nq4nJ. For the
/// particleq1q2is a diquarkand
/// ̄q3 ̄q4an antidiquark, sorted such thatnq1≥nq2,nq3≥nq4,nq1≥nq3,
/// andnq2≥nq4ifnq1=nq3.
/// For the antiparticle, given with a negative sign, ̄q1 ̄q2is an antidiquark
/// andq3q4a diquark,
/// with the same sorting except that eithernq1> nq3ornq2> nq4(so
/// thatflavour-diagonal states are particles). Thenr,nL, andnJnumbers have the
/// same meaningas for ordinary hadrons.
template<>
inline bool isTetraquark(const DecodedPID& p)
{
  return (p.ndigits() == 9 && p(0) == 1 && p(5) == 0 && p.max_digit(1, 3) <= 6
          && p.min_digit(1, 3) > 0 && p.max_digit(1 + 3, 3 + 3) <= 6
          && p.min_digit(1 + 3, 3 + 3) > 0 && (p(3) >= p(4) && p(6) >= p(7))
          && ((p(3) > p(6)) || (p(3) == p(6) && (p(4) >= p(7)))));
}

/// PDG rule 16:
/// Nuclear codes are given as 10-digit numbers±10LZZZAAAI. For a (hyper)nucleus
/// consistingofnpprotons,nnneutrons andnΛΛ’s,A=np+nn+nΛgives the total baryon
/// number, Z=np the total charge andL=nΛthe total number of strange
/// quarks.Igives the isomerlevel, with I= 0 corresponding to the ground state
/// and I >0 to excitations, see [2], wherestates denotedm,n,p,qtranslate to I=
/// 1–4. As examples, the deuteron is 1000010020 and 235U is 1000922350. To
/// avoid ambiguities, nuclear codes should not be applied to a singlehadron,
/// like p,n or Λ0, where quark-contents-based codes already exist.
template<>
inline bool isNucleus(const DecodedPID& p)
{
  if (std::abs(p.pid()) == PROTON)
    return true;
  return (p.ndigits() == 10 && p(0) == 1 && p(1) == 0);
}
/// APID: graviton and all Higgs extensions are BSM
template<>
inline bool isBSM(const DecodedPID& p)
{
  if (p.pid() == GRAVITON || std::abs(p.pid()) == MAVTOP
      || p.pid() == DARKPHOTON)
    return true;
  if (std::abs(p.pid()) > 16 && std::abs(p.pid()) < 19)
    return true;
  if (std::abs(p.pid()) > 31 && std::abs(p.pid()) < 38)
    return true;
  if (std::abs(p.pid()) > 39 && std::abs(p.pid()) < 81)
    return true;
  if (std::abs(p.pid()) > 6 && std::abs(p.pid()) < 9)
    return true;
  if (isSUSY(p))
    return true;
  if (isTechnicolor(p))
    return true;
  if (isExcited(p))
    return true;
  if (isKK(p))
    return true;
  if (isHiddenValley(p))
    return true;
  return false;
}

template<>
inline bool isSUSY(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSUSY(value_digits);
}
template<>
inline bool isDiquark(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isDiquark(value_digits);
}
template<>
inline bool isMeson(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isMeson(value_digits);
}
template<>
inline bool isBaryon(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isBaryon(value_digits);
}
template<>
inline bool isTetraquark(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isTetraquark(value_digits);
}
template<>
inline bool isPentaquark(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isPentaquark(value_digits);
}
template<>
inline bool isNucleus(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isNucleus(value_digits);
}
template<>
inline bool isBSM(const int& p)
{
  if (p == GRAVITON || std::abs(p) == MAVTOP || p == DARKPHOTON)
    return true;
  if (std::abs(p) > 16 && std::abs(p) < 19)
    return true;
  if (std::abs(p) > 31 && std::abs(p) < 38)
    return true;
  if (std::abs(p) > 39 && std::abs(p) < 81)
    return true;
  if (std::abs(p) > 6 && std::abs(p) < 9)
    return true;
  auto value_digits = DecodedPID(p);
  return isBSM(value_digits);
}

template<>
inline bool isHadron(const DecodedPID& p)
{
  return isMeson(p) || isBaryon(p) || isTetraquark(p) || isPentaquark(p);
}
template<>
inline bool isHadron(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isHadron(value_digits);
}
template<>
inline bool isTransportable(const DecodedPID& p)
{
  return isPhoton(p.pid()) || isGeantino(p.pid()) || isHadron(p)
      || isLepton(p.pid()) || p.pid() == DARKPHOTON;
}
template<>
inline bool isTransportable(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isTransportable(value_digits);
}
/// Av: we implement here an ATLAS-sepcific convention: all particles which are
/// 99xxxxx are fine.
template<>
inline bool isValid(const DecodedPID& p)
{
  return isHadron(p) || isTrajectory(p.pid()) || isDiquark(p) || isBSM(p)
      || isNucleus(p) || (std::abs(p.pid()) < 42) || isGenSpecific(p.pid())
      || isGeantino(p.pid()) || isPythia8Specific(p);
}
template<>
inline bool isValid(const int& p)
{
  if (!p)
    return false;
  if (std::abs(p) < 42)
    return true;
  if (isGenSpecific(p))
    return true;
  auto value_digits = DecodedPID(p);
  return isValid(value_digits);
}
template<class T>
inline bool hasQuark(const T& p, const int& q);

template<>
inline bool hasQuark(const DecodedPID& p, const int& q)
{
  if (isQuark(p.pid())) {
    return (std::abs(p.pid()) == q);
  }
  if (isMeson(p)) {
    return *(p.second.rbegin() + 1) == q || *(p.second.rbegin() + 2) == q;
  }
  if (isDiquark(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 3, q);
    return (i != p.second.rbegin() + 3);
  }
  if (isBaryon(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 4, q);
    return (i != p.second.rbegin() + 4);
  }
  if (isTetraquark(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 5, q);
    return (i != p.second.rbegin() + 5);
  }
  if (isPentaquark(p)) {
    auto i = std::find(p.second.rbegin() + 1, p.second.rbegin() + 6, q);
    return (i != p.second.rbegin() + 6);
  }
  if (isNucleus(p) && p.first != PROTON) {
    return q == 3 && p(2) > 0;
  }
  return false;
}

template<>
inline bool hasQuark(const int& p, const int& q)
{
  auto value_digits = DecodedPID(p);
  return hasQuark(value_digits, q);
}

template<class T>
inline int leadingQuark(const T& p)
{
  return leadingQuark(p->pdg_id());
}
template<>
inline int leadingQuark(const DecodedPID& p)
{
  if (isQuark(p.pid())) {
    return std::abs(p.pid());
  }
  if (isMeson(p)) {
    return p.max_digit(1, 3);
  }
  if (isDiquark(p)) {
    return p.max_digit(1, 3);
  }
  if (isBaryon(p)) {
    return p.max_digit(1, 4);
  }
  if (isTetraquark(p)) {
    return p.max_digit(1, 5);
  }
  if (isPentaquark(p)) {
    return p.max_digit(1, 6);
  }
  return 0;
}

template<>
inline int leadingQuark(const int& p)
{
  auto value_digits = DecodedPID(p);
  return leadingQuark(value_digits);
}

template<class T>
inline bool hasStrange(const T& p)
{
  return hasQuark(p, SQUARK);
}
template<class T>
inline bool hasCharm(const T& p)
{
  return hasQuark(p, CQUARK);
}
template<class T>
inline bool hasBottom(const T& p)
{
  return hasQuark(p, BQUARK);
}
template<class T>
inline bool hasTop(const T& p)
{
  return hasQuark(p, TQUARK);
}

template<class T>
inline bool isLightHadron(const T& p)
{
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isHadron(p);
}
template<class T>
inline bool isHeavyHadron(const T& p)
{
  auto lq = leadingQuark(p);
  return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isHadron(p);
}
template<class T>
inline bool isStrangeHadron(const T& p)
{
  return leadingQuark(p) == SQUARK && isHadron(p);
}
template<class T>
inline bool isCharmHadron(const T& p)
{
  return leadingQuark(p) == CQUARK && isHadron(p);
}
template<class T>
inline bool isBottomHadron(const T& p)
{
  return leadingQuark(p) == BQUARK && isHadron(p);
}
template<class T>
inline bool isTopHadron(const T& p)
{
  return leadingQuark(p) == TQUARK && isHadron(p);
}

template<class T>
inline bool isLightMeson(const T& p)
{
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isMeson(p);
}
template<class T>
inline bool isHeavyMeson(const T& p)
{
  auto lq = leadingQuark(p);
  return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isMeson(p);
}
template<class T>
inline bool isStrangeMeson(const T& p)
{
  return leadingQuark(p) == SQUARK && isMeson(p);
}
template<class T>
inline bool isCharmMeson(const T& p)
{
  return leadingQuark(p) == CQUARK && isMeson(p);
}
template<class T>
inline bool isBottomMeson(const T& p)
{
  return leadingQuark(p) == BQUARK && isMeson(p);
}
template<class T>
inline bool isTopMeson(const T& p)
{
  return leadingQuark(p) == TQUARK && isMeson(p);
}

template<class T>
inline bool isCCbarMeson(const T& p)
{
  return isCCbarMeson(p->pdg_id());
}
template<>
inline bool isCCbarMeson(const DecodedPID& p)
{
  return leadingQuark(p) == CQUARK && isMeson(p) && p(1) == CQUARK
      && p(2) == CQUARK;
}
template<>
inline bool isCCbarMeson(const int& p)
{
  return isCCbarMeson(DecodedPID(p));
}

template<class T>
inline bool isBBbarMeson(const T& p)
{
  return isBBbarMeson(p->pdg_id());
}
template<>
inline bool isBBbarMeson(const DecodedPID& p)
{
  return leadingQuark(p) == BQUARK && isMeson(p) && p(1) == BQUARK
      && p(2) == BQUARK;
}
template<>
inline bool isBBbarMeson(const int& p)
{
  return isCCbarMeson(DecodedPID(p));
}

template<class T>
inline bool isLightBaryon(const T& p)
{
  auto lq = leadingQuark(p);
  return (lq == DQUARK || lq == UQUARK || lq == SQUARK) && isBaryon(p);
}
template<class T>
inline bool isHeavyBaryon(const T& p)
{
  auto lq = leadingQuark(p);
  return (lq == CQUARK || lq == BQUARK || lq == TQUARK) && isBaryon(p);
}
template<class T>
inline bool isStrangeBaryon(const T& p)
{
  return leadingQuark(p) == SQUARK && isBaryon(p);
}
template<class T>
inline bool isCharmBaryon(const T& p)
{
  return leadingQuark(p) == CQUARK && isBaryon(p);
}
template<class T>
inline bool isBottomBaryon(const T& p)
{
  return leadingQuark(p) == BQUARK && isBaryon(p);
}
template<class T>
inline bool isTopBaryon(const T& p)
{
  return leadingQuark(p) == TQUARK && isBaryon(p);
}

template<class T>
inline int charge3(const T& p)
{
  return charge3(p->pdg_id());
}
template<class T>
inline double charge(const T& p)
{
  return 1.0 * charge3(p) / 3.0;
}
template<class T>
inline double threeCharge(const T& p)
{
  return charge3(p);
}
template<class T>
inline bool isCharged(const T& p)
{
  return charge3(p) != 0;
}
template<class T>
inline bool isNeutral(const T& p)
{
  return charge3(p) == 0;
}

template<>
inline int charge3(const DecodedPID& p)
{
  auto ap = std::abs(p.pid());
  if (ap < TABLESIZE)
    return p.pid() > 0 ? triple_charge.at(ap) : -triple_charge.at(ap);
  if (ap == K0)
    return 0;
  if (ap == GEANTINO0)
    return 0;
  if (ap == GEANTINOPLUS)
    return p.pid() > 0 ? 3 : -3;
  if (ap == MAVTOP)
    return p.pid() > 0 ? 2 : -2;
  size_t nq = 0;
  int sign = 1;
  int signmult = 1;
  int result = 0;
  bool classified = false;
  if (!classified && isMeson(p)) {
    classified = true;
    nq = 2;
    if ((*(p.second.rbegin() + 2)) == 2 || (*(p.second.rbegin() + 2)) == 4) {
      sign = -1;
    }
    signmult = -1;
  }
  if (!classified && isDiquark(p)) {
    return triple_charge.at(p(0)) + triple_charge.at(p(1));
  }
  if (!classified && isBaryon(p)) {
    classified = true;
    nq = 3;
  }
  if (!classified && isTetraquark(p)) {
    return triple_charge.at(p(3)) + triple_charge.at(p(4))
        - triple_charge.at(p(6)) - triple_charge.at(p(7));
  }
  if (!classified && isPentaquark(p)) {
    return triple_charge.at(p(3)) + triple_charge.at(p(4))
        + triple_charge.at(p(5)) + triple_charge.at(p(6))
        - triple_charge.at(p(7));
  }
  if (!classified && isNucleus(p)) {
    classified = true;
    nq = 0;
    result = 3 * (p(3) * 100 + p(4) * 10 + p(5)) + (-1) * p(2);
  }
  if (!classified && isSUSY(p)) {
    nq = 0;
    auto pp = p.shift(1);
    if (pp.ndigits() > 2)
      pp = pp.shift(1);
    return charge3(pp);
  }
  for (auto r = p.second.rbegin() + 1; r != p.second.rbegin() + 1 + nq; ++r) {
    result += triple_charge.at(*r) * sign;
    sign *= signmult;
  }
  return p.pid() > 0 ? result : -result;
}
template<>
inline int charge3(const int& p)
{
  int ap = std::abs(p);
  if (ap < TABLESIZE)
    return p > 0 ? triple_charge.at(ap) : -triple_charge.at(ap);
  auto value_digits = DecodedPID(p);
  return charge3(value_digits);
}

template<class T>
inline bool isEMInteracting(const T& p)
{
  return isEMInteracting(p->pdg_id());
}
template<>
inline bool isEMInteracting(const int& p)
{
  return (isPhoton(p) || isZ(p) || charge3(p) != 0);
}

template<class T>
inline bool isStrongInteracting(const T& p)
{
  return isStrongInteracting(p->pdg_id());
}
template<>
inline bool isStrongInteracting(const int& p)
{
  return (isGluon(p) || isQuark(p) || isDiquark(p) || isLeptoQuark(p)
          || isHadron(p));
}

template<class T>
inline bool isParton(const T& p)
{
  return isQuark(p) || isGluon(p);
}

namespace SUSY
{
static const int SUSYGLUONCODE = 9;
template<class T>
inline bool isRGlueball(const T& p)
{
  return isRGlueball(p->pdg_id());
}
template<>
inline bool isRGlueball(const DecodedPID& p)
{
  auto pp = p.shift(1).shift(1);
  size_t ng = 0;
  for (size_t i = 1; i + 1 < pp.ndigits(); ++i) {
    if (pp(i) == SUSYGLUONCODE)
      ng++;
  }
  return p(1) == 9 && p(2) == 9 && p.ndigits() == 7 && ng > 0;
}
template<>
inline bool isRGlueball(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isRGlueball(value_digits);
}
template<class T>
inline bool isRHadron(const T& p)
{
  return isRHadron(p->pdg_id());
}
template<>
inline bool isRHadron(const DecodedPID& p)
{
  auto pp = p.shift(1);
  if (pp(1) == 7 || pp(1) == 8)
    return false;
  if (pp.ndigits() > 2)
    pp = pp.shift(1);
  return isSUSY(p) && (isHadron(pp) || isRGlueball(p));
}
template<>
inline bool isRHadron(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isRHadron(value_digits);
}

template<class T>
inline bool isRMeson(const T& p)
{
  return isRMeson(p->pdg_id());
}
template<>
inline bool isRMeson(const DecodedPID& p)
{
  auto pp = p.shift(1);
  if (pp(1) == 7 || pp(1) == 8)
    return false;
  if (pp.ndigits() > 2)
    pp = pp.shift(1);
  return isSUSY(p) && isMeson(pp);
}
template<>
inline bool isRMeson(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isRMeson(value_digits);
}

template<class T>
inline bool isRBaryon(const T& p)
{
  return isRBaryon(p->pdg_id());
}
template<>
inline bool isRBaryon(const DecodedPID& p)
{
  auto pp = p.shift(1);
  if (pp(1) == 7 || pp(1) == 8)
    return false;
  if (pp.ndigits() > 2)
    pp = pp.shift(1);
  return isSUSY(p) && isBaryon(pp);
}
template<>
inline bool isRBaryon(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isRBaryon(value_digits);
}

template<class T>
inline bool isRBottomHadron(const T& p)
{
  return isRBottomHadron(p->pdg_id());
}
template<>
inline bool isRBottomHadron(const DecodedPID& p)
{
  auto pp = p.shift(1);
  if (pp(1) == 7 || pp(1) == 8)
    return false;
  if (pp.ndigits() > 2)
    pp = pp.shift(1);
  return isSUSY(p) && isBottomHadron(pp);
}
template<>
inline bool isRBottomHadron(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isRBottomHadron(value_digits);
}

template<class T>
inline bool isRTopHadron(const T& p)
{
  return isRTopHadron(p->pdg_id());
}
template<>
inline bool isRTopHadron(const DecodedPID& p)
{
  auto pp = p.shift(1);
  if (pp(1) == 7 || pp(1) == 8)
    return false;
  if (pp.ndigits() > 2)
    pp = pp.shift(1);
  return isSUSY(p) && isTopHadron(pp);
}
template<>
inline bool isRTopHadron(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isRTopHadron(value_digits);
}

template<class T>
inline bool isSLepton(const T& p)
{
  return isSLepton(p->pdg_id());
}
template<>
inline bool isSLepton(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isLepton(pp);
}
template<>
inline bool isSLepton(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSLepton(value_digits);
}

template<class T>
inline bool isSTopHadron(const T& p)
{
  return isSTopHadron(p->pdg_id());
}
template<>
inline bool isSTopHadron(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isTopHadron(pp);
}
template<>
inline bool isSTopHadron(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSTopHadron(value_digits);
}

template<class T>
inline bool isSBaryon(const T& p)
{
  return isSBaryon(p->pdg_id());
}
template<>
inline bool isSBaryon(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isBaryon(pp);
}
template<>
inline bool isSBaryon(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSBaryon(value_digits);
}

template<class T>
inline bool isSMeson(const T& p)
{
  return isSMeson(p->pdg_id());
}
template<>
inline bool isSMeson(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isMeson(pp);
}
template<>
inline bool isSMeson(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSMeson(value_digits);
}

template<class T>
inline bool isSBottomHadron(const T& p)
{
  return isSBottomHadron(p->pdg_id());
}
template<>
inline bool isSBottomHadron(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isBottomHadron(pp);
}
template<>
inline bool isSBottomHadron(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSBottomHadron(value_digits);
}

template<class T>
inline bool isSTopMeson(const T& p)
{
  return isSTopMeson(p->pdg_id());
}
template<>
inline bool isSTopMeson(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isTopMeson(pp);
}
template<>
inline bool isSTopMeson(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSTopMeson(value_digits);
}

template<class T>
inline bool isSBottomMeson(const T& p)
{
  return isSBottomMeson(p->pdg_id());
}
template<>
inline bool isSBottomMeson(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isBottomMeson(pp);
}
template<>
inline bool isSBottomMeson(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSBottomMeson(value_digits);
}

template<class T>
inline bool isSTopBaryon(const T& p)
{
  return isSTopBaryon(p->pdg_id());
}
template<>
inline bool isSTopBaryon(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isTopBaryon(pp);
}
template<>
inline bool isSTopBaryon(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSTopBaryon(value_digits);
}

template<class T>
inline bool isSBottomBaryon(const T& p)
{
  return isSBottomBaryon(p->pdg_id());
}
template<>
inline bool isSBottomBaryon(const DecodedPID& p)
{
  auto pp = p.shift(1);
  return isSUSY(p) && isBottomBaryon(pp);
}
template<>
inline bool isSBottomBaryon(const int& p)
{
  auto value_digits = DecodedPID(p);
  return isSBottomBaryon(value_digits);
}

template<class T>
inline bool spin(const T& p)
{
  return spin(p->pdg_id());
}
template<>
inline bool spin(const int& p)
{
  return p % 10;
}

inline std::vector<int> containedQuarks(int p)
{
  auto pp = DecodedPID(p);
  if (isSUSY(pp)) {
    pp = pp.shift(1);
    if (pp.ndigits() > 2)
      pp = pp.shift(1);
  }
  std::vector<int> quarks;
  for (int i = 1; i <= 6; ++i)
    if (hasQuark(pp, i))
      quarks.push_back(i);
  return quarks;
}

}  // namespace SUSY

}  // namespace ParticleData
