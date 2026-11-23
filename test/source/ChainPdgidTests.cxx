// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#include <memory>

#include <TClass.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Core/TFCSParametrization.h"
#include "FastCaloSim/Core/TFCSParametrizationBase.h"
#include "FastCaloSim/Core/TFCSParametrizationChain.h"

namespace
{
// An object created through ROOT, as when it is read from a parametrization
// file. These classes are not exported, so the test owns them through the
// exported base class.
template<typename T>
struct Handle
{
  std::unique_ptr<TFCSParametrizationBase> owner;

  explicit Handle(const char* class_name)
      : owner(static_cast<TFCSParametrizationBase*>(
            TClass::GetClass(class_name)->New()))
  {
  }

  auto operator->() const -> T* { return static_cast<T*>(owner.get()); }
};

// A parametrization valid for all particles, without particle IDs of its own
auto all_particles() -> std::unique_ptr<TFCSParametrizationBase>
{
  auto param = std::make_unique<TFCSParametrizationBase>("all", "all");
  param->set_match_all_pdgid();
  return param;
}

// A parametrization valid for one particle
auto particle(int pdgid) -> Handle<TFCSParametrization>
{
  Handle<TFCSParametrization> param("TFCSParametrization");
  param->set_pdgid(pdgid);
  return param;
}
}  // namespace

// A chain runs all its parametrizations: it is valid for the particles all of
// them are valid for. One valid for all particles does not restrict that.
TEST(ChainPdgidTest, IntersectionWithParametrizationForAllParticles)
{
  const Handle<TFCSParametrizationChain> chain("TFCSParametrizationChain");
  const auto all = all_particles();
  const auto photon = particle(22);

  chain->push_back(all.get());
  EXPECT_TRUE(chain->is_match_all_pdgid());
  EXPECT_TRUE(chain->is_match_pdgid(22));

  chain->push_back(photon.owner.get());
  EXPECT_FALSE(chain->is_match_all_pdgid());
  EXPECT_TRUE(chain->is_match_pdgid(22));
  EXPECT_FALSE(chain->is_match_pdgid(11));

  // Disjoint particles: valid for none
  const auto electron = particle(11);
  chain->push_back(electron.owner.get());
  EXPECT_FALSE(chain->is_match_all_pdgid());
  EXPECT_FALSE(chain->is_match_pdgid(22));
  EXPECT_FALSE(chain->is_match_pdgid(11));
}

// A PDG ID select chain runs the parametrization of the particle: it is valid
// for the particles any of them is valid for, so for all particles if one is.
TEST(ChainPdgidTest, UnionWithParametrizationForAllParticles)
{
  const Handle<TFCSParametrizationChain> chain(
      "TFCSParametrizationPDGIDSelectChain");
  const auto photon = particle(22);
  const auto electron = particle(11);
  const auto all = all_particles();

  chain->push_back(photon.owner.get());
  chain->push_back(electron.owner.get());
  EXPECT_FALSE(chain->is_match_all_pdgid());
  EXPECT_TRUE(chain->is_match_pdgid(22));
  EXPECT_TRUE(chain->is_match_pdgid(11));
  EXPECT_FALSE(chain->is_match_pdgid(211));

  chain->push_back(all.get());
  EXPECT_TRUE(chain->is_match_all_pdgid());
  EXPECT_TRUE(chain->is_match_pdgid(211));
}
