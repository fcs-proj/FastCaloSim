#ifndef ParticleID_h
#define ParticleID_h

namespace HepPDT
{

class ParticleID
{
public:
  ParticleID(const int pdgID);
  inline double charge() const { return m_charge; }

private:
  double m_charge;
};

} // namespace HepPDT

#endif
