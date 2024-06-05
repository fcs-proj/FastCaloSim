/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// BFieldSolenoid.cxx
//
#include <algorithm>
#include <cmath>
#include <string>

#include "BFieldSolenoid.h"
// #include "TTree.h"
using namespace std;

//
// read an ASCII field map from istream
// return 0 if successful
//
int BFieldSolenoid::readMap(istream& input)
{
  const double meter(1000.0);  // meter in mm
  const double gauss(1.0e-7);  // gauss in kT
  const string myname("BFieldSolenoid::readMap()");
  if (m_orig == m_tilt)
    delete m_orig;
  else {
    delete m_orig;
    delete m_tilt;
  }
  m_orig = m_tilt = new BFieldMesh<double>;
  // first line contains version
  int version;
  input >> version;
  if (version != 4) {
    cerr << myname << ": version number is " << version << " instead of 4"
         << endl;
    return 1;
  }
  // second line contains number of bins
  int nz, nr, nphi;
  double unused;
  input >> nz >> nr >> nphi >> unused >> unused >> unused;
  m_orig->reserve(nz, nr, nphi + 1);  // extra phi at 2pi
  // read mesh definitions
  double x;
  for (int i = 0; i < nr; i++) {
    input >> x;
    m_orig->appendMesh(1, x * meter);
  }
  for (int i = 0; i < nz; i++) {
    input >> x;
    m_orig->appendMesh(0, x * meter);
  }
  for (int i = 0; i < nphi; i++) {
    input >> x;
    m_orig->appendMesh(2, x);
  }
  m_orig->appendMesh(2, 2 * M_PI);  // extra phi at 2pi
  double zmin = m_orig->mesh(0, 0);
  double zmax = m_orig->mesh(0, nz - 1);
  double rmin = m_orig->mesh(1, 0);
  double rmax = m_orig->mesh(1, nr - 1);
  double phimin = m_orig->mesh(2, 0);
  double phimax = m_orig->mesh(2, nphi);
  m_orig->setRange(zmin, zmax, rmin, rmax, phimin, phimax);
  // skip one record
  input >> unused >> unused >> unused >> unused >> unused;
  // read field values
  // the index order is opposite to the toroid - have to reorder
  // also convert from gauss to Tesla here.
  vector<double> Bz, Br, Bphi;
  double b;
  for (int k = 0; k < nphi; k++) {
    for (int j = 0; j < nr; j++) {
      for (int i = 0; i < nz; i++) {
        input >> b;
        Bz.push_back(b * gauss);
      }
    }
    for (int j = 0; j < nr; j++) {
      for (int i = 0; i < nz; i++) {
        input >> b;
        Br.push_back(b * gauss);
      }
    }
    for (int j = 0; j < nr; j++) {
      for (int i = 0; i < nz; i++) {
        input >> b;
        Bphi.push_back(b * gauss);
      }
    }
  }
  for (int i = 0; i < nz; i++) {
    for (int j = 0; j < nr; j++) {
      for (int k = 0; k < nphi; k++) {
        int index = i + nz * (j + nr * k);
        BFieldVector<double> field(Bz[index], Br[index], Bphi[index]);
        m_orig->appendField(field);
      }
      // close phi at 2pi
      int index = i + nz * j;
      BFieldVector<double> field(Bz[index], Br[index], Bphi[index]);
      m_orig->appendField(field);
    }
  }
  // build the LUT
  m_orig->buildLUT();

  return 0;
}

//
// wrire the map to a ROOT file
// if tilted = true, write the moved-and-tilted map.
// ohterwise, write the original map.
//
// void
// BFieldSolenoid::writeMap( TFile* rootfile, bool tilted )
//{
//    BFieldMesh<double> *map = tilted ? m_tilt : m_orig;
//    if ( map == 0 ) return; // no map to write
//    if ( rootfile == 0 ) return; // no file
//    if ( rootfile->cd() == false ) return; // could not make it current
//    directory
//    // define the tree
//    TTree* tree = new TTree( "BFieldSolenoid", "BFieldSolenoid version 4" );
//    double zmin = map->zmin();
//    double zmax = map->zmax();
//    double rmin = map->rmin();
//    double rmax = map->rmax();
//    double phimin = map->phimin();
//    double phimax = map->phimax();
//    int nmeshz = map->nmesh(0);
//    int nmeshr = map->nmesh(1);
//    int nmeshphi = map->nmesh(2);
//    double *meshz, *meshr, *meshphi;
//    int nfield = nmeshz*nmeshr*nmeshphi;
//    double *fieldz, *fieldr, *fieldphi;
//    meshz = new double[nmeshz];
//    meshr = new double[nmeshr];
//    meshphi = new double[nmeshphi];
//    fieldz = new double[nfield];
//    fieldr = new double[nfield];
//    fieldphi = new double[nfield];
//    // define the tree branches
//    tree->Branch( "zmin", &zmin, "zmin/D" );
//    tree->Branch( "zmax", &zmax, "zmax/D" );
//    tree->Branch( "rmin", &rmin, "rmin/D" );
//    tree->Branch( "rmax", &rmax, "rmax/D" );
//    tree->Branch( "phimin", &phimin, "phimin/D" );
//    tree->Branch( "phimax", &phimax, "phimax/D" );
//    tree->Branch( "nmeshz", &nmeshz, "nmeshz/I" );
//    tree->Branch( "meshz", meshz, "meshz[nmeshz]/D" );
//    tree->Branch( "nmeshr", &nmeshr, "nmeshr/I" );
//    tree->Branch( "meshr", meshr, "meshr[nmeshr]/D" );
//    tree->Branch( "nmeshphi", &nmeshphi, "nmeshphi/I" );
//    tree->Branch( "meshphi", meshphi, "meshphi[nmeshphi]/D" );
//    tree->Branch( "nfield", &nfield, "nfield/I" );
//    tree->Branch( "fieldz", fieldz, "fieldz[nfield]/D" );
//    tree->Branch( "fieldr", fieldr, "fieldr[nfield]/D" );
//    tree->Branch( "fieldphi", fieldphi, "fieldphi[nfield]/D" );
//    // fill the mesh and field arrays
//    for ( int j = 0; j < nmeshz; j++ ) {
//        meshz[j] = map->mesh(0,j);
//    }
//    for ( int j = 0; j < nmeshr; j++ ) {
//        meshr[j] = map->mesh(1,j);
//    }
//    for ( int j = 0; j < nmeshphi; j++ ) {
//        meshphi[j] = map->mesh(2,j);
//    }
//    for ( int j = 0; j < nfield; j++ ) {
//        const BFieldVector<double> f = map->field(j);
//        fieldz[j] = f.z();
//        fieldr[j] = f.r();
//        fieldphi[j] = f.phi();
//    }
//    // write
//    tree->Fill();
//    rootfile->Write();
//    // clean up
//    delete[] meshz;
//    delete[] meshr;
//    delete[] meshphi;
//    delete[] fieldz;
//    delete[] fieldr;
//    delete[] fieldphi;
//}

//
// read the map from a ROOT file.
// returns 0 if successful.
//
// int
// BFieldSolenoid::readMap( TFile* rootfile )
//{
//    if ( rootfile == 0 ) return 1; // no file
//    if ( rootfile->cd() == false ) return 2; // could not make it current
//    directory if ( m_orig == m_tilt ) delete m_orig; else { delete m_orig;
//    delete m_tilt; } m_orig = m_tilt = new BFieldMesh<double>;
//    // open the tree
//    TTree* tree = (TTree*)rootfile->Get("BFieldSolenoid");
//    if ( tree == 0 ) return 3; // no tree
//    double zmin, zmax, rmin, rmax, phimin, phimax;
//    int nmeshz, nmeshr, nmeshphi;
//    double *meshz, *meshr, *meshphi;
//    int nfield;
//    double *fieldz, *fieldr, *fieldphi;
//    //unsigned char *fbyte;
//    // define the fixed-sized branches first
//    tree->SetBranchAddress( "zmin", &zmin );
//    tree->SetBranchAddress( "zmax", &zmax );
//    tree->SetBranchAddress( "rmin", &rmin );
//    tree->SetBranchAddress( "rmax", &rmax );
//    tree->SetBranchAddress( "phimin", &phimin );
//    tree->SetBranchAddress( "phimax", &phimax );
//    tree->SetBranchAddress( "nmeshz", &nmeshz );
//    tree->SetBranchAddress( "nmeshr", &nmeshr );
//    tree->SetBranchAddress( "nmeshphi", &nmeshphi );
//    tree->SetBranchAddress( "nfield", &nfield );
//    // prepare arrays - need to know the maximum sizes
//    tree->GetEntry(0);
//    meshz = new double[nmeshz];
//    meshr = new double[nmeshr];
//    meshphi = new double[nmeshphi];
//    fieldz = new double[nfield];
//    fieldr = new double[nfield];
//    fieldphi = new double[nfield];
//    // define the variable length branches
//    tree->SetBranchAddress( "meshz", meshz );
//    tree->SetBranchAddress( "meshr", meshr );
//    tree->SetBranchAddress( "meshphi", meshphi );
//    tree->SetBranchAddress( "fieldz", fieldz );
//    tree->SetBranchAddress( "fieldr", fieldr );
//    tree->SetBranchAddress( "fieldphi", fieldphi );
//    // read again, and copy data
//    tree->GetEntry(0);
//    m_orig->setRange( zmin, zmax, rmin, rmax, phimin, phimax );
//    m_orig->reserve( nmeshz, nmeshr, nmeshphi );
//    for ( int j = 0; j < nmeshz; j++ ) {
//        m_orig->appendMesh( 0, meshz[j] );
//    }
//    for ( int j = 0; j < nmeshr; j++ ) {
//        m_orig->appendMesh( 1, meshr[j] );
//    }
//    for ( int j = 0; j < nmeshphi; j++ ) {
//        m_orig->appendMesh( 2, meshphi[j] );
//    }
//    for ( int j = 0; j < nfield; j++ ) {
//        BFieldVector<double> field( fieldz[j], fieldr[j], fieldphi[j] );
//        m_orig->appendField( field );
//    }
//    // clean up
//    tree->Delete();
//    delete[] meshz;
//    delete[] meshr;
//    delete[] meshphi;
//    delete[] fieldz;
//    delete[] fieldr;
//    delete[] fieldphi;
//    // build the LUTs
//    m_orig->buildLUT();
//
//    return 0;
//}

//
// Returns the magnetic field at any position.
//
void BFieldSolenoid::getB(const double* xyz, double* B, double* deriv) const
{
  double z = xyz[2];
  double r = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
  double phi = atan2(xyz[1], xyz[0]);
  if (m_tilt && m_tilt->inside(z, r, phi)) {
    m_tilt->getB(xyz, B, deriv);
  } else {
    B[0] = B[1] = B[2] = 0.0;
    if (deriv)
      for (int i = 0; i < 9; i++)
        deriv[i] = 0.0;
  }
}

//
// Move and tilt the solenoid.
// Modify the m_tilt copy and keep the m_orig copy.
//
void BFieldSolenoid::moveMap(
    double dx, double dy, double dz, double ax, double ay)
{
  if (m_orig == 0) {
    cerr << "BFieldSolenoid::moveMap() : original map has not been read"
         << endl;
    return;
  }
  if (m_tilt != m_orig)
    delete m_tilt;
  //
  // copy the mesh and make it a bit smaller
  //
  const double zlim = 2820.;  // mm
  const double rlim = 1075.;  // mm
  m_tilt = new BFieldMesh<double>(-zlim, zlim, 0.0, rlim, 0.0, 2 * M_PI, 1.0);
  // z
  m_tilt->appendMesh(0, -zlim);
  for (unsigned i = 0; i < m_orig->nmesh(0); i++) {
    if (abs(m_orig->mesh(0, i)) < zlim)
      m_tilt->appendMesh(0, m_orig->mesh(0, i));
  }
  m_tilt->appendMesh(0, zlim);
  // r
  for (unsigned i = 0; i < m_orig->nmesh(1); i++) {
    if (m_orig->mesh(1, i) < rlim)
      m_tilt->appendMesh(1, m_orig->mesh(1, i));
  }
  m_tilt->appendMesh(1, rlim);
  // phi (no change)
  for (unsigned i = 0; i < m_orig->nmesh(2); i++) {
    m_tilt->appendMesh(2, m_orig->mesh(2, i));
  }
  //
  // loop over the new mesh, and compute the field at the
  // corresponding location in the original map
  //
  double sinax(sin(ax));
  double cosax(cos(ax));
  double sinay(sin(ay));
  double cosay(cos(ay));
  for (unsigned i = 0; i < m_tilt->nmesh(0); i++) {
    double z0(m_tilt->mesh(0, i));
    for (unsigned j = 0; j < m_tilt->nmesh(1); j++) {
      double r0(m_tilt->mesh(1, j));
      for (unsigned k = 0; k < m_tilt->nmesh(2); k++) {
        double phi0(m_tilt->mesh(2, k));
        double x0(r0 * cos(phi0));
        double y0(r0 * sin(phi0));
        // shift
        double x1(x0 - dx);
        double y1(y0 - dy);
        double z1(z0 - dz);
        // rotate around x by -ax
        double x2(x1);
        double y2(y1 * cosax + z1 * sinax);
        double z2(z1 * cosax - y1 * sinax);
        // rotate around y by -ay
        double xyz3[3] = {x2 * cosay - z2 * sinay, y2, z2 * cosay + x2 * sinay};
        // get (Bx,By,Bz) in the original frame
        double B[3];
        m_orig->getB(xyz3, B);
        // rotate around y by +ay
        double Bx1(B[0] * cosay + B[2] * sinay);
        double By1(B[1]);
        double Bz1(B[2] * cosay - B[0] * sinay);
        // rotate around x by +ax
        double Bx2(Bx1);
        double By2(By1 * cosax - Bz1 * sinax);
        double Bz2(Bz1 * cosax + By1 * sinax);
        // convert to cylindrical
        double cosphi0(cos(phi0));
        double sinphi0(sin(phi0));
        double Br(Bx2 * cosphi0 + By2 * sinphi0);
        double Bphi(-Bx2 * sinphi0 + By2 * cosphi0);
        BFieldVector<double> field(Bz2, Br, Bphi);
        m_tilt->appendField(field);
      }
    }
  }
  m_tilt->buildLUT();
}
