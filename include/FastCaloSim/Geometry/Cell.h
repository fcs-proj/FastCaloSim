#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

static auto normalize_angle(double angle) -> double
{
  angle = std::fmod(angle + M_PI, 2.0 * M_PI);
  if (angle < 0) {
    angle += 2.0 * M_PI;
  }
  return angle - M_PI;
}

std::pair<double, double> calc_phi_range(double phi, double dphi)
{
  double phi_min = normalize_angle(phi - 0.5 * dphi);
  double phi_max = normalize_angle(phi + 0.5 * dphi);

  // Ensure phi_min is less than phi_max correctly considering boundary crossing
  if (phi_min > phi_max) {
    std::swap(phi_min, phi_max);
  }
  return {phi_min, phi_max};
}

struct Cell
{
  long long id;
  double x;
  double y;
  double z;
  double phi;
  double eta;
  double r;
  double dx;
  double dy;
  double dz;
  double dphi;
  double deta;
  double dr;
  bool isBarrel;
  bool isCartesian;
  bool isCylindrical;
  bool isECCylindrical;
  int layer;

  void print() const
  {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Cell details:" << '\n';
    std::cout << "---------------------------" << '\n';
    std::cout << "ID: " << id << '\n';
    std::cout << "Layer: " << layer << '\n';
    std::cout << "Coordinates:" << '\n';
    std::cout << "  x: " << std::setw(10) << x << '\n';
    std::cout << "  y: " << std::setw(10) << y << '\n';
    std::cout << "  z: " << std::setw(10) << z << '\n';
    std::cout << "  phi: " << std::setw(10) << phi << '\n';
    std::cout << "  eta: " << std::setw(10) << eta << '\n';
    std::cout << "  r: " << std::setw(10) << r << '\n';
    std::cout << "Cell size:" << '\n';

    if (isCartesian) {
      std::cout << "  dx: " << std::setw(10) << dx << '\n';
      std::cout << "  dy: " << std::setw(10) << dy << '\n';
      std::cout << "  dz: " << std::setw(10) << dz << '\n';

      std::cout << "Ranges covered:\n";
      std::cout << "  x:   min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << x - 0.5 * dx << "  max: " << std::setw(10)
                << x + 0.5 * dx << '\n';
      std::cout << "  y:   min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << y - 0.5 * dy << "  max: " << std::setw(10)
                << y + 0.5 * dy << '\n';
      std::cout << "  z:   min: " << std::fixed << std::setprecision(2);
    } else if (isCylindrical) {
      std::cout << "  dphi: " << std::setw(10) << dphi << '\n';
      std::cout << "  deta: " << std::setw(10) << deta << '\n';
      std::cout << "  dr: " << std::setw(10) << dr << '\n';

      std::cout << "Ranges covered:\n";
      auto phi_range = calc_phi_range(phi, dphi);
      std::cout << "  phi: min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << phi_range.first
                << "  max: " << std::setw(10) << phi_range.second << '\n';
      std::cout << "  eta: min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << eta - 0.5 * deta
                << "  max: " << std::setw(10) << eta + 0.5 * deta << '\n';
      std::cout << "  r:   min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << r - 0.5 * dr << "  max: " << std::setw(10)
                << r + 0.5 * dr << '\n';

    } else if (isECCylindrical) {
      std::cout << "  dphi: " << std::setw(10) << dphi << '\n';
      std::cout << "  deta: " << std::setw(10) << deta << '\n';
      std::cout << "  dz: " << std::setw(10) << dz << '\n';

      std::cout << "Ranges covered:\n";
      auto phi_range = calc_phi_range(phi, dphi);
      std::cout << "  phi: min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << phi_range.first
                << "  max: " << std::setw(10) << phi_range.second << '\n';
      std::cout << "  eta: min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << eta - 0.5 * deta
                << "  max: " << std::setw(10) << eta + 0.5 * deta << '\n';
      std::cout << "  z:   min: " << std::fixed << std::setprecision(2)
                << std::setw(10) << z - 0.5 * dz << "  max: " << std::setw(10)
                << z + 0.5 * dz << '\n';
    }

    std::cout << "Flags:" << '\n';
    std::cout << "  isBarrel: " << isBarrel << '\n';
    std::cout << "  isCartesian: " << isCartesian << '\n';
    std::cout << "  isCylindrical: " << isCylindrical << '\n';
    std::cout << "  isECCylindrical: " << isECCylindrical << '\n';
    std::cout << "---------------------------" << '\n' << '\n';
  }
};
