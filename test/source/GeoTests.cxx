#include <chrono>

#include "GeoTests.h"

#include <CLHEP/Random/RanluxEngine.h>
#include <TRandom.h>
#include <gtest/gtest.h>

#include "FastCaloSim/Geometry/CaloCellRDTree.h"
#include "FastCaloSim/Geometry/Cell.h"

// Sample uniformly between two numbers
float sample(float a, float b)
{
  return a + (b - a) * gRandom->Rndm();
}

TEST_F(GeoTests, TestNewGeometrySetup)
{
  CaloCellRDTree tree(geo);
  int n_random_points = 100;

  auto display = geo.Display();
  std::cout << display->AsString() << '\n';

  int layer = 1;

  int nCells = tree.get_num_cells(layer);
  std::cout << "Number of cells " << nCells << std::endl;
  // Loop over all cells
  for (int i = 0; i < nCells; i++) {
    Cell cell = tree.get_cell_at_index(layer, i);

    // Loop over random points
    for (int irnd = 0; irnd < n_random_points; ++irnd) {
      float test_eta, test_phi, test_r, test_x, test_y, test_z;

      // cell described with deta, dphi, dr
      if (cell.isCylindrical) {
        // Sample random points that should be inside the cell
        test_eta =
            sample(cell.eta - 0.49 * cell.deta, cell.eta + 0.49 * cell.deta);
        test_phi = normalize_angle(
            sample(cell.phi - 0.49 * cell.dphi, cell.phi + 0.49 * cell.dphi));
        test_r = sample(cell.r - 0.49 * cell.dr, cell.r + 0.49 * cell.dr);

        test_x = test_r * std::cos(test_phi);
        test_y = test_r * std::sin(test_phi);
        test_z = test_r * std::sinh(test_eta);
      }
      if (cell.isECCylindrical) {
        // Sample random points that should be inside the cell
        test_eta =
            sample(cell.eta - 0.49 * cell.deta, cell.eta + 0.49 * cell.deta);
        test_phi = normalize_angle(
            sample(cell.phi - 0.49 * cell.dphi, cell.phi + 0.49 * cell.dphi));
        test_z = sample(cell.z - 0.49 * cell.dz, cell.z + 0.49 * cell.dz);

        test_r = test_z / std::sinh(test_eta);
        test_x = test_r * std::cos(test_phi);
        test_y = test_r * std::sin(test_phi);
      }
      if (cell.isCartesian) {
        // Sample random points that should be inside the cell
        test_x = sample(cell.x - 0.49 * cell.dx, cell.x + 0.49 * cell.dx);
        test_y = sample(cell.y - 0.49 * cell.dy, cell.y + 0.49 * cell.dy);
        test_z = sample(cell.z - 0.49 * cell.dz, cell.z + 0.49 * cell.dz);

        test_r = std::sqrt(test_x * test_x + test_y * test_y);
        test_phi = normalize_angle(std::atan2(test_y, test_x));
        // theta
        double theta = std::atan2(test_r, test_z);
        test_eta = -std::log(std::tan(theta / 2.0));
      }

      Cell best_cell = tree.get_cell(layer, test_x, test_y, test_z);

      bool is_inside = false;

      if (cell.isCylindrical) {
        is_inside = std::abs(test_eta - best_cell.eta) < 0.5 * best_cell.deta
            && std::abs(normalize_angle(test_phi - best_cell.phi))
                < 0.5 * best_cell.dphi
            && std::abs(test_r - best_cell.r) < 0.5 * best_cell.dr;
      } else if (cell.isECCylindrical) {
        is_inside = std::abs(test_eta - best_cell.eta) < 0.5 * best_cell.deta
            && std::abs(
                   normalize_angle(test_phi - best_cell.phi))  // revisit logic
                < 0.5 * best_cell.dphi
            && std::abs(test_z - best_cell.z) < 0.5 * best_cell.dz;

      } else if (cell.isCartesian) {
        is_inside = std::abs(test_x - best_cell.x) < 0.5 * best_cell.dx
            && std::abs(test_y - best_cell.y) < 0.5 * best_cell.dy
            && std::abs(test_z - best_cell.z) < 0.5 * best_cell.dz;
      }

      if (!is_inside) {
        std::cout << "Test point not inside is test_x=" << test_x
                  << " test_y=" << test_y << " test_z=" << test_z
                  << " test_eta=" << test_eta << " test_phi=" << test_phi
                  << " test_r=" << test_r << '\n';

        std::cout << "Checked cell is " << '\n';
        cell.print();
        std::cout << "Best cell is " << '\n';
        best_cell.print();
        return;
      }
    }
  }
}
