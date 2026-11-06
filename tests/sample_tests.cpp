#include <algorithm>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/kernel.h"
#include "../learnSPH/neighborhood.h"
#include "../learnSPH/simulation.h"

using Catch::Approx;

TEST_CASE("Cubic spline covers all three branches", "[kernel]")
{
	const double alpha = 3.0 / (2.0 * learnSPH::kernel::PI);
	CHECK(learnSPH::kernel::cubic_spline(0.0) == Approx(1.0 / learnSPH::kernel::PI));
	CHECK(learnSPH::kernel::cubic_spline(0.5) == Approx(alpha * (2.0/3.0 - 0.25 + 0.0625)));
	CHECK(learnSPH::kernel::cubic_spline(1.0) == Approx(alpha / 6.0));
	CHECK(learnSPH::kernel::cubic_spline(1.5) == Approx(alpha * 0.125 / 6.0));
	CHECK(learnSPH::kernel::cubic_spline(2.0) == 0.0);
	CHECK(learnSPH::kernel::cubic_spline(3.0) == 0.0);
}

TEST_CASE("Full cubic spline is symmetric, non-negative and compact", "[kernel]")
{
	const Eigen::Vector3d a(0.3, -0.2, 0.7);
	const Eigen::Vector3d b(-0.1, 0.4, 0.8);
	const double h = 0.6;
	CHECK(learnSPH::kernel::cubic_spline(a, b, h) ==
	      Approx(learnSPH::kernel::cubic_spline(b, a, h)));
	CHECK(learnSPH::kernel::cubic_spline(a, b, h) >= 0.0);
	CHECK(learnSPH::kernel::cubic_spline({0, 0, 0}, {2.0*h, 0, 0}, h) == 0.0);
	CHECK(learnSPH::kernel::cubic_spline({0, 0, 0}, {0, 0, 0}, h) ==
	      Approx(1.0 / (learnSPH::kernel::PI * h*h*h)));
}

TEST_CASE("Analytical gradient agrees with central differences", "[kernel]")
{
	const Eigen::Vector3d x_i(0.31, -0.17, 0.22);
	const Eigen::Vector3d x_j(-0.08, 0.12, 0.03);
	const double h = 0.7;
	const double epsilon = 1e-6;
	const Eigen::Vector3d analytical =
		learnSPH::kernel::cubic_spline_gradient(x_i, x_j, h);
	Eigen::Vector3d numerical;
	for (int axis = 0; axis < 3; ++axis) {
		Eigen::Vector3d offset = Eigen::Vector3d::Zero();
		offset[axis] = epsilon;
		numerical[axis] =
			(learnSPH::kernel::cubic_spline(x_i + offset, x_j, h) -
			 learnSPH::kernel::cubic_spline(x_i - offset, x_j, h)) /
			(2.0 * epsilon);
	}
	CHECK((analytical - numerical).norm() < 1e-8);
}

TEST_CASE("Kernel gradient at zero distance is analytically zero", "[kernel]")
{
	const Eigen::Vector3d x(0.2, 0.4, -0.3);
	CHECK(learnSPH::kernel::cubic_spline_gradient(x, x, 0.5).isZero());
}

TEST_CASE("Precomputed kernel matches straightforward implementation", "[kernel]")
{
	const double h = 0.37;
	const Eigen::Vector3d difference(0.11, -0.08, 0.21);
	learnSPH::kernel::CubicSplineKernel optimized(h);
	CHECK(optimized.value(difference) ==
	      Approx(learnSPH::kernel::cubic_spline(difference, Eigen::Vector3d::Zero(), h)));
	CHECK((optimized.gradient(difference) -
	       learnSPH::kernel::cubic_spline_gradient(
	           difference, Eigen::Vector3d::Zero(), h)).norm() < 1e-12);
}

TEST_CASE("Brute force and CompactNSearch return the same neighbors", "[neighbors]")
{
	const auto positions = learnSPH::regular_grid(6, 6, 6, 0.006);
	const double radius = 0.0144;
	const auto brute = learnSPH::brute_force_neighbors(positions, radius);
	CompactNSearch::NeighborhoodSearch search(radius);
	search.add_point_set(positions.front().data(), positions.size(), true, true, true);
	search.find_neighbors();
	const auto& compact = search.point_set(0);
	for (unsigned int i = 0; i < positions.size(); ++i) {
		std::vector<unsigned int> actual;
		for (unsigned int k = 0; k < compact.n_neighbors(0, i); ++k) {
			actual.push_back(compact.neighbor(0, i, k));
		}
		auto expected = brute[i];
		std::sort(actual.begin(), actual.end());
		std::sort(expected.begin(), expected.end());
		CHECK(actual == expected);
	}
}

TEST_CASE("Grid has about fifty neighbors per particle", "[neighbors]")
{
	const auto positions = learnSPH::regular_grid(20, 20, 20, 0.006);
	CompactNSearch::NeighborhoodSearch search(0.0144);
	search.add_point_set(positions.front().data(), positions.size(), true, true, true);
	search.find_neighbors();
	const auto& points = search.point_set(0);
	std::size_t total = 0;
	for (unsigned int i = 0; i < points.n_points(); ++i) {
		total += points.n_neighbors(0, i);
	}
	const double average = static_cast<double>(total) / points.n_points();
	CHECK(average > 45.0);
	CHECK(average < 56.0);
}

TEST_CASE("Semi-implicit Euler updates velocity before position", "[integration]")
{
	std::vector<Eigen::Vector3d> positions{{0.0, 0.0, 10.0}};
	std::vector<Eigen::Vector3d> velocities{{1.0, 0.0, 0.0}};
	const std::vector<Eigen::Vector3d> accelerations{{0.0, 0.0, -9.81}};
	learnSPH::semi_implicit_euler(positions, velocities, accelerations, 0.1);
	CHECK(velocities[0].x() == Approx(1.0));
	CHECK(velocities[0].z() == Approx(-0.981));
	CHECK(positions[0].x() == Approx(0.1));
	CHECK(positions[0].z() == Approx(9.9019));
}

TEST_CASE("CFL time step is clamped and handles rest", "[integration]")
{
	CHECK(learnSPH::cfl_time_step({Eigen::Vector3d::Zero()}, 0.05, 0.5, 0.001, 0.01)
	      == Approx(0.01));
	CHECK(learnSPH::cfl_time_step({{5.0, 0.0, 0.0}}, 0.05, 0.5, 0.001, 0.01)
	      == Approx(0.005));
	CHECK(learnSPH::cfl_time_step({{100.0, 0.0, 0.0}}, 0.05, 0.5, 0.001, 0.01)
	      == Approx(0.001));
}

TEST_CASE("AABB removal keeps particle arrays aligned", "[integration]")
{
	std::vector<Eigen::Vector3d> positions{{0,0,0}, {3,0,0}, {0,1,0}};
	std::vector<Eigen::Vector3d> velocities{{1,0,0}, {2,0,0}, {3,0,0}};
	std::vector<double> ids{10, 20, 30};
	const auto removed = learnSPH::remove_outside_aabb(
		positions, velocities, ids, {-1,-1,-1}, {1,1,1});
	CHECK(removed == 1);
	REQUIRE(positions.size() == 2);
	CHECK(ids[0] == 10);
	CHECK(ids[1] == 30);
	CHECK(velocities[1].x() == 3.0);
}
