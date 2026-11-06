#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../learnSPH/simulation.h"

namespace
{
struct FallResult
{
	bool touched;
	bool tunneled;
	double final_z;
};

FallResult fixed_step_ground_test(const double radius, const double time_step)
{
	std::vector<Eigen::Vector3d> positions{{0.0, 0.0, 10.0}};
	std::vector<Eigen::Vector3d> velocities{{0.0, 0.0, 0.0}};
	const std::vector<Eigen::Vector3d> accelerations{{0.0, 0.0, -9.81}};
	bool touched = false;
	for (double time = 0.0; time < 3.0 && !touched; time += time_step) {
		learnSPH::semi_implicit_euler(positions, velocities, accelerations, time_step);
		if (positions[0].z() >= -radius && positions[0].z() <= radius) {
			touched = true;
			velocities[0].setZero();
		}
	}
	return {touched, positions[0].z() < -radius, positions[0].z()};
}

FallResult cfl_ground_test(const double radius, const double maximum_time_step)
{
	std::vector<Eigen::Vector3d> positions{{0.0, 0.0, 10.0}};
	std::vector<Eigen::Vector3d> velocities{{0.0, 0.0, 0.0}};
	const std::vector<Eigen::Vector3d> accelerations{{0.0, 0.0, -9.81}};
	bool touched = false;
	double time = 0.0;
	while (time < 3.0 && !touched) {
		const double dt = learnSPH::cfl_time_step(
			velocities, radius, 0.5, 0.0001, maximum_time_step);
		learnSPH::semi_implicit_euler(positions, velocities, accelerations, dt);
		time += dt;
		if (positions[0].z() >= -radius && positions[0].z() <= radius) {
			touched = true;
			velocities[0].setZero();
		}
	}
	return {touched, positions[0].z() < -radius, positions[0].z()};
}
}

int main()
{
	std::cout << std::fixed << std::setprecision(8);
	std::cout << "EULER_ERROR\n";
	std::cout << "dt,numerical_z,analytical_z,absolute_error\n";
	for (const double dt : {0.1, 0.05, 0.01, 0.005}) {
		std::vector<Eigen::Vector3d> positions{{0.0, 0.0, 10.0}};
		std::vector<Eigen::Vector3d> velocities{{1.0, -0.5, 2.0}};
		const std::vector<Eigen::Vector3d> accelerations{{0.0, 0.0, -9.81}};
		const int steps = static_cast<int>(std::llround(1.0 / dt));
		for (int step = 0; step < steps; ++step) {
			learnSPH::semi_implicit_euler(positions, velocities, accelerations, dt);
		}
		const double analytical_z = 10.0 + 2.0 - 0.5 * 9.81;
		std::cout << dt << ',' << positions[0].z() << ',' << analytical_z << ','
		          << std::abs(positions[0].z() - analytical_z) << '\n';
	}

	std::cout << "GROUND_FIXED_STEP\n";
	std::cout << "radius,dt,touched,tunneled,final_z\n";
	for (const double radius : {0.05, 0.02}) {
		for (const double dt : {0.01, 0.05, 0.1}) {
			const auto result = fixed_step_ground_test(radius, dt);
			std::cout << radius << ',' << dt << ',' << result.touched << ','
			          << result.tunneled << ',' << result.final_z << '\n';
		}
	}

	std::cout << "GROUND_CFL\n";
	std::cout << "radius,max_dt,touched,tunneled,final_z\n";
	for (const double radius : {0.05, 0.02}) {
		const auto result = cfl_ground_test(radius, 0.01);
		std::cout << radius << ",0.01," << result.touched << ','
		          << result.tunneled << ',' << result.final_z << '\n';
	}
	return 0;
}
