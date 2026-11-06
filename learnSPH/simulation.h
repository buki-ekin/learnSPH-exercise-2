#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

namespace learnSPH
{
	void semi_implicit_euler(
		std::vector<Eigen::Vector3d>& positions,
		std::vector<Eigen::Vector3d>& velocities,
		const std::vector<Eigen::Vector3d>& accelerations,
		double time_step);

	double cfl_time_step(
		const std::vector<Eigen::Vector3d>& velocities,
		double particle_radius,
		double cfl_factor,
		double minimum_time_step,
		double maximum_time_step);

	std::size_t remove_outside_aabb(
		std::vector<Eigen::Vector3d>& positions,
		std::vector<Eigen::Vector3d>& velocities,
		std::vector<double>& scalar_data,
		const Eigen::Vector3d& minimum,
		const Eigen::Vector3d& maximum);
}
