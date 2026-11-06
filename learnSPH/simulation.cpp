#include "simulation.h"

#include <algorithm>
#include <stdexcept>

void learnSPH::semi_implicit_euler(
	std::vector<Eigen::Vector3d>& positions,
	std::vector<Eigen::Vector3d>& velocities,
	const std::vector<Eigen::Vector3d>& accelerations,
	const double time_step)
{
	if (positions.size() != velocities.size() || positions.size() != accelerations.size()) {
		throw std::invalid_argument("Position, velocity and acceleration arrays must have equal size.");
	}
	if (!(time_step > 0.0)) {
		throw std::invalid_argument("Time step must be positive.");
	}
	for (std::size_t i = 0; i < positions.size(); ++i) {
		velocities[i] += time_step * accelerations[i];
		positions[i] += time_step * velocities[i];
	}
}

double learnSPH::cfl_time_step(
	const std::vector<Eigen::Vector3d>& velocities,
	const double particle_radius,
	const double cfl_factor,
	const double minimum_time_step,
	const double maximum_time_step)
{
	if (!(particle_radius > 0.0) || !(cfl_factor > 0.0) ||
		!(minimum_time_step > 0.0) || maximum_time_step < minimum_time_step) {
		throw std::invalid_argument("Invalid CFL parameters.");
	}
	double maximum_speed = 0.0;
	for (const auto& velocity : velocities) {
		maximum_speed = std::max(maximum_speed, velocity.norm());
	}
	if (maximum_speed == 0.0) {
		return maximum_time_step;
	}
	return std::clamp(cfl_factor * particle_radius / maximum_speed,
	                  minimum_time_step,
	                  maximum_time_step);
}

std::size_t learnSPH::remove_outside_aabb(
	std::vector<Eigen::Vector3d>& positions,
	std::vector<Eigen::Vector3d>& velocities,
	std::vector<double>& scalar_data,
	const Eigen::Vector3d& minimum,
	const Eigen::Vector3d& maximum)
{
	if (positions.size() != velocities.size() || positions.size() != scalar_data.size()) {
		throw std::invalid_argument("Particle data arrays must have equal size.");
	}

	std::size_t write = 0;
	for (std::size_t read = 0; read < positions.size(); ++read) {
		const Eigen::Vector3d& p = positions[read];
		const bool inside = (p.array() >= minimum.array()).all() &&
		                    (p.array() <= maximum.array()).all();
		if (!inside) {
			continue;
		}
		if (write != read) {
			positions[write] = positions[read];
			velocities[write] = velocities[read];
			scalar_data[write] = scalar_data[read];
		}
		++write;
	}
	const std::size_t removed = positions.size() - write;
	positions.resize(write);
	velocities.resize(write);
	scalar_data.resize(write);
	return removed;
}
