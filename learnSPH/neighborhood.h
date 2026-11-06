#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

namespace learnSPH
{
	using NeighborList = std::vector<std::vector<unsigned int>>;

	NeighborList brute_force_neighbors(
		const std::vector<Eigen::Vector3d>& positions,
		double search_radius);

	double average_neighbor_count(const NeighborList& neighbors);
	std::vector<Eigen::Vector3d> regular_grid(
		std::size_t nx,
		std::size_t ny,
		std::size_t nz,
		double spacing);
}
