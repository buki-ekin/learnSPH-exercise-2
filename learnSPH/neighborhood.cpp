#include "neighborhood.h"

#include <stdexcept>

learnSPH::NeighborList learnSPH::brute_force_neighbors(
	const std::vector<Eigen::Vector3d>& positions,
	const double search_radius)
{
	if (!(search_radius > 0.0)) {
		throw std::invalid_argument("Search radius must be positive.");
	}

	const double radius_squared = search_radius * search_radius;
	NeighborList result(positions.size());
	for (std::size_t i = 0; i < positions.size(); ++i) {
		for (std::size_t j = 0; j < positions.size(); ++j) {
			if (i == j) {
				continue;
			}
			if ((positions[i] - positions[j]).squaredNorm() < radius_squared) {
				result[i].push_back(static_cast<unsigned int>(j));
			}
		}
	}
	return result;
}

double learnSPH::average_neighbor_count(const NeighborList& neighbors)
{
	if (neighbors.empty()) {
		return 0.0;
	}
	std::size_t total = 0;
	for (const auto& list : neighbors) {
		total += list.size();
	}
	return static_cast<double>(total) / static_cast<double>(neighbors.size());
}

std::vector<Eigen::Vector3d> learnSPH::regular_grid(
	const std::size_t nx,
	const std::size_t ny,
	const std::size_t nz,
	const double spacing)
{
	if (!(spacing > 0.0)) {
		throw std::invalid_argument("Grid spacing must be positive.");
	}
	std::vector<Eigen::Vector3d> result;
	result.reserve(nx * ny * nz);
	for (std::size_t i = 0; i < nx; ++i) {
		for (std::size_t j = 0; j < ny; ++j) {
			for (std::size_t k = 0; k < nz; ++k) {
				result.emplace_back(i * spacing, j * spacing, k * spacing);
			}
		}
	}
	return result;
}
