#include <chrono>
#include <iomanip>
#include <iostream>

#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/neighborhood.h"

namespace
{
using Clock = std::chrono::steady_clock;

template<class Function>
double average_milliseconds(Function&& function, const int measured_runs)
{
	function();
	double total = 0.0;
	for (int run = 0; run < measured_runs; ++run) {
		const auto start = Clock::now();
		function();
		const auto end = Clock::now();
		total += std::chrono::duration<double, std::milli>(end - start).count();
	}
	return total / measured_runs;
}

double compact_average(const CompactNSearch::NeighborhoodSearch& search)
{
	const auto& points = search.point_set(0);
	std::size_t total = 0;
	for (unsigned int i = 0; i < points.n_points(); ++i) {
		total += points.n_neighbors(0, i);
	}
	return static_cast<double>(total) / static_cast<double>(points.n_points());
}
}

int main()
{
	constexpr double particle_radius = 0.003;
	constexpr double spacing = 2.0 * particle_radius;
	constexpr double search_radius = 4.8 * particle_radius;
	std::cout << "method,grid,particles,average_neighbors,milliseconds\n";

	for (const std::size_t side : {10u, 15u, 20u}) {
		auto positions = learnSPH::regular_grid(side, side, side, spacing);
		learnSPH::NeighborList brute_result;
		const double milliseconds = average_milliseconds([&] {
			brute_result = learnSPH::brute_force_neighbors(positions, search_radius);
		}, 3);
		std::cout << "brute_force," << side << ',' << positions.size() << ','
		          << learnSPH::average_neighbor_count(brute_result) << ','
		          << std::fixed << std::setprecision(3) << milliseconds << '\n';
	}

	for (const std::size_t side : {10u, 15u, 20u, 30u, 40u, 50u}) {
		auto positions = learnSPH::regular_grid(side, side, side, spacing);
		CompactNSearch::NeighborhoodSearch search(search_radius, true);
		search.add_point_set(positions.front().data(), positions.size(), true, true, true);
		const double milliseconds = average_milliseconds([&] {
			search.find_neighbors();
		}, 5);
		std::cout << "compact_nsearch," << side << ',' << positions.size() << ','
		          << compact_average(search) << ','
		          << std::fixed << std::setprecision(3) << milliseconds << '\n';
	}
	return 0;
}
