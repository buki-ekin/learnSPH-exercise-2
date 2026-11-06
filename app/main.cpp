#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "../learnSPH/kernel.h"
#include <CompactNSearch/CompactNSearch>
#include <Eigen/Dense>

#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/simulation.h"

namespace
{
std::string frame_path(const std::filesystem::path& folder, const int frame)
{
	std::ostringstream name;
	name << "exercise1_" << std::setfill('0') << std::setw(4) << frame << ".vtk";
	return (folder / name.str()).string();
}

void update_height(const std::vector<Eigen::Vector3d>& positions,
	               std::vector<double>& heights)
{
	for (std::size_t i = 0; i < positions.size(); ++i) {
		heights[i] = positions[i].z();
	}
}
}

int main(int argc, char** argv)
{
	const std::filesystem::path output_folder =
		argc > 1 ? argv[1] : "output/vtk";
	std::filesystem::create_directories(output_folder);

	constexpr double particle_radius = 0.05;
	constexpr double particle_distance = 2.0 * particle_radius;
	// Rest density of water in kg/m^3.
    constexpr double rest_density = 1000.0;

   // Each fluid particle represents one cubic sampling cell.
   constexpr double particle_volume =
    particle_distance * particle_distance * particle_distance;

   // Particle mass in kg.
    constexpr double particle_mass = rest_density * particle_volume;
	constexpr double frame_rate = 30.0;
	constexpr int number_of_frames = 91;
	// Kernel smoothing length required by Exercise 2.
    constexpr double smoothing_length = 1.2 * particle_distance;
	constexpr double compact_support = 2.0 * smoothing_length;

	// Reuse the cubic spline kernel implemented in Exercise 1.
	const learnSPH::kernel::CubicSplineKernel kernel(smoothing_length); 
	std::vector<Eigen::Vector3d> positions;

	// Fill the whole cube with particles. A two-metre cube gives us a
	// denser-looking block while keeping the exercise's 5 cm radius.
	learnSPH::sampling::fluid_box(
		positions,
		Eigen::Vector3d(-1.0, -1.0, 9.0),
		Eigen::Vector3d(1.0, 1.0, 11.0),
		particle_distance);

	std::vector<Eigen::Vector3d> velocities(
		positions.size(), Eigen::Vector3d::Zero());
	std::vector<double> heights(positions.size());
	// Store the estimated density of each fluid particle in kg/m^3.
	std::vector<double> densities(positions.size(), 0.0);
	update_height(positions, heights);

	// CompactNSearch stores the particle positions.
	// Search for neighbors inside the kernel's support radius.
	CompactNSearch::NeighborhoodSearch search(compact_support);	const unsigned int point_set = search.add_point_set(
		positions.front().data(), positions.size(), true, true, true);

	const Eigen::Vector3d gravity(0.0, 0.0, -9.81);
	const Eigen::Vector3d box_min(-5.0, -5.0, -1.0);
	const Eigen::Vector3d box_max(5.0, 5.0, 12.0);
	double simulation_time = 0.0;

	learnSPH::write_particles_to_vtk(
		frame_path(output_folder, 0), positions, heights, velocities);

	for (int frame = 1; frame < number_of_frames; ++frame) {
		const double next_frame_time = frame / frame_rate;

		while (simulation_time + 1e-12 < next_frame_time) {
			// Delete particles outside the safety box.
			const std::size_t removed = learnSPH::remove_outside_aabb(
				positions, velocities, heights, box_min, box_max);
			if (removed > 0) {
				search.resize_point_set(
					point_set, positions.front().data(), positions.size());
			}

			search.find_neighbors();
			// Assignment 1(a): estimate fluid density from neighboring particles.
			densities.resize(positions.size());

			const auto& fluid_neighbors = search.point_set(point_set);

			for (unsigned int i = 0; i < positions.size(); ++i) {

    			// CompactNSearch excludes self, so include the self contribution.
    			double kernel_sum = kernel.value(Eigen::Vector3d::Zero());

    			// Add the contributions of all other fluid neighbors.
   				 for (unsigned int n = 0;
         			n < fluid_neighbors.n_neighbors(point_set, i);
         			++n) {

        			const unsigned int j =
            			fluid_neighbors.neighbor(point_set, i, n);

        			const Eigen::Vector3d separation =
            			positions[i] - positions[j];

       				 kernel_sum += kernel.value(separation);
  				 }

   				 // All fluid particles have the same mass.
   			 densities[i] = particle_mass * kernel_sum;
		}	

			// CFL chooses a safe time step.
			const double cfl_step = learnSPH::cfl_time_step(
				velocities, particle_radius, 0.5, 0.0001, 0.01);
			const double time_step =
				std::min(cfl_step, next_frame_time - simulation_time);

			// Gravity pulls every moving particle down.
			std::vector<Eigen::Vector3d> accelerations(
				positions.size(), gravity);
			for (std::size_t i = 0; i < positions.size(); ++i) {
				if (positions[i].z() <= particle_radius &&
				    velocities[i].z() == 0.0) {
					accelerations[i].setZero();
				}
			}

			learnSPH::semi_implicit_euler(
				positions, velocities, accelerations, time_step);

			// Stop a particle when it touches the floor.
			for (std::size_t i = 0; i < positions.size(); ++i) {
				if (positions[i].z() <= particle_radius) {
					positions[i].z() = particle_radius;
					velocities[i].setZero();
				}
			}

			simulation_time += time_step;
		}

		update_height(positions, heights);
		learnSPH::write_particles_to_vtk(
			frame_path(output_folder, frame), positions, heights, velocities);
	}

	std::cout << "Particles: " << positions.size() << '\n';
	std::cout << "Frames: " << number_of_frames << '\n';
	std::cout << "Output: " << output_folder << '\n';
	return 0;
}
