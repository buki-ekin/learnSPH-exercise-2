#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "../learnSPH/kernel.h"
#include <CompactNSearch/CompactNSearch>
#include <Eigen/Dense>
#include <cmath>
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
        argc > 1 ? argv[1] : "output/assignment5_freefall";

    std::filesystem::create_directories(output_folder);
	// Assignment 5(a): fluid particle radius is 3 cm.
	constexpr double particle_radius = 0.03;		
	constexpr double particle_distance = 2.0 * particle_radius;
	// Rest density of water in kg/m^3.
    constexpr double rest_density = 1000.0;

   // Each fluid particle represents one cubic sampling cell.
   constexpr double particle_volume =
    particle_distance * particle_distance * particle_distance;

   // Particle mass in kg.
    // Assignment 5(a)(i): reduce nominal fluid mass by 10%
    // to limit initial pressure caused by sampling errors.
    constexpr double particle_mass =
        0.9 * rest_density * particle_volume;

    constexpr double pressure_stiffness = 1000.0;

    // Assignment 5(a)(ii): fluid kinematic viscosity in m^2/s.
    constexpr double fluid_viscosity = 0.1;
    constexpr double frame_rate = 30.0;
    constexpr int number_of_frames = 91;
		// Kernel smoothing length required by Exercise 2.
    constexpr double smoothing_length = 1.2 * particle_distance;
	constexpr double compact_support = 2.0 * smoothing_length;

	// Reuse the cubic spline kernel implemented in Exercise 1.
	const learnSPH::kernel::CubicSplineKernel kernel(smoothing_length); 
	std::vector<Eigen::Vector3d> positions;

	// Assignment 5(a): sample a one-metre cube with 6 cm spacing.
// Keep particle centres at least one radius inside the box.
const double box_length = 1.0;

for (int ix = 0; ; ++ix) {
    const double x = particle_radius + ix * particle_distance;
    if (x > box_length - particle_radius + 1e-12) {
        break;
    }

    for (int iy = 0; ; ++iy) {
        const double y = particle_radius + iy * particle_distance;
        if (y > box_length - particle_radius + 1e-12) {
            break;
        }

        for (int iz = 0; ; ++iz) {
            const double z = particle_radius + iz * particle_distance;
            if (z > box_length - particle_radius + 1e-12) {
                break;
            }

            positions.emplace_back(x, y, z);
        }
    }
}

	       // Assignment 5(a)(iv): a single fluid block starts from rest.
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
    // Assignment 5(a)(iv): gravity acts in the negative z direction.
    const Eigen::Vector3d gravity(0.0, 0.0, -9.81);

    // Allow approximately 44 metres of free fall over three seconds.
    // This is a diagnostic box, not a physical boundary.
    const Eigen::Vector3d box_min(-5.0, -5.0, -60.0);
    const Eigen::Vector3d box_max(5.0, 5.0, 5.0);
	
	    // Save initial positions for the stationary baseline check.
    const std::vector<Eigen::Vector3d> initial_positions = positions;

    double simulation_time = 0.0;

	learnSPH::write_particles_to_vtk(
		frame_path(output_folder, 0), positions, heights, velocities);

	for (int frame = 1; frame < number_of_frames; ++frame) {
		const double next_frame_time = frame / frame_rate;

		while (simulation_time + 1e-12 < next_frame_time) {
			            // Diagnose the first escaping particle before deleting anything.
            for (std::size_t i = 0; i < positions.size(); ++i) {
                const bool outside =
                    (positions[i].array() < box_min.array()).any()
                    || (positions[i].array() > box_max.array()).any();

                if (outside) {
                    double maximum_speed_now = 0.0;

                    for (const Eigen::Vector3d& velocity : velocities) {
                        maximum_speed_now = std::max(
                            maximum_speed_now, velocity.norm());
                    }

                    std::cerr
                        << "\nCollision boundary check\n"
                        << "Time [s]: " << simulation_time << '\n'
                        << "Particle index: " << i << '\n'
                        << "Position [m]: "
                        << positions[i].transpose() << '\n'
                        << "Velocity [m/s]: "
                        << velocities[i].transpose() << '\n'
                        << "Maximum speed [m/s]: "
                        << maximum_speed_now << '\n'
                        << "Box minimum [m]: "
                        << box_min.transpose() << '\n'
                        << "Box maximum [m]: "
                        << box_max.transpose() << '\n';

                    return 1;
                }
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
				velocities, particle_radius, 0.5, 0.0001, 0.001);
			const double time_step =
				std::min(cfl_step, next_frame_time - simulation_time);
                       // Assignment 5(a)(i): compute pressure from current densities.
            std::vector<double> pressures(positions.size(), 0.0);

            for (std::size_t i = 0; i < positions.size(); ++i) {
                pressures[i] = std::max(
                    0.0,
                    pressure_stiffness * (densities[i] - rest_density));
            }

            // Fluid pressure and viscosity are active; no boundaries yet.            // Gravity remains zero.
            std::vector<Eigen::Vector3d> accelerations(
                positions.size(), gravity);

            for (unsigned int i = 0; i < positions.size(); ++i) {
                const double pressure_factor_i =
                    pressures[i] / (densities[i] * densities[i]);

                for (unsigned int n = 0;
                     n < fluid_neighbors.n_neighbors(point_set, i);
                     ++n) {
                    const unsigned int j =
                        fluid_neighbors.neighbor(point_set, i, n);

                    const double pressure_factor_j =
                        pressures[j] / (densities[j] * densities[j]);

                    const Eigen::Vector3d separation =
                        positions[i] - positions[j];

                           const Eigen::Vector3d gradient =
                        kernel.gradient(separation);

                    // Pressure acceleration.
                    accelerations[i] -=
                        particle_mass
                        * (pressure_factor_i + pressure_factor_j)
                        * gradient;

                    // Assignment 5(a)(ii): fluid viscosity acceleration.
                    const double geometric_factor =
                        separation.dot(gradient)
                        / (separation.squaredNorm()
                           + 0.01 * smoothing_length * smoothing_length);

                                      accelerations[i] +=
                        2.0 * fluid_viscosity
                        * (particle_mass / densities[j])
                        * geometric_factor
                        * (velocities[i] - velocities[j]);
                } // End neighbor loop.
            } // End particle loop.

            learnSPH::semi_implicit_euler(
                positions, velocities, accelerations, time_step);

            // Detect numerical failure before the next neighbor search.
            for (std::size_t i = 0; i < positions.size(); ++i) {
                if (!positions[i].allFinite() ||
                    !velocities[i].allFinite()) {
                    std::cerr << "Non-finite particle state at t = "
                              << simulation_time + time_step
                              << " s.\n";
                    return 1;
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
    // Assignment 5(a)(iii): report particle count, displacement and speed.
    if (positions.size() != initial_positions.size()) {
        std::cerr << "FAIL: particle count changed.\n";
        return 1;
    }

    double maximum_displacement = 0.0;
    double maximum_speed = 0.0;

    for (std::size_t i = 0; i < positions.size(); ++i) {
        maximum_displacement = std::max(
            maximum_displacement,
            (positions[i] - initial_positions[i]).norm());

        maximum_speed = std::max(
            maximum_speed,
            velocities[i].norm());
    }

    std::cout << "Maximum displacement [m]: "
              << maximum_displacement << '\n';
    std::cout << "Maximum speed [m/s]: "
              << maximum_speed << '\n';

    return 0;

	
}
