#include <algorithm>
#include <array>
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
// Compute effective volumes using boundary neighbors, including self.
std::vector<double> compute_boundary_volumes(
    std::vector<Eigen::Vector3d>& boundary_positions,
    const double smoothing_length)
{
    std::vector<double> volumes(boundary_positions.size(), 0.0);

    if (boundary_positions.empty()) {
        return volumes;
    }

    const learnSPH::kernel::CubicSplineKernel kernel(smoothing_length);

    CompactNSearch::NeighborhoodSearch search(
        2.0 * smoothing_length);

    const unsigned int boundary_set = search.add_point_set(
        boundary_positions.front().data(),
        static_cast<unsigned int>(boundary_positions.size()),
        false, true, true);

    search.find_neighbors();

    const auto& neighbors = search.point_set(boundary_set);

    for (unsigned int k = 0; k < boundary_positions.size(); ++k) {
        // CompactNSearch excludes self, so add its contribution explicitly.
        double kernel_sum = kernel.value(Eigen::Vector3d::Zero());

        for (unsigned int n = 0;
             n < neighbors.n_neighbors(boundary_set, k);
             ++n) {
            const unsigned int l =
                neighbors.neighbor(boundary_set, k, n);

            const Eigen::Vector3d separation =
                boundary_positions[k] - boundary_positions[l];

            kernel_sum += kernel.value(separation);
        }

        volumes[k] = 1.0 / kernel_sum;
    }

    return volumes;
}

// Convert effective boundary volumes to masses using fluid rest density.
std::vector<double> compute_boundary_masses(
    const std::vector<double>& boundary_volumes,
    const double rest_density)
{
    std::vector<double> masses(boundary_volumes.size());

    for (std::size_t k = 0; k < boundary_volumes.size(); ++k) {
        masses[k] = rest_density * boundary_volumes[k];
    }

    return masses;
}
int main(int argc, char** argv)
{
            const std::filesystem::path output_folder =
        argc > 1 ? argv[1] : "output/assignment5_floor_friction";

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
	 // Assignment 5(a)(vi): friction with the stationary boundary.
    constexpr double boundary_viscosity = 1.0;
    constexpr double frame_rate = 30.0;
    // Inspect the initial floor impact over one second.
    constexpr int number_of_frames = 31;
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
    // Assignment 5(a)(v): sample a static floor below the fluid block.
    // The initial gap avoids placing fluid directly against the boundary.
    const double floor_z = -0.30;
    const double boundary_spacing = 0.8 * particle_distance;

    const std::vector<Eigen::Vector3d> boundary_vertices = {
        Eigen::Vector3d(-2.0, -2.0, floor_z),
        Eigen::Vector3d( 3.0, -2.0, floor_z),
        Eigen::Vector3d( 3.0,  3.0, floor_z),
        Eigen::Vector3d(-2.0,  3.0, floor_z)
    };

    const std::vector<std::array<int, 3>> boundary_triangles = {
        {0, 1, 2},
        {0, 2, 3}
    };

    std::vector<Eigen::Vector3d> boundary_positions;

    learnSPH::sampling::triangle_mesh(
        boundary_positions,
        boundary_vertices,
        boundary_triangles,
        boundary_spacing);

    if (boundary_positions.empty()) {
        std::cerr << "Floor sampling produced no particles.\n";
        return 1;
    }

    // Static boundary volumes and masses are computed only once.
    const std::vector<double> boundary_volumes =
        compute_boundary_volumes(
            boundary_positions, smoothing_length);

    const std::vector<double> boundary_masses =
        compute_boundary_masses(
            boundary_volumes, rest_density);

    // Register the floor in the same search as the fluid.
    const unsigned int wall_set = search.add_point_set(
        boundary_positions.front().data(),
        static_cast<unsigned int>(boundary_positions.size()),
        false, true, true);

    std::cout << "Floor particles: "
              << boundary_positions.size() << '\n';
    std::cout << "Floor height [m]: "
              << floor_z << '\n';
    std::cout << "Boundary spacing [m]: "
              << boundary_spacing << '\n';

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
                // Fluid contribution, including the self contribution.
                densities[i] = particle_mass * kernel_sum;

                // Assignment 5(a)(v): add the boundary density contribution.
                for (unsigned int n = 0;
                     n < fluid_neighbors.n_neighbors(wall_set, i);
                     ++n) {
                    const unsigned int k =
                        fluid_neighbors.neighbor(wall_set, i, n);

                    const Eigen::Vector3d separation =
                        positions[i] - boundary_positions[k];

                    densities[i] +=
                        boundary_masses[k] * kernel.value(separation);
                }		}	

			            // Compare floor impact using a smaller maximum time step.
            const double cfl_step = learnSPH::cfl_time_step(
                velocities, particle_radius, 0.5, 0.00001, 0.00025);
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
                               } // End fluid neighbor loop.

                // Assignment 5(a)(v): pressure interaction with the floor.
                               // Boundary pressure and viscosity are active.
                			for (unsigned int n = 0;
                     			n < fluid_neighbors.n_neighbors(wall_set, i);
                     ++n) { 
                    const unsigned int k =
                        fluid_neighbors.neighbor(wall_set, i, n);

                    const Eigen::Vector3d separation =
                        positions[i] - boundary_positions[k];

                                       const Eigen::Vector3d gradient =
                        kernel.gradient(separation);

                    // Boundary pressure contribution.
                    accelerations[i] -=
                        boundary_masses[k]
                        * pressure_factor_i
                        * gradient;

                    // Assignment 5(a)(vi): viscosity at a stationary wall.
                    const double boundary_geometric_factor =
                        separation.dot(gradient)
                        / (separation.squaredNorm()
                           + 0.01 * smoothing_length * smoothing_length);

                    accelerations[i] +=
                        2.0 * boundary_viscosity
                        * boundary_volumes[k]
                        * boundary_geometric_factor
                        * velocities[i];
                }
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
            // Detect penetration beneath the interior of the floor.
            // Exclude an edge margin where particles can leave the floor.
            for (std::size_t i = 0; i < positions.size(); ++i) {
                const bool above_floor_interior =
                    positions[i].x() > -2.0 + compact_support &&
                    positions[i].x() <  3.0 - compact_support &&
                    positions[i].y() > -2.0 + compact_support &&
                    positions[i].y() <  3.0 - compact_support;

                if (above_floor_interior &&
                    positions[i].z() < floor_z - particle_radius) {
                    std::cerr
                        << "\nFloor penetration check\n"
                        << "Time [s]: "
                        << simulation_time + time_step << '\n'
                        << "Particle index: " << i << '\n'
                        << "Position [m]: "
                        << positions[i].transpose() << '\n'
                        << "Velocity [m/s]: "
                        << velocities[i].transpose() << '\n';

                    return 1;
                }
            }

            simulation_time += time_step;		}

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
