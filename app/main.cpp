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
        argc > 1 ? argv[1] : "output/assignment5_dambreak_reference";

    std::filesystem::create_directories(output_folder);
		// Assignment 5(a): fluid particle radius is 3 cm.
    // Assignment 5(b): particle radius is 5 mm.
    constexpr double particle_radius = 0.005;
	constexpr double particle_distance = 2.0 * particle_radius;
	// Rest density of water in kg/m^3.
    constexpr double rest_density = 1000.0;

   // Each fluid particle represents one cubic sampling cell.
   constexpr double particle_volume =
    particle_distance * particle_distance * particle_distance;

   // Particle mass in kg.
    // Assignment 5(b): nominal mass from rest density and sampling volume.
    constexpr double particle_mass =
        rest_density * particle_volume;

    constexpr double pressure_stiffness = 1000.0;

    // Assignment 5(a)(ii): fluid kinematic viscosity in m^2/s.
    constexpr double fluid_viscosity = 0.0025;
	 // Assignment 5(a)(vi): friction with the stationary boundary.
    constexpr double boundary_viscosity = 0.0;
    // Assignment 5(b): export every 0.02 s to include exactly t = 0.38 s.
    constexpr double frame_rate = 50.0;    // Inspect the initial floor impact over one second.
    // Run through the observed velocity spike: 0.4 seconds.
    // Frames 0 through 19 cover t = 0 through 0.38 s.
    constexpr int number_of_frames = 20;			// Kernel smoothing length required by Exercise 2.
    constexpr double smoothing_length = 1.2 * particle_distance;
	constexpr double compact_support = 2.0 * smoothing_length;

	// Reuse the cubic spline kernel implemented in Exercise 1.
	const learnSPH::kernel::CubicSplineKernel kernel(smoothing_length); 
	std::vector<Eigen::Vector3d> positions;

	// Assignment 5(a): sample a one-metre cube with 6 cm spacing.
// Keep particle centres at least one radius inside the box.
    // Assignment 5(b): fluid dimensions are 0.15 x 0.25 x 0.50 m.
    // Place particles at the centres of 1 cm sampling cells.
    constexpr int fluid_count_x = 15;
    constexpr int fluid_count_y = 25;
    constexpr int fluid_count_z = 50;

    for (int ix = 0; ix < fluid_count_x; ++ix) {
        for (int iy = 0; iy < fluid_count_y; ++iy) {
            for (int iz = 0; iz < fluid_count_z; ++iz) {
                positions.emplace_back(
                    (ix + 0.5) * particle_distance,
                    (iy + 0.5) * particle_distance,
                    (iz + 0.5) * particle_distance);
            }
        }
    }

    std::cout << "Initial fluid particles: "
              << positions.size() << '\n';

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
    // Assignment 5(b): open-top tank.
    // Nominal dimensions: 0.15 x 0.80 x 1.00 m.
    // Offset each boundary outwards by one particle diameter.
    const double boundary_gap = particle_distance;
    const double boundary_spacing = 0.8 * particle_distance;

    const double tank_x_min = -boundary_gap;
    const double tank_x_max = 0.15 + boundary_gap;

    const double tank_y_min = -boundary_gap;
    const double tank_y_max = 0.80 + boundary_gap;

    const double floor_z = -boundary_gap;
    const double tank_z_max = 1.00 + boundary_gap;

    // Bottom vertices: 0, 1, 2, 3.
    // Corresponding upper vertices: 4, 5, 6, 7.
    const std::vector<Eigen::Vector3d> boundary_vertices = {
        Eigen::Vector3d(tank_x_min, tank_y_min, floor_z),
        Eigen::Vector3d(tank_x_max, tank_y_min, floor_z),
        Eigen::Vector3d(tank_x_max, tank_y_max, floor_z),
        Eigen::Vector3d(tank_x_min, tank_y_max, floor_z),

        Eigen::Vector3d(tank_x_min, tank_y_min, tank_z_max),
        Eigen::Vector3d(tank_x_max, tank_y_min, tank_z_max),
        Eigen::Vector3d(tank_x_max, tank_y_max, tank_z_max),
        Eigen::Vector3d(tank_x_min, tank_y_max, tank_z_max)
    };

    // Two triangles per face; no top face.
    const std::vector<std::array<int, 3>> boundary_triangles = {
        // Bottom.
        {0, 1, 2}, {0, 2, 3},

        // Wall at minimum x.
        {0, 3, 7}, {0, 7, 4},

        // Wall at maximum x.
        {1, 5, 6}, {1, 6, 2},

        // Wall at minimum y.
        {0, 4, 5}, {0, 5, 1},

        // Wall at maximum y.
        {3, 2, 6}, {3, 6, 7}
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

    std::cout << "Boundary particles: "              << boundary_positions.size() << '\n';
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
                // Diagnose the velocity spike with a smaller maximum time step.
    const double cfl_step = learnSPH::cfl_time_step(
        velocities, particle_radius, 0.5, 0.000001, 0.00005);
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
            // Assignment 5(b): diagnose penetration through the tank.
            // The top is open, so rising above the rim is allowed.
            const double penetration_tolerance = particle_radius;

            for (std::size_t i = 0; i < positions.size(); ++i) {
                const Eigen::Vector3d& p = positions[i];

                // Exclude the rim region from side-wall diagnostics.
                const bool below_rim =
                    p.z() < tank_z_max - compact_support;

                const bool outside_side =
                    p.x() < tank_x_min - penetration_tolerance ||
                    p.x() > tank_x_max + penetration_tolerance ||
                    p.y() < tank_y_min - penetration_tolerance ||
                    p.y() > tank_y_max + penetration_tolerance;

                const bool inside_floor_footprint =
                    p.x() >= tank_x_min &&
                    p.x() <= tank_x_max &&
                    p.y() >= tank_y_min &&
                    p.y() <= tank_y_max;

                const bool below_floor =
                    p.z() < floor_z - penetration_tolerance;

                if ((below_rim && outside_side) ||
                    (inside_floor_footprint && below_floor)) {
                    std::cerr
                        << "\nTank escape or penetration detected\n"
                        << "Time [s]: "
                        << simulation_time + time_step << '\n'
                        << "Particle index: " << i << '\n'
                        << "Position [m]: "
                        << p.transpose() << '\n'
                        << "Velocity [m/s]: "
                        << velocities[i].transpose() << '\n';

                    return 1;
                }
            }

            simulation_time += time_step;		}

        update_height(positions, heights);

        // Densities belong to the last force evaluation before this output.
        const auto maximum_density_it =
            std::max_element(densities.begin(), densities.end());

        const double maximum_density_now = *maximum_density_it;

        const double maximum_pressure_now = std::max(
            0.0,
            pressure_stiffness * (maximum_density_now - rest_density));

        // Locate the fastest particle in the current output frame.
        std::size_t fastest_particle = 0;

        for (std::size_t i = 1; i < velocities.size(); ++i) {
            if (velocities[i].squaredNorm() >
                velocities[fastest_particle].squaredNorm()) {
                fastest_particle = i;
            }
        }

        std::cout
            << "Density/pressure check at t=" << simulation_time
            << "  max_rho=" << maximum_density_now
            << "  max_pressure=" << maximum_pressure_now
            << '\n'
            << "Fastest particle position [m]: "
            << positions[fastest_particle].transpose()
            << '\n'
            << "Fastest particle velocity [m/s]: "
            << velocities[fastest_particle].transpose()
            << '\n';
        // Track when strong upward motion first appears.
        double frame_max_z = positions.front().z();
        double frame_max_speed = 0.0;
        std::size_t frame_above_rim = 0;

        for (std::size_t i = 0; i < positions.size(); ++i) {
            frame_max_z = std::max(frame_max_z, positions[i].z());
            frame_max_speed = std::max(
                frame_max_speed, velocities[i].norm());

            if (positions[i].z() > tank_z_max) {
                ++frame_above_rim;
            }
        }

        std::cout << "t=" << simulation_time
                  << "  max_z=" << frame_max_z
                  << "  max_speed=" << frame_max_speed
                  << "  above_rim=" << frame_above_rim
                  << '\n';

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

    // Diagnose the final spatial extent of the dam-break simulation.
    Eigen::Vector3d minimum_position = positions.front();
    Eigen::Vector3d maximum_position = positions.front();

    std::size_t particles_above_rim = 0;
    std::size_t particles_outside_footprint = 0;

    for (const Eigen::Vector3d& p : positions) {
        minimum_position = minimum_position.cwiseMin(p);
        maximum_position = maximum_position.cwiseMax(p);

        if (p.z() > tank_z_max) {
            ++particles_above_rim;
        }

        if (p.x() < tank_x_min || p.x() > tank_x_max ||
            p.y() < tank_y_min || p.y() > tank_y_max) {
            ++particles_outside_footprint;
        }
    }

    std::cout << "\nDam-break extent check\n";
    std::cout << "Minimum position [x y z, m]: "
              << minimum_position.transpose() << '\n';
    std::cout << "Maximum position [x y z, m]: "
              << maximum_position.transpose() << '\n';
    std::cout << "Particles above the rim: "
              << particles_above_rim << '\n';
    std::cout << "Particles outside the tank footprint: "
              << particles_outside_footprint << '\n';

    return 0;

	
}
