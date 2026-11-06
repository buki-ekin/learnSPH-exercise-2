#include <iomanip>
#include <string>
#include <array>
#include "../learnSPH/sampling.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>
#include <filesystem>
#include "../learnSPH/io.h"
#include "../learnSPH/kernel.h"

// Assignment 2(a): compute effective boundary particle volumes.
// Sum over boundary neighbors only, including the particle itself.
// For static boundaries, these volumes only need to be computed once.
std::vector<double> compute_boundary_volumes(
    std::vector<Eigen::Vector3d>& boundary_positions,
    const double smoothing_length)
{
    std::vector<double> volumes(boundary_positions.size(), 0.0);

    // An empty boundary set has no volumes to compute.
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
        // CompactNSearch excludes self, so add W(0) explicitly.
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

        // A 3D kernel has units of 1/m^3.
        // The reciprocal of its sum therefore gives a volume in m^3.
        volumes[k] = 1.0 / kernel_sum;
    }

    return volumes;
}

// Assignment 2(a): convert boundary volumes to effective masses.
// Use the same rest density as the fluid.
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

int main()
{
    // Assignment 1(b): sample a one-litre cube with 20^3 particles.
    constexpr int particles_per_axis = 20;
    constexpr double box_length = 0.1; // metres
    constexpr double rest_density = 1000.0; // kg/m^3

    constexpr double particle_distance =
        box_length / particles_per_axis;

    constexpr double particle_radius = particle_distance / 2.0;
    constexpr double smoothing_length = 1.2 * particle_distance;
    constexpr double compact_support = 2.0 * smoothing_length;

    constexpr double particle_mass =
        rest_density * particle_distance
                     * particle_distance
                     * particle_distance;

    std::vector<Eigen::Vector3d> positions;

    // Place each particle at the centre of its sampling cell.
    for (int i = 0; i < particles_per_axis; ++i) {
        for (int j = 0; j < particles_per_axis; ++j) {
            for (int k = 0; k < particles_per_axis; ++k) {
                positions.emplace_back(
                    (i + 0.5) * particle_distance,
                    (j + 0.5) * particle_distance,
                    (k + 0.5) * particle_distance);
            }
        }
    }

    std::cout << "Particle count: " << positions.size() << '\n';
    std::cout << "Particle radius [m]: " << particle_radius << '\n';
    std::cout << "Particle mass [kg]: " << particle_mass << '\n';
    std::cout << "Total mass [kg]: "
              << positions.size() * particle_mass << '\n';
    // Assignment 1(b): find neighbors within the kernel support.
    const learnSPH::kernel::CubicSplineKernel kernel(smoothing_length);

    CompactNSearch::NeighborhoodSearch search(compact_support);

    const unsigned int point_set = search.add_point_set(
        positions.front().data(),
        static_cast<unsigned int>(positions.size()),
        true, true, true);

    search.find_neighbors();

    const auto& fluid_neighbors = search.point_set(point_set);
    std::vector<double> densities(positions.size(), 0.0);

    // Compute density, including each particle's own contribution.
    for (unsigned int i = 0; i < positions.size(); ++i) {
        double kernel_sum = kernel.value(Eigen::Vector3d::Zero());

        for (unsigned int n = 0;
             n < fluid_neighbors.n_neighbors(point_set, i);
             ++n) {
            const unsigned int j =
                fluid_neighbors.neighbor(point_set, i, n);

            const Eigen::Vector3d separation =
                positions[i] - positions[j];

            kernel_sum += kernel.value(separation);
        }

        densities[i] = particle_mass * kernel_sum;
    }

    // Summarize the computed density distribution.
    double minimum_density = densities.front();
    double maximum_density = densities.front();
    double density_sum = 0.0;
    unsigned int near_rest_density = 0;

    for (const double density : densities) {
        if (density < minimum_density) {
            minimum_density = density;
        }
        if (density > maximum_density) {
            maximum_density = density;
        }

        density_sum += density;

        // Count particles within 5 percent of the rest density.
        if (density >= 0.95 * rest_density &&
            density <= 1.05 * rest_density) {
            ++near_rest_density;
        }
    }

    std::cout << "Minimum density [kg/m^3]: "
              << minimum_density << '\n';
    std::cout << "Maximum density [kg/m^3]: "
              << maximum_density << '\n';
    std::cout << "Mean density [kg/m^3]: "
              << density_sum / densities.size() << '\n';
    std::cout << "Particles within 5% of rest density: "
              << near_rest_density << " / "
                << densities.size() << '\n';
        // Assignment 1(b): group densities into 50 kg/m^3 intervals.
    constexpr double bin_width = 50.0;

    const std::size_t bin_count =
        static_cast<std::size_t>(maximum_density / bin_width) + 1;

    std::vector<unsigned int> histogram(bin_count, 0);

    for (const double density : densities) {
        const std::size_t bin =
            static_cast<std::size_t>(density / bin_width);

        ++histogram[bin];
    }

    // Print exact counts and a compact text visualization.
    // Each star represents approximately 100 particles.
    std::cout << "\nDensity histogram [kg/m^3]\n";
    std::cout << "Intervals include the lower bound, exclude the upper.\n";
    std::cout << "One * represents approximately 100 particles.\n\n";

    for (std::size_t bin = 0; bin < histogram.size(); ++bin) {
        const double lower = bin * bin_width;
        const double upper = lower + bin_width;
        const unsigned int count = histogram[bin];

        // Keep small non-empty bins visible with at least one star.
        const std::size_t stars =
            count == 0 ? 0 : (count + 99) / 100;

        std::cout << "[" << std::setw(4) << lower
                  << ", " << std::setw(4) << upper << ") : "
                  << std::setw(4) << count << "  "
                  << std::string(stars, '*') << '\n';
    }
    
        // Assignment 2(b): represent a horizontal floor with two triangles.
    // Boundary spacing is 0.8 times the fluid particle diameter.
    const double boundary_spacing = 0.8 * particle_distance;
    const double margin = compact_support;

    // The lowest fluid centres are at z = particle_distance / 2.
    // This places the floor one fluid diameter below those centres.
    const double floor_z = -0.5 * particle_distance;

    const std::vector<Eigen::Vector3d> boundary_vertices = {
        Eigen::Vector3d(-margin, -margin, floor_z),
        Eigen::Vector3d(box_length + margin, -margin, floor_z),
        Eigen::Vector3d(
            box_length + margin, box_length + margin, floor_z),
        Eigen::Vector3d(-margin, box_length + margin, floor_z)
    };

    const std::vector<std::array<int, 3>> boundary_triangles = {
        {0, 1, 2},
        {0, 2, 3}
    };

    // Sample the triangle surface using the supplied framework.
    std::vector<Eigen::Vector3d> boundary_positions;

    learnSPH::sampling::triangle_mesh(
        boundary_positions,
        boundary_vertices,
        boundary_triangles,
        boundary_spacing);

    // Store effective volumes and masses for the static boundary.
    const std::vector<double> boundary_volumes =
        compute_boundary_volumes(
            boundary_positions, smoothing_length);

    const std::vector<double> boundary_masses =
        compute_boundary_masses(
            boundary_volumes, rest_density);

    std::cout << "\nAssignment 2(b): boundary sampling\n";
    std::cout << "Boundary particle count: "
              << boundary_positions.size() << '\n';
    std::cout << "Boundary spacing [m]: "
              << boundary_spacing << '\n';

    // Check that sampling produced particles before accessing values.
    if (boundary_volumes.empty()) {
        std::cerr << "Error: boundary sampling produced no particles.\n";
        return 1;
    }

    double minimum_volume = boundary_volumes.front();
    double maximum_volume = boundary_volumes.front();
    double total_boundary_mass = 0.0;

    for (std::size_t k = 0; k < boundary_volumes.size(); ++k) {
        if (boundary_volumes[k] < minimum_volume) {
            minimum_volume = boundary_volumes[k];
        }
        if (boundary_volumes[k] > maximum_volume) {
            maximum_volume = boundary_volumes[k];
        }

        total_boundary_mass += boundary_masses[k];
    }

    std::cout << "Minimum boundary volume [m^3]: "
              << minimum_volume << '\n';
    std::cout << "Maximum boundary volume [m^3]: "
              << maximum_volume << '\n';
    std::cout << "Total effective boundary mass [kg]: "
              << total_boundary_mass << '\n';
    
        // Assignment 3(a): find boundary neighbors for each fluid particle.
    CompactNSearch::NeighborhoodSearch coupled_search(compact_support);

    const unsigned int fluid_set = coupled_search.add_point_set(
        positions.front().data(),
        static_cast<unsigned int>(positions.size()),
        true, true, true);

    const unsigned int wall_set = coupled_search.add_point_set(
        boundary_positions.front().data(),
        static_cast<unsigned int>(boundary_positions.size()),
        false, true, true);

    coupled_search.find_neighbors();

    const auto& fluid_point_set = coupled_search.point_set(fluid_set);

    // Start from the fluid-only densities calculated in Assignment 1.
    // Keep the original values unchanged for comparison.
    std::vector<double> densities_with_boundary = densities;

    unsigned int affected_particles = 0;
    double maximum_boundary_contribution = 0.0;

    for (unsigned int i = 0; i < positions.size(); ++i) {
        double boundary_contribution = 0.0;

        for (unsigned int n = 0;
             n < fluid_point_set.n_neighbors(wall_set, i);
             ++n) {
            const unsigned int k =
                fluid_point_set.neighbor(wall_set, i, n);

            const Eigen::Vector3d separation =
                positions[i] - boundary_positions[k];

            // Effective boundary mass is rest_density * boundary_volume.
            boundary_contribution +=
                boundary_masses[k] * kernel.value(separation);
        }

        densities_with_boundary[i] += boundary_contribution;

        if (boundary_contribution > 0.0) {
            ++affected_particles;
        }

        maximum_boundary_contribution = std::max(
            maximum_boundary_contribution, boundary_contribution);
    }

    // Report the density distribution after adding the boundary.
    double minimum_with_boundary = densities_with_boundary.front();
    double maximum_with_boundary = densities_with_boundary.front();
    double sum_with_boundary = 0.0;

    for (const double density : densities_with_boundary) {
        minimum_with_boundary = std::min(minimum_with_boundary, density);
        maximum_with_boundary = std::max(maximum_with_boundary, density);
        sum_with_boundary += density;
    }

    std::cout << "\nAssignment 3(a): density with boundary contribution\n";
    std::cout << "Particles receiving a boundary contribution: "
              << affected_particles << '\n';
    std::cout << "Maximum boundary contribution [kg/m^3]: "
              << maximum_boundary_contribution << '\n';
    std::cout << "Minimum density [kg/m^3]: "
              << minimum_with_boundary << '\n';
    std::cout << "Maximum density [kg/m^3]: "
              << maximum_with_boundary << '\n';
    std::cout << "Mean density [kg/m^3]: "
              << sum_with_boundary / densities_with_boundary.size()
              << '\n';  
        // Assignment 3(b): export fluid densities for visualization.
    const std::filesystem::path output_folder =
        "output/assignment3";

    std::filesystem::create_directories(output_folder);

    // Export the fluid-only estimate as a reference.
    learnSPH::write_particles_to_vtk(
        (output_folder / "fluid_density_without_boundary.vtk").string(),
        positions,
        densities);

    // Export the density estimate including the floor contribution.
    learnSPH::write_particles_to_vtk(
        (output_folder / "fluid_density_with_boundary.vtk").string(),
        positions,
        densities_with_boundary);

    // Export boundary positions and their corrected volumes.
    learnSPH::write_particles_to_vtk(
        (output_folder / "boundary_volumes.vtk").string(),
        boundary_positions,
        boundary_volumes);

    std::cout << "\nAssignment 3(b): VTK files saved to "
              << output_folder.string() << '\n';                  

return 0;

}