#include <iomanip>
#include <string>

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/kernel.h"

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
    return 0;
    
}