#pragma once
#include <Eigen/Dense>

namespace learnSPH
{
	namespace kernel
	{
		constexpr double PI = 3.14159265358979323846;

		// Dimensionless cubic spline s(q) and its derivative ds/dq.
		double cubic_spline(const double q);
		double cubic_grad_spline(const double q);

		// Full three-dimensional kernel W(x_i-x_j; h) and its gradient
		// with respect to x_i.
		double cubic_spline(const Eigen::Vector3d& x_i,
		                    const Eigen::Vector3d& x_j,
		                    double h);
		Eigen::Vector3d cubic_spline_gradient(const Eigen::Vector3d& x_i,
		                                      const Eigen::Vector3d& x_j,
		                                      double h);

		// Optimized version for repeated evaluations with constant h. All
		// powers and normalization factors involving h are precomputed once.
		class CubicSplineKernel
		{
		public:
			explicit CubicSplineKernel(double smoothing_length);

			double value(const Eigen::Vector3d& difference) const;
			Eigen::Vector3d gradient(const Eigen::Vector3d& difference) const;
			double smoothing_length() const { return h_; }
			double support_radius() const { return 2.0 * h_; }

		private:
			double h_;
			double inv_h_;
			double inv_h3_;
			double inv_h4_;
		};
	};
};
