#include "kernel.h"

#include <cassert>
#include <stdexcept>

double learnSPH::kernel::cubic_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	if (q < 1.0) {
		return alpha * (2.0/3.0 - q*q + 0.5*q*q*q);
	}
	else if (q < 2.0) {
		const double two_minus_q = 2.0 - q;
		return alpha * (two_minus_q * two_minus_q * two_minus_q) / 6.0;
	}
	else {
		return 0.0;
	}
}

double learnSPH::kernel::cubic_grad_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	if (q < 1.0) {
		return alpha * (-2.0*q + 1.5*q*q);
	}
	if (q < 2.0) {
		const double two_minus_q = 2.0 - q;
		return -0.5 * alpha * two_minus_q * two_minus_q;
	}
	return 0.0;
}

double learnSPH::kernel::cubic_spline(const Eigen::Vector3d& x_i,
	                                  const Eigen::Vector3d& x_j,
	                                  const double h)
{
	return CubicSplineKernel(h).value(x_i - x_j);
}

Eigen::Vector3d learnSPH::kernel::cubic_spline_gradient(
	const Eigen::Vector3d& x_i,
	const Eigen::Vector3d& x_j,
	const double h)
{
	return CubicSplineKernel(h).gradient(x_i - x_j);
}

learnSPH::kernel::CubicSplineKernel::CubicSplineKernel(const double smoothing_length)
	: h_(smoothing_length),
	  inv_h_(0.0),
	  inv_h3_(0.0),
	  inv_h4_(0.0)
{
	if (!(smoothing_length > 0.0)) {
		throw std::invalid_argument("The smoothing length h must be positive.");
	}
	inv_h_ = 1.0 / smoothing_length;
	inv_h3_ = inv_h_ * inv_h_ * inv_h_;
	inv_h4_ = inv_h3_ * inv_h_;
}

double learnSPH::kernel::CubicSplineKernel::value(
	const Eigen::Vector3d& difference) const
{
	const double q = difference.norm() * inv_h_;
	return inv_h3_ * cubic_spline(q);
}

Eigen::Vector3d learnSPH::kernel::CubicSplineKernel::gradient(
	const Eigen::Vector3d& difference) const
{
	const double distance = difference.norm();
	if (distance == 0.0) {
		return Eigen::Vector3d::Zero();
	}

	const double q = distance * inv_h_;
	if (q >= 2.0) {
		return Eigen::Vector3d::Zero();
	}

	return inv_h4_ * cubic_grad_spline(q) * difference / distance;
}

