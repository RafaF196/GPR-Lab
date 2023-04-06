#include <iostream>
#include <Eigen/Dense>
#include "MongePatch.h"


// Given a point P, its normal, and its closest neighbors (including itself) 
// compute a quadratic Monge patch that approximates the neighborhood of P.
// The resulting patch will be used to compute the principal curvatures of the 
// surface at point P.

void MongePatch::init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest)
{
	// Defining the coordinate system aligned to the normal and P
	glm::vec3 w, u, v;
	w = -normal;
	u = glm::vec3(1, 0, 0) * w;
	v = glm::cross(w, u);
	v = glm::normalize(v);

	// Transform all neighbors to this system
	std::vector<glm::vec3> trans_neigh;
	glm::vec3 aux;
	for (glm::vec3 p_i : closest) {
		aux = glm::vec3(glm::dot(u, p_i - P), glm::dot(v, p_i - P), glm::dot(w, p_i - P));
		trans_neigh.push_back(aux);
	}

	// Fitting the function using least squares
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
	Eigen::VectorXd s = Eigen::VectorXd::Zero(6);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(6);

	Eigen::VectorXd q_i(6);
	for (glm::vec3 p : trans_neigh) {
		q_i << p.x*p.x, p.x*p.y, p.y*p.y, p.x, p.y, 1;
		A += q_i * q_i.transpose();
		b += p.z * q_i;
	}

	// According to this: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
	// This seems like the best decomposition method (fast and accurate)
	s = A.colPivHouseholderQr().solve(b); 

	// Defining Hessian Hw
	Eigen::Matrix2d Hw;
	Hw << 2*s(0), s(1),
		  s(1), 2*s(2);

	// Obtaining the result (eigenvalues)
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(Hw);
	eigensolver.compute(Hw);
	Eigen::Vector2d res = eigensolver.eigenvalues();

	minX = res(0);
	maxX = res(1);

}

// Return the values of the two principal curvatures for this patch

void MongePatch::principalCurvatures(float &kmin, float &kmax) const
{
	kmin = minX;
	kmax = maxX;
}


