#include <iostream>
#include <algorithm>
#include <set>
#include "IterativeClosestPoint.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <glm/gtc/matrix_transform.hpp>


#define PI 3.14159265


void IterativeClosestPoint::setClouds(PointCloud *pointCloud1, PointCloud *pointCloud2)
{
	cloud1 = pointCloud1;
	cloud2 = pointCloud2;
	knn.setPoints(&(cloud1->getPoints()));
}

// This method should mark the border points in cloud 1. It also changes their color (for example to red).
// You will need to add an attribute to this class that stores this property for all points in cloud 1. 

void IterativeClosestPoint::markBorderPoints()
{
	int knn_value = 10; // Arbitrary value, looks good enough
	vector<size_t> V;
	vector<float> dist; // not used
	vector<glm::vec3> new_p = vector<glm::vec3>(knn_value);

	Eigen::Matrix3f tempM;
	Eigen::Vector3f v1, v2, v3, center, eigenpoint, aux;

	vector<Eigen::Vector3f> transformed_knn = vector<Eigen::Vector3f>(knn_value);
	vector<float> angles = vector<float>(knn_value);

	float angle_threshold = 1.22f; // Value fixed after some testing on the bunny scans

	vector<glm::vec3> points = cloud1->getPoints();

	for (int p = 0; p < points.size(); p++) {

		// Obtain the k nearest neighbors to the point
		knn.getKNearestNeighbors(points[p], knn_value, V, dist);

		// Compute the centroid
		glm::vec3 centroid = glm::vec3(0.0f);
		for (int i = 0; i < knn_value; i++) {
			centroid += points[V[i]];
		}
		centroid /= (float)knn_value;

		// Translate the points
		for (int i = 0; i < knn_value; i++) {
			new_p[i] = points[V[i]] - centroid;
		}

		// Compute the C matrix
		Eigen::Matrix3f C, Ctemp;
		C << 0, 0, 0, 0, 0, 0, 0, 0, 0;
		for (int i = 0; i < knn_value; i++) {
			Ctemp << new_p[i].x * new_p[i].x, new_p[i].x* new_p[i].y, new_p[i].x* new_p[i].z,
				new_p[i].y* new_p[i].x, new_p[i].y* new_p[i].y, new_p[i].y* new_p[i].z,
				new_p[i].z* new_p[i].x, new_p[i].z* new_p[i].y, new_p[i].z* new_p[i].z;
			C += Ctemp;
		}

		// Obtaining the eigenvalues and eigenvectors
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(C);
		tempM << eigensolver.eigenvectors();

		v3 << tempM(0, 0), tempM(1, 0), tempM(2, 0);
		v2 << tempM(0, 1), tempM(1, 1), tempM(2, 1);
		v1 << tempM(0, 2), tempM(1, 2), tempM(2, 2);
		center << centroid.x, centroid.y, centroid.z;

		// Transforming the k-nearest neighbors to the new frame.
		for (int i = 0; i < knn_value; i++) {
			eigenpoint << points[V[i]].x, points[V[i]].y, points[V[i]].z;
			aux << (eigenpoint - center).dot(v1), (eigenpoint - center).dot(v2), 0; // Projected onto the z-plane
			transformed_knn[i] = aux;
		}

		// Computing the angles for each point and sorting them
		for (int i = 0; i < knn_value; i++) {
			angles[i] = atan2(transformed_knn[i].y(), transformed_knn[i].x());
		}
		sort(angles.begin(), angles.end());

		// Calculating the maximum difference
		float maxDeltaAlpha = 0.0f;
		for (int i = 0; i < knn_value - 1; i++) {
			if (angles[i + 1] - angles[i] > maxDeltaAlpha) maxDeltaAlpha = angles[i + 1] - angles[i];
		}
		if ((angles[0]+2*PI) - angles[knn_value-1] > maxDeltaAlpha) 
			maxDeltaAlpha = (angles[0] + 2 * PI) - angles[knn_value - 1];

		if (maxDeltaAlpha > angle_threshold) border_vertices.push_back(true);
		else border_vertices.push_back(false);

	}

}

// Returns the border vertices. Used to color them differently.

vector<bool> IterativeClosestPoint::getBorderVertices() {
	return border_vertices;
}


// This method should compute the closest point in cloud 1 for all non border points in cloud 2. 
// This correspondence will be useful to compute the ICP step matrix that will get cloud 2 closer to cloud 1.
// Store the correspondence in this class as the following method is going to need it.
// As it is evident in its signature this method also returns the correspondence. The application draws this if available.

vector<int> *IterativeClosestPoint::computeCorrespondence()
{
	correspondence = new vector<int>();
	set<size_t> cor_set;

	vector<glm::vec3> p2 = cloud2->getPoints();

	for (int i = 0; i < border_vertices.size(); i++) {
		if (border_vertices[i]) cor_set.insert(i);
	}

	vector<size_t> knn1;
	vector<float> dist; // not used

	for (int i = 0; i < p2.size(); i++) {
		knn.getKNearestNeighbors(p2[i], 1, knn1, dist);
		if (cor_set.find(knn1[0]) != cor_set.end()) correspondence->push_back(-1);
		else correspondence->push_back(knn1[0]);
	}

	return correspondence;

}


// This method should compute the rotation and translation of an ICP step from the correspondence
// information between clouds 1 and 2. Both should be encoded in the returned 4x4 matrix.
// To do this use the SVD algorithm in Eigen.

glm::mat4 IterativeClosestPoint::computeICPStep()
{
	vector<glm::vec3> p1 = cloud1->getPoints(), p2 = cloud2->getPoints();
	vector<Eigen::Vector3f> pi, qi, adjPi, adjQi;

	Eigen::Vector3f aux1, aux2;
	int index;

	// Obtaining p_i and q_i point sets (correspondence dependent)
	for (int i = 0; i < p2.size(); i++) {
		index = correspondence->at(i);
		if (index > -1) {
			aux1 = Eigen::Vector3f(p2[i].x, p2[i].y, p2[i].z);
			aux2 = Eigen::Vector3f(p1[correspondence->at(i)].x, p1[correspondence->at(i)].y,
				p1[correspondence->at(i)].z);
			qi.push_back(aux1);
			pi.push_back(aux2);
		}
	}

	// Adjust p_i and q_i to their respective centroid
	Eigen::Vector3f centP(0.0f, 0.0f, 0.0f), centQ(0.0f, 0.0f, 0.0f);;
	for (int i = 0; i < pi.size(); i++) {
		centP += pi[i];
	}
	centP /= (float)pi.size();
	for (int i = 0; i < pi.size(); i++) {
		adjPi.push_back(pi[i] - centP);
	}

	for (int i = 0; i < qi.size(); i++) {
		centQ += qi[i];
	}
	centQ /= (float)qi.size();
	for (int i = 0; i < qi.size(); i++) {
		adjQi.push_back(qi[i] - centQ);
	}

	// Generate matrices P and Q
	Eigen::Matrix3Xf P(3, adjPi.size()), Q(3, adjQi.size());
	for (int i = 0; i < pi.size(); i++) {
		P.col(i) = pi[i];
		Q.col(i) = qi[i];
	}

	// Computing the covariance matrix S of both point sets
	Eigen::Matrix3f S = Q * P.transpose();

	// Performing the singular value decomposition of S
	Eigen::JacobiSVD<Eigen::Matrix3Xf> SVD(S, Eigen::ComputeFullU | Eigen::ComputeFullV);

	// Computing the optimal rotation matrix
	Eigen::Matrix3f R = SVD.matrixV() * SVD.matrixU().transpose();

	// No reflection correction needed

	Eigen::Vector3f translation = centP - R * centQ; // translation vector

	// Generating the 4x4 transformation matrix
	Eigen::Matrix4f Rtrans;
	Rtrans << R(0, 0), R(0, 1), R(0, 2), translation.x(),
		R(1, 0), R(1, 1), R(1, 2), translation.y(),
		R(2, 0), R(2, 1), R(2, 2), translation.z(),
		0, 0, 0, 1;

	frobenius_norm = Rtrans.norm(); // Frobenius norm calculation

	glm::mat4 res;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			res[i][j] = Rtrans(j, i); // transpose the transformation matrix before returning it
		}
	}

	return res;

}


// This method should perform the whole ICP algorithm with as many steps as needed.
// It should stop when maxSteps are performed, when the Frobenius norm of the transformation matrix of
// a step is smaller than a small threshold, or when the correspondence does not change from the 
// previous step.

vector<int> *IterativeClosestPoint::computeFullICP(unsigned int maxSteps)
{
	vector<int>* previous_cor = computeCorrespondence(); // Initial correspondance
	vector<int>* actual_cor;

	for (int i = 0; i < maxSteps; i++) {
		cloud2->transform(computeICPStep()); // Transformation
		actual_cor = computeCorrespondence(); // New correspondance

		// First stopping criteria: When the norm of t and the Frobenius norm of R − I are smaller than a threshold.
		if (frobenius_norm < 0.0001f) i = maxSteps;

		// Second stopping criteria: When the correspondence computed is the same as it was in the previous iteration
		bool is_correspondance_equal = true;
		for (int i = 0; i < actual_cor->size(); i++) {
			if (actual_cor->at(i) != previous_cor->at(i)) {
				is_correspondance_equal = false;
				break;
			}
		}
		if (is_correspondance_equal) i = maxSteps;

		previous_cor = actual_cor;
	}

	return previous_cor;
}
