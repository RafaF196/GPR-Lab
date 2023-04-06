#include "NormalEstimator.h"
#include "NearestNeighbors.h"
#include <Eigen/Dense>

// This method has to compute a normal per point in the 'points' vector and put it in the 
// 'normals' vector. The 'normals' vector already has the same size as the 'points' vector. 
// There is no need to push_back new elements, only to overwrite them ('normals[i] = ...;')
// In order to compute the normals you should use PCA. The provided 'NearestNeighbors' class
// wraps the nanoflann library that computes K-nearest neighbors effciently. 


void NormalEstimator::computePointCloudNormals(const vector<glm::vec3> &points, vector<glm::vec3> &normals)
{
	int knn_value = 10;
	vector <size_t> V;
	vector <float> dist; // not used
	vector<glm::vec3> new_p = vector<glm::vec3>(knn_value);

	Eigen::Matrix3f tempM;

	// Initialize a NearestNeighbors instance
	NearestNeighbors knn = NearestNeighbors();
	knn.setPoints(&points);
	
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
		// Since the eigenvalues (and corresponding eigenvectors) are already sorted
		// by the solver, we can take the first eigenvector straightaway

		// Obtaining the result (the eigenvector is already normalized, no need to do it again)
		if (tempM(2, 0) < 0) { // If z component of the normal is negative, invert the normal
			normals[p].x = -1*tempM(0, 0);
			normals[p].y = -1*tempM(1, 0);
			normals[p].z = -1*tempM(2, 0);
		} else {
			normals[p].x = tempM(0, 0);
			normals[p].y = tempM(1, 0);
			normals[p].z = tempM(2, 0);
		}

	}

}

