#include "SimpleDistance.h"


/* Initialize everything to be able to compute the implicit distance of [Hoppe92] 
   at arbitrary points that are close enough to the point cloud.
 */

void SimpleDistance::init(const PointCloud *pointCloud, float samplingRadius)
{
	sampling_radius = samplingRadius;
	points = pointCloud;

	knn = NearestNeighbors();
	knn.setPoints(&points->getPoints());
}


/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the implicit function of [Hoppe92]
   at point P.
 */

bool SimpleDistance::operator()(const glm::vec3 &P, float &value) const
{
	vector<size_t> knn1;
	vector<float> dist; // not used

	// i = index of tangent plane whose center is closest to P
	knn.getKNearestNeighbors(P, 1, knn1, dist);
	size_t i = knn1[0];

	glm::vec3 o_i = points->getPoints()[i];
	glm::vec3 n_i = points->getNormals()[i];

	// Compute z as the projection of P onto Tp(P)
	float f_p = glm::dot(P-o_i, n_i);
	glm::vec3 z = o_i - (f_p)*n_i;

	if (glm::distance(z, P) < sampling_radius) {
		value = f_p;
		return true;
	}

	return false; // The point is too far from the point cloud (value is undefined)
}






