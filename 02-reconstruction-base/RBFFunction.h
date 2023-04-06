#ifndef _RBF_FUNCTION_INCLUDE
#define _RBF_FUNCTION_INCLUDE


#include "ImplicitFunction.h"
#include "PointCloud.h"
#include "NearestNeighbors.h"


class RBFFunction : public ImplicitFunction
{

public:
	void init(const PointCloud *pointCloud, float standardDeviation, float supportRadius);

	bool operator()(const glm::vec3 &P, float &value) const;
	
private:
	const PointCloud* points;
	float standard_deviation, support_radius;

	NearestNeighbors knn;
	vector<glm::vec3>* p_def;
	vector<double> c_def;

};


#endif // _RBF_FUNCTION_INCLUDE
