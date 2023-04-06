#define NOMINMAX

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "RBFFunction.h"


// Function that returns the value of the Gaussian RBF for a given r
float GaussianRBF(float r, float c) {
    return exp( -pow(r, 2) / (2*pow(c, 2)) );
}

/* Initialize everything to be able to compute the implicit distance to the reconstructed
   point cloud at arbitrary points that are close enough to the point cloud. As should be
   obvious by the name of the class, the distance has to be computed using RBFs.
 */

void RBFFunction::init(const PointCloud *pointCloud, float standardDeviation, float supportRadius)
{
    points = pointCloud;
    standard_deviation = standardDeviation;
    support_radius = supportRadius;

    vector<glm::vec3> P = points->getPoints();
    vector<glm::vec3> N = points->getNormals();

    // Compute the value of d as 0.5% of the diagonal of the model's bounding box (yields good results)
    float x_min = 10000, y_min = 10000, z_min = 10000, x_max = -10000, y_max = -10000, z_max = -10000;
    for (int i = 0; i < P.size(); i++) {
        x_min = min(x_min, P[i].x);
        y_min = min(y_min, P[i].y);
        z_min = min(z_min, P[i].z);
        x_max = max(x_max, P[i].x);
        y_max = max(y_max, P[i].y);
        z_max = max(z_max, P[i].z);
    }
    float d = 0.005f * glm::distance(glm::vec3(x_max, y_max, z_max), glm::vec3(x_min, y_min, z_min));

    // Create artificial points and add them to a definitive point cloud
    // Store the values of v too to solve the system later
    p_def = new vector<glm::vec3>();
    vector<float> v;
    for (int i = 0; i < P.size(); i++) {
        p_def->push_back(P[i]);
        p_def->push_back(P[i] + d*N[i]);
        p_def->push_back(P[i] - d*N[i]);

        v.push_back(0);
        v.push_back(d);
        v.push_back(-d);
    }

    knn = NearestNeighbors();
    knn.setPoints(p_def);

    // Generate the A matrix using a Sparse Matrix to reduce size
    Eigen::SparseMatrix<double> A(p_def->size(), p_def->size());
    vector<Eigen::Triplet<float>> A_Vector;

    vector<pair<size_t, float>> knn_radius;
    glm::vec3 p1, p2;

    for (int i = 0; i < p_def->size(); i++) {
        p1 = p_def->at(i);
        knn.getNeighborsInRadius(p1, support_radius, knn_radius);

        for (int j = 0; j < knn_radius.size(); j++) {
            p2 = p_def->at(knn_radius[j].first);
            float val = GaussianRBF(glm::distance(p1, p2), standard_deviation);
            A_Vector.push_back(Eigen::Triplet<float>(i, knn_radius[j].first, val));
        }
    }

    A.setFromTriplets(A_Vector.begin(), A_Vector.end());

    // Apply a regularization by A' = A + Lambda*I
    float lambda = 0.2f; // Arbitrary value for Lambda (looks good on test results)

    Eigen::SparseMatrix<double> LambdaxI(p_def->size(), p_def->size());
    LambdaxI.setIdentity();
    LambdaxI *= lambda;
    A = A + LambdaxI;

    Eigen::VectorXd v_aux(v.size()), c;
    for (int i = 0; i < v.size(); i++) {
        v_aux[i] = v[i];
    }

    // Solving the system
    // According to this: http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
    // BiCGSTAB is the best option for this situation (Iterative solver for a square matrix)
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> systemsolver;
    systemsolver.compute(A);
    c = systemsolver.solve(v_aux);

    for (int i = 0; i < c.size(); i++) {
        c_def.push_back(c[i]);
    }

}


/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the RBF implicit distance at point P.
 */

bool RBFFunction::operator()(const glm::vec3 &P, float &value) const
{
    vector<pair<size_t, float>> knn_radius;
    knn.getNeighborsInRadius(P, support_radius, knn_radius);
    float res = 0;

    if (knn_radius.size() != 0) {
        // Add all the gaussian values for all the points inside the support radius of P
        for (int i = 0; i < knn_radius.size(); i++) {
            res += GaussianRBF(glm::distance(P, p_def->at(knn_radius[i].first)), 
                standard_deviation) * c_def[knn_radius[i].first];
        }
        value = res;
        return true;
    }

	return false; // The point is too far from the point cloud (value is undefined)
}
