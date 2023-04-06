#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "LaplacianSmoothing.h"


void LaplacianSmoothing::setMesh(TriangleMesh *newMesh)
{
	mesh = newMesh;
}

// Applies the Laplacian to given set of points using a different multiplier (used for the different operators)
vector<glm::vec3> LaplacianSmoothing::applyLaplacian(vector<glm::vec3> P, float multiplier) {
    vector<glm::vec3> P_prime;
    vector<unsigned int> neighbors;
    glm::vec3 pi, pj;
    glm::vec3 laplacian = glm::vec3(0, 0, 0);

    for (int i = 0; i < P.size(); i++) {
        laplacian = glm::vec3(0, 0, 0);
        pi = P[i];
        mesh->getNeighbors(i, neighbors);
        
        for (int j = 0; j < neighbors.size(); j++) {
            pj = P[neighbors[j]];
            laplacian += pj - pi;
        }

        laplacian /= neighbors.size();
        P_prime.push_back(pi + multiplier * laplacian);
    }

    return P_prime;
}

/* This method should apply nIterations iterations of the laplacian vector multiplied by lambda 
   to each of the vertices. */

void LaplacianSmoothing::iterativeLaplacian(int nIterations, float lambda)
{
    vector<glm::vec3> P_prime;

    for (int i = 0; i < nIterations; i++) {
        vector<glm::vec3> P = mesh->getVertices(); // Need to get updated points every iteration
        P_prime = applyLaplacian(P, lambda);

        // Apply the changes only after all the points have been modified
        for (int j = 0; j < P.size(); j++) {
            mesh->getVertices()[j] = P_prime[j];
        }
    }
}

/* This method should apply nIterations iterations of the bilaplacian operator using lambda 
   as a scaling factor. */

void LaplacianSmoothing::iterativeBilaplacian(int nIterations, float lambda)
{
    vector<glm::vec3> P_prime, P_prime2;

    for (int i = 0; i < nIterations; i++) {
        vector<glm::vec3> P = mesh->getVertices(); // Need to get updated points every iteration
        P_prime = applyLaplacian(P, lambda);
        P_prime2 = applyLaplacian(P_prime, -lambda); // Once the first transformation is done, do the second one

        // Apply the changes only after all the points have been modified
        for (int j = 0; j < P.size(); j++) {
            mesh->getVertices()[j] = P_prime2[j];
        }
    }
}

/* This method should apply nIterations iterations of Taubin's operator using lambda 
   as a scaling factor, and computing the corresponding nu value. */

void LaplacianSmoothing::iterativeLambdaNu(int nIterations, float lambda)
{
    vector<glm::vec3> P_prime, P_prime2;

    float nu = lambda / (0.1f * lambda - 1); // Using the formula from the notes

    for (int i = 0; i < nIterations; i++) {
        vector<glm::vec3> P = mesh->getVertices(); // Need to get updated points every iteration
        P_prime = applyLaplacian(P, lambda);
        P_prime2 = applyLaplacian(P_prime, nu); // Once the first transformation is done, do the second one

        // Apply the changes only after all the points have been modified
        for (int j = 0; j < P.size(); j++) {
            mesh->getVertices()[j] = P_prime2[j];
        }
    }
}

/* This method should compute new vertices positions by making the laplacian zero, while
   maintaing the vertices marked as constraints fixed. */

void LaplacianSmoothing::globalLaplacian(const vector<bool> &constraints)
{
    
}

/* This method has to optimize the vertices' positions in the least squares sense, 
   so that the laplacian is close to zero and the vertices remain close to their 
   original locations. The constraintWeight parameter is used to control how close 
   the vertices have to be to their original positions. */

void LaplacianSmoothing::globalBilaplacian(const vector<bool> &constraints, float constraintWeight)
{

}








