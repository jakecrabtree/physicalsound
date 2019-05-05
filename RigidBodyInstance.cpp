#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "RigidBodyTemplate.h"
#include <Eigen/Geometry>
#include <iostream>
#include "CollisionDetection.h"

using namespace Eigen;
using namespace std;

RigidBodyInstance::RigidBodyInstance(const RigidBodyTemplate &rbtemplate,
    const Eigen::Vector3d &c, const Eigen::Vector3d &theta,
    const Eigen::Vector3d &cvel, const Eigen::Vector3d &w,
    double density)
    : c(c), theta(theta), cvel(cvel), w(w), density(density), rbtemplate_(rbtemplate)
{
    AABB = buildAABB(this);
	V = getTemplate().getVerts();
	for(int i = 0; i < V.rows(); i++) {
		V.row(i) = (c + VectorMath::rotationMatrix(theta)* V.row(i).transpose()).transpose();
	}
	Vdot = Eigen::MatrixXd(V.rows(), V.cols());
	Vdot.setZero();
	for(int i = 0; i  < Vdot.rows(); i++) {
		Vdot.row(i) = Vector3d(0, i % 3 == 0 ? 0.2 : -0.2 , 0);
	}
	Vmass = VectorXd(V.rows());
	Vmass.setZero();
}

RigidBodyInstance::~RigidBodyInstance()
{    
}
