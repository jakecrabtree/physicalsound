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
	V = getTemplate().getVerts();
	for(int i = 0; i < V.rows(); i++) {
		V.row(i) = (c + VectorMath::rotationMatrix(theta)* V.row(i).transpose()).transpose();
	}
	Vdot = Eigen::MatrixXd(V.rows(), V.cols());
	Vdot.setZero();
	for(int i = 0; i < Vdot.rows(); i++) {
		Vdot.row(i) = cvel;
	}
    AABB = buildAABB(this);
}

Vector3d RigidBodyInstance::partialX(int tet, int i){
	Vector3d p0 = V.row(getTemplate().getTets()(tet, 0));
	Vector3d p1 = V.row(getTemplate().getTets()(tet, 1));
	Vector3d p2 = V.row(getTemplate().getTets()(tet, 2));
	Vector3d p3 = V.row(getTemplate().getTets()(tet, 3));
	Matrix<double, 3, 4> p;
	p << p0[0], p1[0], p2[0], p3[0],
		 p0[1], p1[1], p2[1], p3[1],
		 p0[2], p1[2], p2[2], p3[2];
	Vector4d kd;
	kd.setZero();
	kd[i] = 1;
	return p * getTemplate().getBetas()[tet] * kd;
}

Vector3d RigidBodyInstance::partialXdot(int tet, int i){
	Vector3d v0 = Vdot.row(getTemplate().getTets()(tet, 0));
	Vector3d v1 = Vdot.row(getTemplate().getTets()(tet, 1));
	Vector3d v2 = Vdot.row(getTemplate().getTets()(tet, 2));
	Vector3d v3 = Vdot.row(getTemplate().getTets()(tet, 3));
	Matrix<double, 3, 4> v;
	v << v0[0], v1[0], v2[0], v3[0],
		 v0[1], v1[1], v2[1], v3[1],
		 v0[2], v1[2], v2[2], v3[2];
	Vector4d kd;
	kd.setZero();
	kd[i] = 1;
	return v * getTemplate().getBetas()[tet] * kd;
}

double RigidBodyInstance::kd(int i, int j){
	return (i == j) ? 1 : 0;
}


void RigidBodyInstance::strainRateTensor(int tet, Eigen::Matrix3d& nu, Eigen::Matrix3d& Xpartials, Eigen::Matrix3d& Xdotpartials){
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			nu(i,j) = Xpartials.col(i).dot(Xdotpartials.col(j));
			nu(i,j) += Xpartials.col(j).dot(Xdotpartials.col(i));
		}
	}
}


void RigidBodyInstance::strainTensor(int tet, Eigen::Matrix3d& epsilon, Eigen::Matrix3d& Xpartials){
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			epsilon(i,j) = Xpartials.col(i).dot(Xpartials.col(j)) - kd(i,j);
		}
	}
}

void RigidBodyInstance::stressTensor(int tet, Eigen::Matrix3d& sigma){
	sigma.setZero();
	Matrix3d Xpartials;
	Matrix3d Xdotpartials;
	Xpartials.col(0) = partialX(tet, 0);
	Xpartials.col(1) = partialX(tet, 1);
	Xpartials.col(2) = partialX(tet, 2);
	Xdotpartials.col(0) = partialXdot(tet, 0);
	Xdotpartials.col(1) = partialXdot(tet, 1);
	Xdotpartials.col(2) = partialXdot(tet, 2);
	Matrix3d nu;
	strainRateTensor(tet, nu, Xpartials, Xdotpartials);
	Matrix3d epsilon;
	strainTensor(tet, epsilon, Xpartials);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			for (int k = 0; k < 3; ++k){
				sigma(i,j) += lambda * epsilon(k,k) * kd(i,j) + 2.0 * mu * epsilon(i,j);
				sigma(i,j) += phi * nu(k,k) * kd(i,j) + 2.0 * psi * nu(i,j);
			}
		}
	}
}

Vector3d RigidBodyInstance::elasticForce(int tet, int i){
	Vector3d force;
	force.setZero();
	Matrix3d sigma;
	stressTensor(tet, sigma);
	Vector4i vertIndices = getTemplate().getTets().row(tet);
	for (int j = 0; j < vertIndices.size(); ++j){
		Vector3d pj = V.row(vertIndices[j]);
		double mag = 0;
		for (int k = 0; k < 3; ++k){
			for (int l = 0; l < 3; ++l){
				mag += getTemplate().getBetas()[tet](j,l) * getTemplate().getBetas()[tet](i,k) * sigma(k,l);
			}
		}
		force += mag * pj;
	}
	force *= (-1.0/2.0) * getTemplate().getTvol()[tet];
	return force;
}


RigidBodyInstance::~RigidBodyInstance()
{    
}
