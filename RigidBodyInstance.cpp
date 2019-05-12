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
    double density, double young, double poisson, double phi, double psi)
    : c(c), theta(theta), cvel(cvel), w(w), density(density), rbtemplate_(rbtemplate), phi(phi), psi(psi)
{
	
	double v = (double) poisson;
	double E = (double) young;
	lambda = /*2.65e6;//*/v * E / ((1 + v) * (1 - 2 * v));
	mu = /*3.97e6;//*/E / (2 * (1 + v));
	//phi = 264;
	//psi = 367;
	//density = 1013;
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
	p.col(0) = p0;
	p.col(1) = p1;
	p.col(2) = p2;
	p.col(3) = p3;
	/*p << p0[0], p1[0], p2[0], p3[0],
		 p0[1], p1[1], p2[1], p3[1],
		 p0[2], p1[2], p2[2], p3[2];*/
	Vector4d kd;
	kd.setZero();
	kd[i] = 1;
	return p * getTemplate().getBetas()[tet] * kd;
}

void RigidBodyInstance::partialXMatrix(int tet, Eigen::Matrix3d& Xpartials) {
	Vector3d p0 = V.row(getTemplate().getTets()(tet, 0));
	Vector3d p1 = V.row(getTemplate().getTets()(tet, 1));
	Vector3d p2 = V.row(getTemplate().getTets()(tet, 2));
	Vector3d p3 = V.row(getTemplate().getTets()(tet, 3));
	Matrix<double, 3, 4> p;
	p.col(0) = p0;
	p.col(1) = p1;
	p.col(2) = p2;
	p.col(3) = p3;
	/*p << p0[0], p1[0], p2[0], p3[0],
		 p0[1], p1[1], p2[1], p3[1],
		 p0[2], p1[2], p2[2], p3[2];*/
	Xpartials = (p * getTemplate().getBetas()[tet]).block<3, 3>(0, 0);
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

void RigidBodyInstance::partialXdotMatrix(int tet, Eigen::Matrix3d& Xdotpartials) {
	Vector3d v0 = Vdot.row(getTemplate().getTets()(tet, 0));
	Vector3d v1 = Vdot.row(getTemplate().getTets()(tet, 1));
	Vector3d v2 = Vdot.row(getTemplate().getTets()(tet, 2));
	Vector3d v3 = Vdot.row(getTemplate().getTets()(tet, 3));
	Matrix<double, 3, 4> v;
	v << v0[0], v1[0], v2[0], v3[0],
		 v0[1], v1[1], v2[1], v3[1],
		 v0[2], v1[2], v2[2], v3[2];
	Xdotpartials = (v * getTemplate().getBetas()[tet]).block<3, 3>(0, 0);
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
    partialXMatrix(tet, Xpartials);
    partialXdotMatrix(tet, Xdotpartials);
	//Xpartials.col(0) = partialX(tet, 0);
	//Xpartials.col(1) = partialX(tet, 1);
	//Xpartials.col(2) = partialX(tet, 2);
	//Xdotpartials.col(0) = partialXdot(tet, 0);
	//Xdotpartials.col(1) = partialXdot(tet, 1);
	//Xdotpartials.col(2) = partialXdot(tet, 2);
	Matrix3d nu;
	strainRateTensor(tet, nu, Xpartials, Xdotpartials);
	Matrix3d epsilon;
	strainTensor(tet, epsilon, Xpartials);
    double etrace = epsilon.trace();
    double ntrace = nu.trace();
    Matrix3d I;
    I.setIdentity();
	/*for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			for (int k = 0; k < 3; ++k){
				sigma(i,j) += lambda * epsilon(k,k) * kd(i,j) + 2.0 * mu * epsilon(i,j);
				sigma(i,j) += phi * nu(k,k) * kd(i,j) + 2.0 * psi * nu(i,j);
			}
		}
	}*/
    sigma += lambda * etrace * I + 6.0 * mu * epsilon;
    sigma += phi * ntrace * I + 6.0 * psi * nu;
}

Vector3d RigidBodyInstance::elasticForce(int tet, int i){
	Vector3d force;
	force.setZero();
	Matrix3d sigma;
	stressTensor(tet, sigma);
	Vector4i vertIndices = getTemplate().getTets().row(tet);
	for (int j = 0; j < 4; ++j){
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

Eigen::Matrix3d RigidBodyInstance::cuteLilFunction(int tet) const {
	double x1 = V(getTemplate().getTets()(tet, 0), 0);
	double x2 = V(getTemplate().getTets()(tet, 1), 0);
	double x3 = V(getTemplate().getTets()(tet, 2), 0);
	double x4 = V(getTemplate().getTets()(tet, 3), 0);	
	double y1 = V(getTemplate().getTets()(tet, 0), 1);
	double y2 = V(getTemplate().getTets()(tet, 1), 1);
	double y3 = V(getTemplate().getTets()(tet, 2), 1);
	double y4 = V(getTemplate().getTets()(tet, 3), 1);
	double z1 = V(getTemplate().getTets()(tet, 0), 2);
	double z2 = V(getTemplate().getTets()(tet, 1), 2);
	double z3 = V(getTemplate().getTets()(tet, 2), 2);
	double z4 = V(getTemplate().getTets()(tet, 3), 2);
	Eigen::Matrix3d linTrans;
	linTrans(0, 0) = x1 - x4;
	linTrans(0, 1) = x2 - x4;
	linTrans(0, 2) = x3 - x4; 
	linTrans(1, 0) = y1 - y4;
	linTrans(1, 1) = y2 - y4;
	linTrans(1, 2) = y3 - y4; 
	linTrans(2, 0) = z1 - z4;
	linTrans(2, 1) = z2 - z4;
	linTrans(2, 2) = z3 - z4;
    return linTrans;
}

//p must be in template coordinates
double RigidBodyInstance::distance(Vector3d p, int tet) const
{
    double x4 = V(getTemplate().getTets()(tet, 3), 0);	
    double y4 = V(getTemplate().getTets()(tet, 3), 1);
    double z4 = V(getTemplate().getTets()(tet, 3), 2);
    Eigen::Matrix3d linTrans = cuteLilFunction(tet);

	Vector3d bary = linTrans.inverse() * (p - Eigen::Vector3d(x4, y4, z4));
	double last = 1 - bary[0] - bary[1] - bary[2];
    double ret = bary[0] * getTemplate().getDistance(getTemplate().getTets()(tet, 0)) + bary[1] * getTemplate().getDistance(getTemplate().getTets()(tet, 1)) + bary[2] * getTemplate().getDistance(getTemplate().getTets()(tet, 2)) + last * getTemplate().getDistance(getTemplate().getTets()(tet, 3));
    return ret;
}

Vector3d RigidBodyInstance::Ddistance(int tet) const
{
    Eigen::Matrix3d linTrans = cuteLilFunction(tet);
    double d0 = getTemplate().getDistance(getTemplate().getTets()(tet, 0));
    double d1 = getTemplate().getDistance(getTemplate().getTets()(tet, 1));
    double d2 = getTemplate().getDistance(getTemplate().getTets()(tet, 2));
    double d3 = getTemplate().getDistance(getTemplate().getTets()(tet, 3));
    Vector3d result(0, 0, 0);
    result[0] = d0 - d3;
    result[1] = d1 - d3;
    result[2] = d2 - d3;

    Vector3d ret = (result.transpose() * linTrans.inverse());//.transpose(); //todo hmmmmmmm
    return ret;
}

void RigidBodyInstance::computeFacePressures(Eigen::VectorXd& pressures){
	pressures.setZero();
	pressures.resize(getTemplate().getFaces().rows());
	for (int i = 0; i < pressures.size(); ++i){
		int i0 = getTemplate().getFaces()(i, 0);
		int i1 = getTemplate().getFaces()(i, 1);
		int i2 = getTemplate().getFaces()(i, 2);
		Vector3d p0 = V.row(i0);
		Vector3d p1 = V.row(i1);
		Vector3d p2 = V.row(i2);
		Vector3d norm = ((p1-p0).cross(p2-p0)).normalized();

		Vector3d vel = (Vdot.row(i0) + Vdot.row(i1) + Vdot.row(i2))/3.0;
		double pressure = acousticImpedance * vel.dot(norm);
		pressures[i] = pressure;
	}	
}



