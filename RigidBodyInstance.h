#ifndef RIGIDBODYINSTANCE_H
#define RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>
#include <vector>

class RigidBodyTemplate;
struct AABBNode;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta, const Eigen::Vector3d &cvel, const Eigen::Vector3d &w, double density);
    ~RigidBodyInstance();

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;	

    Eigen::MatrixX3d V;
	Eigen::MatrixX3d Vdot;

    double density;

    AABBNode *AABB;
    
    const RigidBodyTemplate &getTemplate() const {return rbtemplate_;}
    double lambda = 1;
    double mu = 1;
    double phi = 1;
    double psi = 1;

    Eigen::Vector3d elasticForce(int tet, int i);

    
private:
    const RigidBodyTemplate &rbtemplate_;
    Eigen::Vector3d partialX(int tet, int i);
    Eigen::Vector3d partialXdot(int tet, int i);
    Eigen::Matrix3d strainTensor(int tet);
    Eigen::Matrix3d strainRateTensor(int tet);
    Eigen::Matrix3d stressTensor(int tet);
    double kd(int i, int j);
};

#endif // RIGIDBODYINSTANCE_H
