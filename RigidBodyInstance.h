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
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta, const Eigen::Vector3d &cvel, const Eigen::Vector3d &w, double density, double young, double mu, double phi, double psi);
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
    double lambda;
    double mu;
    double phi;
    double psi;

    double acousticImpedance = 415;

    Eigen::Vector3d elasticForce(int tet, int i);

    double distance(Eigen::Vector3d p, int tet) const;
    Eigen::Vector3d Ddistance(int tet) const;
	Eigen::Matrix3d cuteLilFunction(int tet) const;	
    void computeFacePressures(Eigen::VectorXd& pressures);

    std::vector<Eigen::VectorXd> facePressures;
    std::vector<double> facePressureDelays;
    int currFacePressure = 0;

    
private:
    const RigidBodyTemplate &rbtemplate_;
    Eigen::Vector3d partialX(int tet, int i);
    void partialXMatrix(int tet, Eigen::Matrix3d& Xpartials);
    Eigen::Vector3d partialXdot(int tet, int i);
    void partialXdotMatrix(int tet, Eigen::Matrix3d& Xdotpartials);
    void strainTensor(int tet, Eigen::Matrix3d& epsilon, Eigen::Matrix3d& Xpartials);
    void strainRateTensor(int tet, Eigen::Matrix3d& nu, Eigen::Matrix3d& Xpartials, Eigen::Matrix3d& Xdotpartials);
    void stressTensor(int tet, Eigen::Matrix3d& sigma);
    double kd(int i, int j);
};

#endif // RIGIDBODYINSTANCE_H
