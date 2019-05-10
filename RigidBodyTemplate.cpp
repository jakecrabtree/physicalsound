#include "RigidBodyTemplate.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <Eigen/Sparse>
#include "Distance.h"

using namespace std;
using namespace Eigen;

RigidBodyTemplate::RigidBodyTemplate(const std::string &meshFilename, double scale) : volume_(0)
{
    inertiaTensor_.setZero();
    Eigen::MatrixXd mV;
    Eigen::MatrixXi mF;
    igl::readOBJ(meshFilename, mV, mF);

    mV *= scale;
	
    igl::copyleft::tetgen::tetrahedralize(mV, mF, "pq10.414a0.21", V, T, F);
    /*V.resize(4, 3);
	V.row(0) = Eigen::Vector3d(0, 0, 0);
	V.row(1) = Eigen::Vector3d(1, 0, 0);
	V.row(2) = Eigen::Vector3d(0, 1, 0);
	V.row(3) = Eigen::Vector3d(0, 0, 1);
	F.resize(4, 3);
	F.row(0) = Eigen::Vector3i(0, 1, 2);
	F.row(1) = Eigen::Vector3i(0, 2, 3);
	F.row(2) = Eigen::Vector3i(0, 3, 1);
	F.row(3) = Eigen::Vector3i(1, 3, 2);
	T.resize(1, 4);
	T.row(0) = Eigen::Vector4i(0, 1, 2, 3);*/
	computeFaces();
    initialize();
}

RigidBodyTemplate::RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX4i &tets) : volume_(0)
{
    V = verts;
    T = tets;
    computeFaces();
    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate()
{    
}

void RigidBodyTemplate::initialize()
{
    computeVolume();
    Vector3d cm = computeCenterOfMass();
    for(int i=0; i<V.rows(); i++)
        V.row(i) -= cm;

    computeInertiaTensor();    
    computeDistances();
	computeTetVols();
	computePointVolumes();
    computeBeta();
}

void RigidBodyTemplate::computeFaces()
{

    struct triple
    {
        triple(int aa, int bb, int cc) : a(aa), b(bb), c(cc)
        {
            if(a < b)
                std::swap(a,b);
            if(a < c)
                std::swap(a,c);
            if(b < c)
                std::swap(b,c);
        }

        int a, b, c;
        bool operator<(const triple &other) const
        {
            if(a < other.a)
                return true;
            else if(a > other.a)
                return false;
            if(b < other.b)
                return true;
            else if(b > other.b)
                return false;
            return c < other.c;
        }
    };

    int ntets = (int)T.rows();
    MatrixX3i allfaces(4*ntets, 3);
    Matrix<int, 4, 3> faceidx;
    faceidx << 0, 1, 3,
        3, 1, 2,
        3, 2, 0,
        0, 2, 1;

    for(int i=0; i<ntets; i++)
    {
        Vector4i tet = T.row(i);
        for(int face=0; face<4; face++)
        {
            for(int k=0; k<3; k++)
                allfaces(4*i+face, k) = tet[faceidx(face,k)];
        }
    }

    map<triple, vector<int> > faces;
    for(int i=0; i<4*ntets; i++)
    {
        triple t(allfaces(i,0), allfaces(i,1), allfaces(i,2));
        faces[t].push_back(i/4);
    }

    int nfaces=0;
    for(map<triple, vector<int> >::iterator it = faces.begin(); it != faces.end(); ++it)
        if(it->second.size() == 1)
            nfaces++;

    F.resize(nfaces,3);
    int idx=0;

    for(int i=0; i<4*ntets; i++)
    {
        triple t(allfaces(i,0), allfaces(i,1), allfaces(i,2));
        if(faces[t].size() == 1)
        {
            F.row(idx) = allfaces.row(i);
            idx++;
        }        
    }    
}

void RigidBodyTemplate::computeVolume()
{
    volume_ = 0;
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
            centroid += pts[j];
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        centroid /= 3.0;
        volume_ += centroid.dot(normal) * area / 3.0;
    } 
}

Vector3d RigidBodyTemplate::computeCenterOfMass()
{
    Vector3d cm(0, 0, 0);
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        Vector3d term(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                    term[j] += pts[k][j] * pts[l][j];
            }
            term[j] *= area*normal[j] / 12.0;
        }

        cm += term;
    }

    return cm / volume_;
}

void RigidBodyTemplate::computeInertiaTensor()
{
    Vector3d quads(0, 0, 0);
    Vector3d mixed(0, 0, 0);
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();


        Vector3d term(0, 0, 0);
        Vector3d mixterm(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                {
                    for (int m = l; m < 3; m++)
                        term[j] += pts[k][j] * pts[l][j] * pts[m][j];
                }
            }
            term[j] *= area*normal[j] / 30.0;
        }
        double mix = 0;
        for (int j = 0; j < 3; j++)
        {
            mix += 6.0*pts[j][0] * pts[j][1] * pts[j][2];
            for (int k = 0; k < 3; k++)
            {
                mix += 2.0*pts[j][k] * pts[j][(k + 1) % 3] * pts[(j + 1) % 3][(k + 2) % 3];
                mix += 2.0*pts[j][k] * pts[j][(k + 1) % 3] * pts[(j + 2) % 3][(k + 2) % 3];
            }
            mix += pts[j][0] * pts[(j + 1) % 3][1] * pts[(j + 2) % 3][2];
            mix += pts[j][2] * pts[(j + 1) % 3][1] * pts[(j + 2) % 3][0];
        }
        for (int j = 0; j < 3; j++)
            mixterm[j] = mix*area*normal[j] / 60.0;

        quads += term;
        mixed += mixterm;
    }

    inertiaTensor_ << quads[1] + quads[2], -mixed[2], -mixed[1],
        -mixed[2], quads[0] + quads[2], -mixed[0],
        -mixed[1], -mixed[0], quads[0] + quads[1];
}

void RigidBodyTemplate::computeDistances()
{
    int nverts = (int)V.rows();
    int nfaces = (int)F.rows();
    distances.resize(nverts);
    for (int i = 0; i < nverts; i++)
    {
        double dist = numeric_limits<double>::infinity();
        for (int j = 0; j < nfaces; j++)
        {
            double dummy;
            if (Distance::vertexPlaneDistanceLessThan(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dist))
            {
                Vector3d distvec = Distance::vertexFaceDistance(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dummy, dummy, dummy);
                dist = min(dist, distvec.norm());
            }
        }
        distances[i] = dist;        
    }
}

void RigidBodyTemplate::computeTetVols()
{
    Tvol.resize(T.rows());
    for (int r = 0; r < Tvol.size(); ++r){
        Vector3d m0 = V.row(T(r, 0));        
        Vector3d m1 = V.row(T(r, 1));
        Vector3d m2 = V.row(T(r, 2));
        Vector3d m3 = V.row(T(r, 3));
        Tvol[r] = std::abs((1.0/6.0) * ((m1 - m0).cross(m2 - m0)).dot(m3-m0));
    }
}

void RigidBodyTemplate::computePointVolumes() {
    Vvol.resize(V.rows());
    Vvol.setZero();
    for (int r = 0; r < Tvol.size(); ++r){
        int p0 = T(r, 0);        
        int p1 = T(r, 1);
        int p2 = T(r, 2);
        int p3 = T(r, 3);
        Vvol[p0] += Tvol[r]/4.0;
        Vvol[p1] += Tvol[r]/4.0;
        Vvol[p2] += Tvol[r]/4.0;
        Vvol[p3] += Tvol[r]/4.0;
    }
}

void RigidBodyTemplate::computeBeta(){
    betas.clear();
    for (int r = 0; r < Tvol.size(); ++r){
        Vector3d m0 = V.row(T(r, 0));        
        Vector3d m1 = V.row(T(r, 1));
        Vector3d m2 = V.row(T(r, 2));
        Vector3d m3 = V.row(T(r, 3));
        Matrix4d beta;
        beta << m0[0], m1[0], m2[0], m3[0],
                m0[1], m1[1], m2[1], m3[1],
                m0[2], m1[2], m2[2], m3[2],
                1, 1, 1, 1;
        beta = beta.inverse();
        betas.push_back(beta);
    }
}

Eigen::Matrix3d RigidBodyTemplate::cuteLilFunction(int tet) const {
	double x1 = V(T(tet, 0), 0);
	double x2 = V(T(tet, 1), 0);
	double x3 = V(T(tet, 2), 0);
	double x4 = V(T(tet, 3), 0);	
	double y1 = V(T(tet, 0), 1);
	double y2 = V(T(tet, 1), 1);
	double y3 = V(T(tet, 2), 1);
	double y4 = V(T(tet, 3), 1);
	double z1 = V(T(tet, 0), 2);
	double z2 = V(T(tet, 1), 2);
	double z3 = V(T(tet, 2), 2);
	double z4 = V(T(tet, 3), 2);
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
double RigidBodyTemplate::distance(Vector3d p, int tet) const
{
    double x4 = V(T(tet, 3), 0);	
    double y4 = V(T(tet, 3), 1);
    double z4 = V(T(tet, 3), 2);
    Eigen::Matrix3d linTrans = cuteLilFunction(tet);

	Vector3d bary = linTrans.inverse() * (p - Eigen::Vector3d(x4, y4, z4));
	double last = 1 - bary[0] - bary[1] - bary[2];
    double ret = bary[0] * distances[T(tet, 0)] + bary[1] * distances[T(tet, 1)] + bary[2] * distances[T(tet, 2)] + last * distances[T(tet, 3)];
    return ret;
}

Vector3d RigidBodyTemplate::Ddistance(int tet) const
{
    Eigen::Matrix3d linTrans = cuteLilFunction(tet);
    double d0 = distances[T(tet, 0)];
    double d1 = distances[T(tet, 1)];
    double d2 = distances[T(tet, 2)];
    double d3 = distances[T(tet, 3)];
    Vector3d result(0, 0, 0);
    result[0] = d0 - d3;
    result[1] = d1 - d3;
    result[2] = d2 - d3;

    Vector3d ret = (result.transpose() * linTrans.inverse());//.transpose(); //todo hmmmmmmm
    return ret;
}
