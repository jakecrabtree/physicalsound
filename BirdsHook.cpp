#include "BirdsHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include "CollisionDetection.h"
#include <Eigen/Geometry>

using namespace Eigen;

BirdsHook::BirdsHook() : PhysicsHook(), sceneFile_("box.scn")
{
    birdTemplate_ = NULL;
    launch_ = false;
}

void BirdsHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Filename", sceneFile_);
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &params_.timeStep, 0, 0, 3);
        ImGui::DragFloat("Newton Tolerance", &params_.NewtonTolerance, 0.01, 1e-16, 1e-1, "%.3e", 10);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputFloat("Gravity G", &params_.gravityG, 0, 0, 3);
        ImGui::Checkbox("Penalty Forces Enabled", &params_.penaltyEnabled);
        ImGui::InputFloat("Penalty Stiffness", &params_.penaltyStiffness, 0, 0, 3);
        ImGui::Checkbox("Impulses Enabled", &params_.impulsesEnabled);
        ImGui::InputFloat("CoR", &params_.CoR, 0, 0, 3);
    }    
}

void BirdsHook::updateRenderGeometry()
{
    int totverts = 0;
    int totfaces = 0;

    // floor

    totverts += 5;
    totfaces += 4;

    for (RigidBodyInstance *rbi : bodies_)
    {
        totverts += rbi->getTemplate().getVerts().rows();
        totfaces += rbi->getTemplate().getFaces().rows();
    }
    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;

    double floory = -1.0;
    // floor
    renderQ.row(0) << 0, floory, 0;
    renderQ.row(1) << 1e6, floory, 1e6;
    renderQ.row(2) << -1e6, floory, 1e6;
    renderQ.row(3) << -1e6, floory, -1e6;
    renderQ.row(4) << 1e6, floory, -1e6;
    voffset += 5;

    renderF.row(0) << 0, 2, 1;
    renderF.row(1) << 0, 3, 2;
    renderF.row(2) << 0, 4, 3;
    renderF.row(3) << 0, 1, 4;
    foffset += 4;

    for (RigidBodyInstance *rbi : bodies_)
    {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++) {
            //renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
			renderQ.row(voffset + i) = (rbi->V.row(i));
		}
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }
}


void BirdsHook::initSimulation()
{
    time_ = 0;    
    loadScene();
    updateRenderGeometry();
}

void BirdsHook::tick()
{    
    launchMutex_.lock();

    double launchVel = 100.0;

    if (launch_)
    {
        Eigen::Vector3d cvel(0, 0, 0);
        Eigen::Vector3d w(0, 0, 0);
        RigidBodyInstance *rbi = new RigidBodyInstance(*birdTemplate_, launchPos_, cvel, launchVel * launchDir_, w, 1.0);
        bodies_.push_back(rbi);
    
        launch_ = false;
    }

    launchMutex_.unlock();
}

void BirdsHook::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    int numVerts = 0;
    for(int i=0; i<bodies_.size(); i++)
    {
        numVerts += bodies_[i]->getTemplate().getVerts().rows();
    }
    Fc.resize(3*numVerts);
   // Ftheta.resize(3*bodies_.size());
    Fc.setZero();
   // Ftheta.setZero();    
    if(params_.gravityEnabled)
    {
        int count = 0;
        for(int i=0; i<bodies_.size(); i++)
        {
            int numRows = bodies_[i]->getTemplate().getVerts().rows();
            for (int p = 0; p < numRows; ++p){
                double m = bodies_[i]->density * bodies_[i]->getTemplate().getVvol()[p];
                Fc[3*count++ + 1] -= params_.gravityG*m;
            }
        }
    }
    int count = 0;
    for(int i = 0; i < bodies_.size(); i++) {
        int numRows = bodies_[i]->getTemplate().getVerts().rows();
        int numTets = bodies_[i]->getTemplate().getTets().rows();
        for(int t = 0; t < numTets; t++) {
            for(int j = 0; j < 4; j++) {
                Vector3d toAdd = bodies_[i]->elasticForce(t,j);
                int vidx = bodies_[i]->getTemplate().getTets()(t,j);
                Fc[count + 3 * vidx] += toAdd[0];
                Fc[count + 3 * vidx + 1] += toAdd[1];
                Fc[count + 3 * vidx + 2] += toAdd[2];
            }
        }
        count += 3 * numRows;
    }
}

bool BirdsHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir, int button)
{
    if (button != 2)
        return false;

    launchMutex_.lock();
    launch_ = true;
    Eigen::Matrix4f view = viewer.core.view;
    Eigen::Vector4f eye = view.inverse() * Eigen::Vector4f(0, 0, 0, 1.0f);
    for (int i = 0; i < 3; i++)
    {

        launchPos_[i] = eye[i] + dir[i];
        launchDir_[i] = dir[i];
    }

    launchMutex_.unlock();
    return true;
}

bool BirdsHook::simulateOneStep()
{
	time_ += params_.timeStep;
	int nbodies = (int)bodies_.size();
	for(int bodyidx = 0; bodyidx < nbodies; bodyidx++) {
		RigidBodyInstance &body = *bodies_[bodyidx];
		body.V += params_.timeStep * body.Vdot;
	}
    Eigen::VectorXd cForce;
    Eigen::VectorXd thetaForce;
    computeForces(cForce, thetaForce);
	
	//Collision stuff
	std::set<Collision> collisions;
	collisionDetection(bodies_, collisions);
	std::vector<int> pre;
	pre.push_back(0);
	for(int bodyidx = 0; bodyidx < nbodies; bodyidx++) {
		pre.push_back(pre[bodyidx] + bodies_[bodyidx]->getTemplate().getVerts().rows());
	}
	for(const Collision& collision: collisions) {
		if(collision.body2 == -1) {
			int vidx = collision.collidingVertex;
			double vy = bodies_[collision.body1]->V(vidx, 1);
			double dist = vy + 1;
			cForce[3 * pre[collision.body1] + 3 * vidx + 1] += -1 * params_.penaltyStiffness * dist;
		} 
	}

    int counter = 0;
    for(int bodyidx = 0; bodyidx < nbodies; bodyidx++) {
		RigidBodyInstance &body = *bodies_[bodyidx];
        int numRows = body.getTemplate().getVerts().rows();
        for (int p = 0; p < numRows; ++p){
            body.Vdot.row(p) += params_.timeStep*cForce.segment<3>(3 * counter++) / body.density / body.getTemplate().getVvol()[p];
        }
	}    
	return false;
}

/*bool BirdsHook::simulateOneStep()
{   
    time_ += params_.timeStep;
    int nbodies = (int)bodies_.size();

    std::vector<Vector3d> oldthetas;
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.c += params_.timeStep*body.cvel;
        Matrix3d Rhw = VectorMath::rotationMatrix(params_.timeStep*body.w);
        Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);

        Vector3d oldtheta = body.theta;
        body.theta = VectorMath::axisAngle(Rtheta*Rhw);  
        if (body.theta.dot(oldtheta) < 0 && oldtheta.norm() > M_PI/2.0)
        {
            double oldnorm = oldtheta.norm();
            oldtheta = (oldnorm - 2.0*M_PI)*oldtheta/oldnorm;
        }
        oldthetas.push_back(oldtheta);
    }

    std::set<Collision> collisions;
    collisionDetection(bodies_, collisions);

    Eigen::VectorXd cForce(3 * nbodies);
    Eigen::VectorXd thetaForce(3 * nbodies);
    computeForces(cForce, thetaForce);
    
    // TODO compute and add penalty forces
    // TODO apply collision impulses
    
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {        
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d Mi = body.getTemplate().getInertiaTensor();

        body.cvel += params_.timeStep*cForce.segment<3>(3 * bodyidx) / body.density / body.getTemplate().getVolume();

        Vector3d newwguess(body.w);

        int iter = 0;
        for (iter = 0; iter < params_.NewtonMaxIters; iter++)
        {
            Vector3d term1 = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * newwguess;
            Vector3d term2 = (VectorMath::TMatrix(params_.timeStep*body.w).inverse()*VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * body.w;
            Vector3d term3 = -params_.timeStep * thetaForce.segment<3>(3 * bodyidx);
            Vector3d fval = term1 + term2 + term3;
            if (fval.norm() / body.density / Mi.trace() <= params_.NewtonTolerance)
                break;

            Matrix3d Df = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

            Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
        body.w = newwguess;

    }

    return false;
}*/

void BirdsHook::loadScene()
{
    for (RigidBodyInstance *rbi : bodies_)
        delete rbi;
    for (RigidBodyTemplate *rbt : templates_)
        delete rbt;
    bodies_.clear();
    templates_.clear();

    std::string prefix;
    std::string scenefname = std::string("scenes/") + sceneFile_;
    std::ifstream ifs(scenefname);
    if (!ifs)
    {
        // run from the build directory?
        prefix = "../";
        scenefname = prefix + scenefname;        
        ifs.open(scenefname);
        if(!ifs)
            return;
    }
        

    int nbodies;
    ifs >> nbodies;
    for (int body = 0; body < nbodies; body++)
    {
        std::string meshname;
        ifs >> meshname;
        meshname = prefix + std::string("meshes/") + meshname;
        double scale;
        ifs >> scale;
        RigidBodyTemplate *rbt = new RigidBodyTemplate(meshname, scale);
        double rho;
        ifs >> rho;
        Eigen::Vector3d c, theta, cvel, w;
        for (int i = 0; i < 3; i++)
            ifs >> c[i];
        for (int i = 0; i < 3; i++)
            ifs >> theta[i];
        for (int i = 0; i < 3; i++)
            ifs >> cvel[i];
        for (int i = 0; i < 3; i++)
            ifs >> w[i];
        RigidBodyInstance *rbi = new RigidBodyInstance(*rbt, c, theta, cvel, w, rho);
        templates_.push_back(rbt);
        bodies_.push_back(rbi);
    }

    // bird mesh    
    std::string birdname = prefix + std::string("meshes/bird2.obj");
    delete birdTemplate_;
    birdTemplate_ = new RigidBodyTemplate(birdname, 0.1);
}
