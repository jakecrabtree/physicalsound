#include "PhysicsHook.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerData.h>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include "SimParameters.h"
#include <set>
#include "CollisionDetection.h"
#include <AudioPlayer.h>
#include <iostream>
#include <fstream>

class RigidBodyTemplate;
class RigidBodyInstance;

class BirdsHook : public PhysicsHook
{
public:
    BirdsHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    virtual void initSimulation(int mode);

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir, int button);

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        viewer.data().set_mesh(renderQ, renderF);
		viewer.data().set_colors(renderC);
    }

    void setCamera(Eigen::Vector3f* cameraPos){
        cam = cameraPos;
    }

private:
    void loadScene();
    void computeForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta);    

    std::mutex launchMutex_;
    bool launch_;
    Eigen::Vector3d launchPos_;
    Eigen::Vector3d launchDir_;
	
	int lastTime = -1;
    double time_;
    SimParameters params_;
    std::string sceneFile_;	

    std::vector<RigidBodyTemplate *> templates_;
    std::vector<RigidBodyInstance *> bodies_;

    RigidBodyTemplate *birdTemplate_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
	Eigen::MatrixXd renderC;
	
	std::vector<Eigen::MatrixXd> playbackData;
	int mode = 0;	
	std::ofstream ofs;
	AudioPlayer aud;
    Eigen::Vector3f* cam;
};
