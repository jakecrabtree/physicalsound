#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters()
    {
        timeStep = 0.0001;
        NewtonMaxIters = 50;
        NewtonTolerance = 1e-8;
        
        gravityEnabled = true;
        gravityG = 9.8;                
        penaltyEnabled = true;
        penaltyStiffness = 100.0;
        impulsesEnabled = true;
        CoR = 0.5;
		young = 165;
		poisson = .33;
		maxFrame = 240;
    }

    float timeStep;
    float NewtonTolerance;
    int NewtonMaxIters;
    
    bool gravityEnabled;
    float gravityG;    
    bool penaltyEnabled;
    float penaltyStiffness;
    bool impulsesEnabled;
    float CoR;
	float young;
	float poisson;
	int maxFrame;
};

#endif
