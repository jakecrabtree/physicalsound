Project By: Jake Crouch and Jake Crabtree

How to Build/Run:
Go into project directory
Create build directory
Go into build directory
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./birds_bin

How to use:

There are two basic modes of operation, simulating the sound/visuals and playing them back.

To create a new simulation:
Type in a scene name (no file suffix) and hit reset simulation
Set max frames to be the number of frames you want to simulate for
Set Penalty Stiffness to be the value you want
Set timestep to reasonable value (at least < 10^-6 for good sounds, might be able to get away with < 10^-5 if lucky)
Run simulation
When Max Frames have been reached the program will exit and have written two files to renders (a .fac with face data and color data and a .ren with vertex data and audio data)
They will be saved with the name that is in the scene name dialog box (so if you didn't change it after loading the scene just the scene name)

To view a render:
Type in a render name (no file suffix) and hit reset render (Make sure you hit reset render and not reset simulation, else it may override your render when you start running)
Hit Run Simulation
If sound not heard, try raising Volume Multiplier to an amount so that you can hear the sound (you can reset render/run to repeat)

How Scene Files Work:
4 Parameters were added per mesh, so after mesh name there is now the Young's Modulus, Poisson Ratio, Phi, Psi (and then it goes on to scale and density, position, ...).

Render Files:
We provide sample render files, a metalbox.ren and a metalspheres.ren
