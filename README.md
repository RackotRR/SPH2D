# SPH2D
Implementation for 2D SPH simulation of water

![video_vx_0](https://user-images.githubusercontent.com/60754292/225715442-05b63d00-a091-4c36-adef-35da8d234468.gif)

Two-dimensional water simulation with smoothed particle hydrodynamics (SPH) method. This code is based on fortran code from ["Liu G.R., Liu M.B. - Smoothed Particle Hydrodynamics A Meshfree Particle Method 2003" book](https://www.worldscientific.com/worldscibooks/10.1142/5340).

### General
Here are two versions of the program: based on **OpenMP** (SPH2D_OMP) and based on **OpenCL** (SPH2D_CL). 
Generally they are the same except NNPS (nearest neighbour particle search) modules:
- SPH2D_OMP uses counting sort here;
- SPH2D_OCL uses bitonic sort (so here the maximum number of particles should be a power of 2) and binary search.
Though SPH is meshfree method, NNPS here uses mesh in order to accelerate pairs searching (method's described in [my report](https://github.com/RackotRR/SPH2D/files/10994545/_._2021_._._._.pdf)).

Here are implemented dynamic boundaries based on Lennard-Jones potential. In several experiments I used them to simulate piston movement and generate waves by its means. Here is [my other report](https://github.com/RackotRR/SPH2D/files/10994558/_._2022_._._._._._._.pdf) on this topic.

All the input generated in `Input.cpp` and `Params.h` files. Also you can manually fill in Params.json file to run experiment without recompilation. SPH2D_OCL will efficiently load it in runtime, create `clparams.h` header and compile its programs.

### Using
You can use this project however you want. I hope it just can be helpful or just interesting to see.

Project requires C++20. Tested compilation by g++ and msvc.
Here are two main executables: SPH2D_OMP and SPH2D_CL. Installation of SPH2D_CL also copies its OpenCL code to destination folder.
There were experiments for several millions particles powered by SPH2D_CL.

There are also several tools: 
- FuncAtPoint: finds specified function value at point and plots it;
- WaterProfile: takes AnalysisParams.json file and plots water profile by space or by time with specified transformations;
- [SPH2D_Drawer](https://github.com/RackotRR/SPH2D_Drawer): my visualizetion tool with heatmap and output customization by commands.

### Output
- If something goes wrong, there's log. It contains all experiment params, device info and experiment flow.
- Experiments take a lot of time, so program provides time estimate every `Params::step_time_estimate` steps based on previous iterations.
- Every `params.save_step` steps program dumps all particles state into a file (in separate thread) with format string, so you can read it later.
- To find out all the parameters used to perform experiment program also dumps generated `clparams.h` file, copy of that is used in calculations by SPH2D_OCL.

Here are some examples of water simulations:
![video_vy_0](https://user-images.githubusercontent.com/60754292/225718496-eba40340-dff1-415e-94d0-2c0ef6b56693.gif)
![r_39](https://user-images.githubusercontent.com/60754292/225718440-62defc68-e626-4ef6-8eb2-ceebafcd79b2.png)
![vy_39](https://user-images.githubusercontent.com/60754292/225718463-5491cd74-7171-4de0-bb5f-23520666eebb.png)
