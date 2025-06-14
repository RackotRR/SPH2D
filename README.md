# RRSPH
Implementation for 2D and 3D SPH simulation of water

![video_vx_0](https://user-images.githubusercontent.com/60754292/225715442-05b63d00-a091-4c36-adef-35da8d234468.gif)

Two-dimensional water simulation with weakly compressible smoothed particle hydrodynamics (WCSPH) method. This code is based on fortran code from ["Liu G.R., Liu M.B. - Smoothed Particle Hydrodynamics A Meshfree Particle Method 2003" book](https://www.worldscientific.com/worldscibooks/10.1142/5340).

### General
Here are two versions of the program: based on **OpenMP** (RRSPH_OMP) and based on **OpenCL** (RRSPH_CL). 
Generally they use the same algorithms except nearest neighbour particle search (NNPS) modules:
- RRSPH_OMP uses counting sort here;
- RRSPH_CL uses bitonic sort (so here the maximum number of particles should be a power of 2) and binary search.
Though SPH is meshfree method, NNPS here uses mesh in order to accelerate pairs searching (method's described in [my report](https://github.com/RackotRR/SPH2D/files/10994545/_._2021_._._._.pdf)).

Here are implemented dynamic boundaries based on Lennard-Jones potential. In several experiments I used them to simulate piston movement and generate waves by its means. Here is [my other report](https://github.com/RackotRR/SPH2D/files/10994558/_._2022_._._._._._._.pdf) on this topic.

All the input is divided into particles data and model params. Particle data is generated with Python scripts or with PicGen tool. Model params can be filled in manually (based on default-generated params) or by the means of SPH2DParamsGenerator tools (Windows only). 

RRSPH_CL can efficiently load it in runtime, create `clparams.h` header and compile its programs with a lot of compile-time substitutions. 

### Install and run

Configure, build and install project with cmake:

```
git clone https://github.com/RackotRR/SPH2D.git
git submodule update --init
cd RRSPH
cmake -S . -B {binary_path} -DCMAKE_INSTALL_PREFIX={install_path}
cmake --build {binary_path}
cmake --install {binary_path}
cd {install_path}/RRSPH
```

You have to install your project in order to use RRSPH_CL!
You can meet problems with SDL2 based projects on Windows. If so, specify path for its cache variables (see SPH2D_Drawer/cmake or SPH2D_PicGen/cmake) or update your PATH system variable.

If succeed, run particles generation script:

```
python ./Scripts/{script_name}
> Enter experiment name:
{experiment_name}
```

Then you have to copy Model params in order to be able to run an experiment:

```
Windows: copy default_experiment_params\ModelParams.json .\{experiment_name}\
Linux: cp default_experiment_params/ModelParams.json ./{experiment_name}/
```

Now you can start experiment with RRSPH_CL or RRSPH_OMP:
```
./RRSPH_CL
> Found experiments:
> [0] {experiment_name}: (0/1) data/dump layers
> Type experiment number you want to load:
0
```

### Using
You can use this project however you want. I hope it just can be helpful or just interesting to see.

It's recommended to compile the project with C++20 or later (you can use C++17, but some features can be disabled). GCC and msvc compilers supported. You can build and run it on Windows or Linux (macOS isn't tested).
Here are two main executables: RRSPH_OMP and RRSPH_CL. Installation of RRSPH_CL also copies its OpenCL code to destination folder.
There were experiments for several million particles powered by RRSPH_CL.
If you don't have OpenCL package, it won't be compiled and installed, so you'll have only RRSPH_OMP. You can also use RRSPH_OMP without OpenMP in single-threaded way.

There are also several tools: 
- FuncAtPoint: finds specified function value at point and plots it;
- WaterProfile: takes AnalysisParams.json file and plots water profile by space or by time with specified transformations;
- PartToGridConverter: converts particle data format into grid format by smoothing;
- [SPH2D_Drawer](https://github.com/RackotRR/SPH2D_Drawer): my visualization tool with heatmap and output customization by commands (requires SDL2 package installed);
- SPH2D_PicGen: tool for experiment generation from pixel image (requires SDL2 package installed);
- [SPH2DParamsGenerator](https://github.com/RackotRR/SPH2DParamsGenerator): my Model params editor (Windows only).

### Output
- If something goes wrong, there's log. It contains all experiment params, device info and experiment flow.
- Experiments take a lot of time, so program provides time estimate based on previous iterations.
- Every `params.save_time` program saves selected particles state into csv file (in separate thread) so you can read it later.
- Every `params.dump_time` program dumps all particles state into csv file (in separate thread) so you can start from this state later.

Here are some examples of water simulations:
![video_vy_0](https://user-images.githubusercontent.com/60754292/225718496-eba40340-dff1-415e-94d0-2c0ef6b56693.gif)
![r_39](https://user-images.githubusercontent.com/60754292/225718440-62defc68-e626-4ef6-8eb2-ceebafcd79b2.png)
![vy_39](https://user-images.githubusercontent.com/60754292/225718463-5491cd74-7171-4de0-bb5f-23520666eebb.png)
