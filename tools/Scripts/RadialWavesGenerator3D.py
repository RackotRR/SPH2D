import os
from shutil import copyfile
import sys
import json
import numpy as np
import csv

# script fields:
#mandatory:
# - experiment_name
# - x_mingeom
# - x_maxgeom
# - y_mingeom
# - y_maxgeom
# - z_mingeom
# - z_maxgeom
# - delta
# - ntotal
# - nfluid
# - nvirt
#optional:
# - rho0
# - nwm_particles_start
# - nwm_particles_end
# - nwm
# - boundary_layers_num
# - use_chess_order
# - depth
# - x_fluid_min
# - x_fluid_max
# - y_fluid_min
# - y_fluid_max
# - z_fluid_min
# - z_fluid_max
# - x_fluid_particles
# - y_fluid_particles
# - z_fluid_particles
# - x_boundary_left
# - x_boundary_right
# - x_boundary_center
# - y_boundary_bottom
# - y_boundary_top
# - z_boundary_near
# - z_boundary_far
# - boundary_delta
# - boundary_separation

class nwm:
    no_waves = 0
    relaxation_zone = 1
    dynamic_boundaries = 2
    impulse_method = 3
    wall_disappear = 4

class particle_type:
    boundary = -2
    non_existent = 0
    fluid = 2
    

def fill_in_common_params():
    global params
    params["dim"] = 3
    params["x_fluid_min"] = -1.0
    params["y_fluid_min"] = 0.0    
    params["x_fluid_max"] = 1.0
    params["y_fluid_max"] = 0.5
    params["z_fluid_min"] = -1.0
    params["z_fluid_max"] = 1.0
    
    params["delta"] = 0.01
    params["boundary_delta"] = params["delta"] * 1.0
    
    boundary_separation = 2 * params["delta"]
    params["boundary_separation"] = boundary_separation
    params["x_boundary_left"] = params["x_fluid_min"] - boundary_separation
    params["y_boundary_bottom"] = params["y_fluid_min"] - boundary_separation
    params["x_boundary_right"] = params["x_fluid_max"] + boundary_separation
    params["y_boundary_top"] = 1.0 + boundary_separation
    params["z_boundary_near"] = params["z_fluid_min"] - boundary_separation
    params["z_boundary_far"] = params["z_fluid_max"] + boundary_separation

    params["nwm"] = nwm.dynamic_boundaries

    params["x_mingeom"] = -6.0
    params["y_mingeom"] = -1.0
    params["z_mingeom"] = -6.0
    params["x_maxgeom"] = 6.0
    params["y_maxgeom"] = 2.0
    params["z_maxgeom"] = 6.0

    params["depth"] = params["y_fluid_max"] - params["y_fluid_min"]

    params["rho0"] = 1000.0
    params["target_eos_sound"] = 40.0
    params["nwm_radius"] = 0.15
    params["nwm_inner_radius"] = 0.09

    params["enable_fluid"] = True


def is_in_circle(x : float, z : float, R):
    return x * x + z * z < R * R

def p_hydrostatic(y):
    h = params["y_fluid_max"] - y
    if h > 0:
        return params["rho0"] * 9.8 * h
    else:
        return 0
    
def rho_eos(y):
    """ density from equation of state """
    import math

    p = p_hydrostatic(y)
    rho0 = params["rho0"]
    gamma = 7.0
    c = params["target_eos_sound"]
    B = (c ** 2) * rho0 / gamma
    return rho0 * math.pow(p / B + 1, 1.0 / gamma)
        
# - expects:
# params.x_fluid_particles
# params.y_fluid_particles
# params.x_fluid_min
# params.y_fluid_min
# params.delta
def generate_fluid_particles(r, v, rho, p, itype, start_i, fill : bool):
    if (params["enable_fluid"] == False):
       return 0

    x_fluid_min = params["x_fluid_min"]
    y_fluid_min = params["y_fluid_min"]
    z_fluid_min = params["z_fluid_min"]
    delta = params["delta"]
    x_fluid_particles = int((params["x_fluid_max"] - params["x_fluid_min"]) / delta)
    y_fluid_particles = int((params["y_fluid_max"] - params["y_fluid_min"]) / delta)
    z_fluid_particles = int((params["z_fluid_max"] - params["z_fluid_min"]) / delta)
    
    print("max x_fluid_particles: " + str(x_fluid_particles))
    print("max y_fluid_particles: " + str(y_fluid_particles))
    print("z_fluid_particles: " + str(z_fluid_particles))

    delta = params["delta"]
    i = start_i

    for x_i in range(x_fluid_particles):
        for y_i in range(y_fluid_particles):
            for z_i in range(z_fluid_particles):
                x = x_fluid_min + x_i * delta
                y = y_fluid_min + y_i * delta
                z = z_fluid_min + z_i * delta
               
                if is_in_circle(x, z, params["nwm_radius"]):
                    continue

                if fill:
                    r[0, i] = x
                    r[1, i] = y
                    r[2, i] = z

                    v[0, i] = 0.0
                    v[1, i] = 0.0
                    v[2, i] = 0.0
                    
                    rho[i] = rho_eos(y)
                    p[i] = 0.0
                    itype[i] = particle_type.fluid

                i += 1

    return i

# expects
# params.x_boundary_left
# params.x_boundary_right
# params.y_boundary_bottom
# params.y_boundary_top
# params.boundary_delta
# params.boundary_layers
# params.use_chess_order
# fills in:
# params.nwm_particles_start
# params.nwm_particles_end
# params.nwm
def generate_virt_particles(r, v, rho, p, itype, start_i, fill : bool):
    x_boundary_left = params["x_boundary_left"]
    x_boundary_right = params["x_boundary_right"]
    y_boundary_top = params["y_boundary_top"]
    y_boundary_bottom = params["y_boundary_bottom"]
    z_boundary_far = params["z_boundary_far"]
    z_boundary_near = params["z_boundary_near"]

    boundary_delta = params["boundary_delta"]

    x_range = (x_boundary_right - x_boundary_left) / boundary_delta
    y_range = (y_boundary_top - y_boundary_bottom) / boundary_delta
    z_range = (z_boundary_far - z_boundary_near) / boundary_delta

    i = start_i

    if params["nwm"] == nwm.dynamic_boundaries:
        params["nwm_particles_start"] = i
    x_fluid_min = params["x_fluid_min"]
    y_fluid_min = params["y_fluid_min"]
    z_fluid_min = params["z_fluid_min"]
    x_fluid_particles = int((params["x_fluid_max"] - params["x_fluid_min"]) / boundary_delta)
    y_fluid_particles = int((y_boundary_top - y_boundary_bottom) / boundary_delta)
    z_fluid_particles = int((params["z_fluid_max"] - params["z_fluid_min"]) / boundary_delta)
    R_outer = params["nwm_radius"]
    R_inner = params["nwm_inner_radius"]
    for x_i in range(x_fluid_particles):
        for y_i in range(y_fluid_particles):
            for z_i in range(z_fluid_particles):
                x = x_fluid_min + x_i * boundary_delta
                y = y_boundary_bottom + y_i * boundary_delta + boundary_delta
                z = z_fluid_min + z_i * boundary_delta

                if is_in_circle(x, z, R_outer) and not is_in_circle(x, z, R_inner):
                    if fill:
                        r[0, i] = x
                        r[1, i] = y
                        r[2, i] = z
                    i += 1
    if params["nwm"] == nwm.dynamic_boundaries:
        params["nwm_particles_end"] = i
    print("nwm particles:", i - start_i)
        
    # left wall
    for z_i in range(int(z_range) + 1):
        for y_i in range(int(y_range)):
            if fill:
                r[0, i] = x_boundary_left
                r[1, i] = y_boundary_bottom + y_i * boundary_delta
                r[2, i] = z_boundary_near + z_i * boundary_delta
            i = i + 1

    # right wall
    for z_i in range(int(z_range) + 1):
        for y_i in range(int(y_range)):
            if fill:
                r[0, i] = x_boundary_right
                r[1, i] = y_boundary_bottom + y_i * boundary_delta
                r[2, i] = z_boundary_near + z_i * boundary_delta
            i = i + 1
                
    # near wall
    for x_i in range(1, int(x_range)):
        for y_i in range(int(y_range)):
            if fill:
                r[0, i] = x_boundary_left + x_i * boundary_delta
                r[1, i] = y_boundary_bottom + y_i * boundary_delta
                r[2, i] = z_boundary_near
            i = i + 1
                
    # far wall
    for x_i in range(1, int(x_range)):
        for y_i in range(int(y_range)):
            if fill:
                r[0, i] = x_boundary_left + x_i * boundary_delta 
                r[1, i] = y_boundary_bottom + y_i * boundary_delta
                r[2, i] = z_boundary_far
            i = i + 1

    # ground
    for z_i in range(int(z_range)):
        for x_i in range(int(x_range)):
            if fill:
                r[0, i] = x_boundary_left + x_i * boundary_delta
                r[1, i] = y_boundary_bottom
                r[2, i] = z_boundary_near + z_i * boundary_delta
            i = i + 1

    nvirt = i
    
    # fill in other particles data
    if fill:
        for i in range(start_i, nvirt):
            v[0, i] = 0.0
            v[1, i] = 0.0
            v[2, i] = 0.0
            rho[i] = rho_eos(r[1, i])
            p[i] = 0.0
            itype[i] = particle_type.boundary

    return nvirt


def generate_particles_data():
    params["nvirt"] = nvirt = generate_virt_particles(None, None, None, None, None, 0, False)
    params["nfluid"] = nfluid = generate_fluid_particles(None, None, None, None, None, 0, False)
    params["ntotal"] = ntotal = nvirt + nfluid
    
    print("ntotal:" + str(ntotal))
    print("nvirt:" + str(nvirt))
    print("nfluid:" + str(nfluid))

    r = np.zeros([3, ntotal], dtype=float)
    v = np.zeros([3, ntotal], dtype=float)
    rho = np.zeros(ntotal, dtype=float)
    p = np.zeros(ntotal, dtype=float)
    itype = np.zeros(ntotal, dtype=int)

    generate_fluid_particles(r, v, rho, p, itype, 0, True)
    generate_virt_particles(r, v, rho, p, itype, nfluid, True)
    
    dump_dir = os.path.join(experiment_dir, "dump")
    if not os.path.exists(dump_dir):
        os.mkdir(dump_dir)
    dump_name = "0.csv"
    file = open(os.path.join(dump_dir, dump_name), "w", newline='')
    writer = csv.writer(file)
    header = ["x", "y", "z", "itype", "vx", "vy", "vz", "rho", "p"]
    writer.writerow(header)
    for i in range(ntotal):
        writer.writerow([r[0, i], r[1, i], r[2, i], itype[i], v[0, i], v[1, i], v[2, i], rho[i], p[i]])
    file.close()
        
def generate_project():
    fill_in_common_params()
    generate_particles_data()

    copyfile(script_path, os.path.join(experiment_dir, "Generator.py"))
    
    file = open(os.path.join(experiment_dir, "ParticleParams.json"), "w")
    json.dump(params, file, indent=2)
    file.close()
    

def isdigit(line : str):
    for c in line:
        if c.isnumeric():
            return True
    return False

def find_experiment_directory(directory : str):
    if not os.path.exists(directory):
        return directory

    index = directory.rfind("_")
    if index + 3 == len(directory):
        number_str : str = directory[index + 1:]
        if isdigit(number_str):
            number = int(number_str)
            next_experiment_name = directory[0 : index + 1] + "{:02}".format(number + 1)
            return find_experiment_directory(next_experiment_name)
        
    return find_experiment_directory(directory + "_01")


print("Enter experiment name: ")
experiment_name = input("> ")
script_path = os.path.abspath(sys.argv[0])
current_dir = os.path.dirname(script_path)
parent_dir = os.path.split(current_dir)[0]
experiment_dir = os.path.join(parent_dir, experiment_name)

if os.path.exists(experiment_dir):
    experiment_dir = find_experiment_directory(experiment_dir)
os.mkdir(experiment_dir)
experiment_name = os.path.basename(experiment_dir)

params = {}
params["experiment_name"] = experiment_name
generate_project()

