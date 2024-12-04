import os
from shutil import copyfile
import sys
import json
import numpy
import csv

class particle_type:
    boundary = -2
    non_existent = 0
    fluid = 2

def fill_in_common_params():
    # spacing between fluid particles
    params["delta"] = 0.005

    # Simulation domain
    params["x_mingeom"] = -1
    params["y_mingeom"] = -1
    params["x_maxgeom"] = 80
    params["y_maxgeom"] = 50

    # Fill in rectangle with fluid particles: 
    # [x_fluid_min; x_fluid_max)
    # [y_fluid_min; y_fluid_max)
    params["x_fluid_min"] = params["delta"] * 2
    params["y_fluid_min"] = beach_y(0)
    params["x_fluid_max"] = 100
    params["y_fluid_max"] = params["delta"] * 100
    
    # spacing between boundary particles
    params["boundary_delta"] = params["delta"]
    # spacing between boundary and fluid particles (makes sense when SBT=REPULSIVE)
    params["boundary_separation"] = params["delta"]

    # Left boundary wall: x_boundary <= x_boundary_left (piston position)
    params["x_boundary_left"] = 0
    # Boundary ground: y_boundary <= y_boundary_bottom
    params["y_boundary_bottom"] = params["y_fluid_min"] - params["boundary_separation"]    
    # Boundary walls max y coordinate
    params["y_boundary_top"] = 5
    # Instead of right wall here is beach_y function

    # Additional boundary layers will go away from fluid particles
    params["boundary_layers_num"] = 1
    # Chess order for boundary layers
    params["use_chess_order"] = True

    # Boundary ground parameters
    params["x_boundary_start"] = -0.5
    params["x_boundary_end"] = 75

    # Left wall will be marked as piston
    params["use_nwm_particles"] = True

    # Optional depth parameter
    params["depth"] = params["y_fluid_max"] - params["y_fluid_min"]

    # Target fluid density
    params["rho0"] = 1000.0

    # Target sound velocity
    params["target_eos_sound"] = 50.11347231455704
    
def beach_y(x):
    y0 = 0
    X_BEACH = 70
    X_MAX = 75
    Y_MAX = 2
    if x < X_BEACH: # strait line
        return y0
    else: # linear beach
        x0 = X_BEACH
        x1 = X_MAX
        y1 = Y_MAX
        return y0 + (x - x0) / (x1 - x0) * (y1 - y0)

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
        
def generate_fluid_particles(fill : bool, r, v, rho, itype, start_i):
    x_fluid_min = params["x_fluid_min"]
    x_fluid_max = params["x_fluid_max"]
    y_fluid_max = params["y_fluid_max"]

    delta = params["delta"]
    i = start_i

    X = numpy.arange(x_fluid_min, x_fluid_max, delta)

    # ground
    for x in X:
        y = beach_y(x) + 2 * delta 
        while y < y_fluid_max:
            if fill:
                r[0, i] = x
                r[1, i] = y
                rho[i] = rho_eos(y)
            y += delta
            i = i + 1

    if fill:
        v[0, start_i:i] = 0.0
        v[1, start_i:i] = 0.0
        itype[start_i:i] = particle_type.fluid

    return i - start_i

def generate_virt_particles(fill : bool, r, v, rho, itype, start_i):
    x_boundary_left = params["x_boundary_left"]
    y_boundary_top = params["y_boundary_top"]
    y_boundary_bottom = params["y_boundary_bottom"]

    boundary_delta = params["boundary_delta"]
    layers = params["boundary_layers_num"]
    use_chess_order = params["use_chess_order"]

    y_range = (y_boundary_top - y_boundary_bottom) / boundary_delta

    i = start_i

    # left wall
    if params["use_nwm_particles"]:
        params["nwm_particles_start"] = i
    for y_i in range(int(y_range)):
        for layer in range(layers):
            if fill:
                offset_x = -layer * boundary_delta
                offset_y = 0.0
                if use_chess_order:
                    offset_y = -(layer % 2) * boundary_delta * 0.5
                r[0, i] = x_boundary_left + offset_x
                r[1, i] = y_boundary_bottom + y_i * boundary_delta + offset_y
                rho[i] = rho_eos(r[1, i])
            i = i + 1
    if params["use_nwm_particles"]:
        params["nwm_particles_end"] = i

    # ground and right wall
    X = numpy.arange(params["x_boundary_start"], params["x_boundary_end"], boundary_delta)
    for x in X:
        for layer in range(layers):
            if fill:
                offset_x = 0.0
                if use_chess_order:
                    offset_x = -(layer % 2) * boundary_delta * 0.5

                x += offset_x
                y = beach_y(x) - layer * boundary_delta + offset_y
                r[0, i] = x
                r[1, i] = y
                rho[i] = rho_eos(y)
            i = i + 1
    
    # fill in other particles data
    if fill:
        v[0, start_i:i] = 0.0
        v[1, start_i:i] = 0.0
        itype[start_i:i] = particle_type.boundary

    return i - start_i

def save_dump(r, itype, v, rho):
    dump_dir = os.path.join(experiment_dir, "dump")
    if not os.path.exists(dump_dir):
        os.mkdir(dump_dir)
    dump_name = "0.csv"
    file = open(os.path.join(dump_dir, dump_name), "w", newline='')
    writer = csv.writer(file)
    header = ["x", "y", "itype", "vx", "vy", "rho"]
    writer.writerow(header)
    for i in range(r.shape[1]):
        writer.writerow([r[0, i], r[1, i], itype[i], v[0, i], v[1, i], rho[i]])
    file.close()

def generate_particles_data():
    # count particles
    params["nvirt"] = nvirt = generate_virt_particles(False, None, None, None, None, 0)
    params["nfluid"] = nfluid = generate_fluid_particles(False, None, None, None, None, 0)
    params["ntotal"] = ntotal = nvirt + nfluid
    
    print("ntotal: ", ntotal)
    print("nvirt: ", nvirt)
    print("nfluid: ", nfluid)

    # allocate memory
    r = numpy.zeros([2, ntotal], dtype=float)
    v = numpy.zeros([2, ntotal], dtype=float)
    rho = numpy.zeros(ntotal, dtype=float)
    itype = numpy.zeros(ntotal, dtype=int)

    # fill in arrays
    generate_fluid_particles(True, r, v, rho, itype, 0)
    generate_virt_particles(True, r, v, rho, itype, nfluid)
    
    save_dump(r, itype, v, rho)
        
def save_params():
    file = open(os.path.join(experiment_dir, "ParticleParams.json"), "w")
    json.dump(params, file, indent=2)
    file.close()

def generate_project():
    fill_in_common_params()
    generate_particles_data()

    copyfile(script_path, os.path.join(experiment_dir, "Generator.py"))
    save_params()
    

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
print("save as", experiment_name)
