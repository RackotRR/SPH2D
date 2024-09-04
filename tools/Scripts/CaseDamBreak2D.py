import os
from shutil import copyfile
import sys
import json
import numpy
import math
import csv
import string

# script fields:
#mandatory:
# - experiment_name
# - x_mingeom
# - x_maxgeom
# - y_mingeom
# - y_maxgeom
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
# - x_fluid_particles
# - y_fluid_particles
# - x_boundary_left
# - x_boundary_right
# - x_boundary_center
# - y_boundary_bottom
# - y_boundary_top
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
    params["dim"] = 2

    params["x_fluid_min"] = 0.0
    params["y_fluid_min"] = 0.0    
    params["x_fluid_max"] = 1.0
    params["y_fluid_max"] = 2.0
    
    params["delta"] = 0.0025
    params["boundary_delta"] = params["delta"] * 0.5
    
    params["x_fluid_particles"] = int((params["x_fluid_max"] - params["x_fluid_min"]) / params["delta"])
    params["y_fluid_particles"] = int((params["y_fluid_max"] - params["y_fluid_min"]) / params["delta"])

    boundary_separation = params["delta"]
    params["boundary_separation"] = boundary_separation
    params["x_boundary_left"] = params["x_fluid_min"] - 3 * boundary_separation
    params["y_boundary_bottom"] = params["y_fluid_min"] - 3 * boundary_separation
    params["x_boundary_right"] = 4.0 + boundary_separation
    params["y_boundary_top"] = 3.0 + boundary_separation
    params["x_boundary_center"] = params["x_fluid_max"] + 3 * boundary_separation
    params["boundary_layers_num"] = 3
    params["use_chess_order"] = True

    params["nwm"] = nwm.no_waves

    params["x_mingeom"] = -1.0
    params["y_mingeom"] = -1.0
    params["x_maxgeom"] = 4.5
    params["y_maxgeom"] = 3.5

    params["depth"] = params["y_fluid_max"] - params["y_fluid_min"]

    params["rho0"] = 1000.0
        
# - expects:
# params.x_fluid_particles
# params.y_fluid_particles
# params.x_fluid_min
# params.y_fluid_min
# params.delta
def generate_fluid_particles(r, v, rho, p, itype, start_i):
    x_fluid_min = params["x_fluid_min"]
    y_fluid_min = params["y_fluid_min"]
    x_fluid_particles = params["x_fluid_particles"]
    y_fluid_particles = params["y_fluid_particles"]
    
    print("x_fluid_particles: " + str(x_fluid_particles))
    print("y_fluid_particles: " + str(y_fluid_particles))

    rho0 = params["rho0"]
    delta = params["delta"]
    i = start_i

    for x_i in range(x_fluid_particles):
        for y_i in range(y_fluid_particles):
            r[0, i] = x_fluid_min + x_i * delta
            r[1, i] = y_fluid_min + y_i * delta
            v[0, i] = 0.0
            v[1, i] = 0.0
            rho[i] = rho0
            p[i] = 0.0
            itype[i] = particle_type.fluid

            i += 1

    return i

# expects
# params.x_boundary_left
# params.x_boundary_right
# params.x_boundary_center
# params.y_boundary_bottom
# params.y_boundary_top
# params.boundary_delta
# params.boundary_layers
# params.use_chess_order
# fills in:
# params.nwm_particles_start
# params.nwm_particles_end
# params.nwm
def generate_virt_particles(fill : bool, r, v, rho, p, itype, start_i):
    nvirt = 0
    x_boundary_left = params["x_boundary_left"]
    x_boundary_right = params["x_boundary_right"]
    x_boundary_center = params["x_boundary_center"]
    y_boundary_top = params["y_boundary_top"]
    y_boundary_bottom = params["y_boundary_bottom"]

    boundary_delta = params["boundary_delta"]
    layers = params["boundary_layers_num"]
    use_chess_order = params["use_chess_order"]

    x_range = (x_boundary_right - x_boundary_left) / boundary_delta
    y_range = (y_boundary_top - y_boundary_bottom) / boundary_delta

    i = start_i

    # left wall
    if params["nwm"] == nwm.dynamic_boundaries:
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
            i = i + 1
    if params["nwm"] == nwm.dynamic_boundaries:
        params["nwm_particles_end"] = i

    # right wall
    for y_i in range(int(y_range)):
        for layer in range(layers):
            if fill:
                offset_x = layer * boundary_delta
                offset_y = 0.0
                if use_chess_order:
                    offset_y = -(layer % 2) * boundary_delta * 0.5
                r[0, i] = x_boundary_right + offset_x
                r[1, i] = y_boundary_bottom + y_i * boundary_delta + offset_y
            i = i + 1

    # ground
    for x_i in range(int(x_range)):
        for layer in range(layers):
            if fill:
                offset_x = 0.0
                offset_y = -layer * boundary_delta
                if use_chess_order:
                    offset_x = -(layer % 2) * boundary_delta * 0.5
                r[0, i] = x_boundary_left + x_i * boundary_delta + offset_x
                r[1, i] = y_boundary_bottom + offset_y
            i = i + 1


    # center wall
    if params["nwm"] == nwm.wall_disappear:
        params["nwm_particles_start"] = i
        for y_i in range(int(y_range)):
            for layer in range(layers):
                if fill:
                    offset_x = layer * boundary_delta
                    offset_y = 0.0
                    if use_chess_order:
                        offset_y = -(layer % 2) * boundary_delta * 0.5
                    r[0, i] = x_boundary_center + offset_x
                    r[1, i] = y_boundary_bottom + y_i * boundary_delta + offset_y
                i = i + 1
        params["nwm_particles_end"] = i

    nvirt = i
    
    rho0 = params["rho0"]

    # fill in other particles data
    if fill:
        for i in range(start_i, nvirt):
            v[0, i] = 0.0
            v[1, i] = 0.0
            rho[i] = rho0
            p[i] = 0.0
            itype[i] = particle_type.boundary

    return nvirt


def generate_particles_data():
    x_fluid_particles = params["x_fluid_particles"]
    y_fluid_particles = params["y_fluid_particles"]
    params["nvirt"] = nvirt = generate_virt_particles(False, None, None, None, None, None, 0)
    params["nfluid"] = nfluid = x_fluid_particles * y_fluid_particles
    params["ntotal"] = ntotal = nvirt + nfluid
    
    print("ntotal:" + str(ntotal))
    print("nvirt:" + str(nvirt))
    print("nfluid:" + str(nfluid))

    r = numpy.zeros([2, ntotal], dtype=float)
    v = numpy.zeros([2, ntotal], dtype=float)
    rho = numpy.zeros(ntotal, dtype=float)
    p = numpy.zeros(ntotal, dtype=float)
    itype = numpy.zeros(ntotal, dtype=int)

    generate_fluid_particles(r, v, rho, p, itype, 0)
    generate_virt_particles(True, r, v, rho, p, itype, nfluid)
    
    dump_dir = os.path.join(experiment_dir, "dump")
    if not os.path.exists(dump_dir):
        os.mkdir(dump_dir)
    dump_name = "0.csv"
    file = open(os.path.join(dump_dir, dump_name), "w", newline='')
    writer = csv.writer(file)
    header = ["x", "y", "itype", "vx", "vy", "rho", "p"]
    writer.writerow(header)
    for i in range(ntotal):
        writer.writerow([r[0, i], r[1, i], itype[i], v[0, i], v[1, i], rho[i], p[i]])
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

