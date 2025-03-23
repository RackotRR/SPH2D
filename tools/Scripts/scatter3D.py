import sys
import os
import json
import pandas as pd
import numpy as np

from vispy import scene
from vispy.color import Color
from vispy.scene import visuals
import vispy.io as io

def input_experiment_dir() -> str:
    script_path = os.path.abspath(sys.argv[0])
    current_dir = os.path.dirname(script_path)
    parent_dir = os.path.split(current_dir)[0]

    directory = os.fsencode(parent_dir)
    experiments = []
    i = 0
    for dir in os.listdir(directory):
        dirname = os.fsdecode(os.path.join(directory, dir))
        if not os.path.isdir(dirname):
            continue
        if not os.path.exists(os.path.join(dirname, "ParticleParams.json")):
            continue

        print(f"[{i}] {os.fsdecode(dir)}")
        experiments.append(dirname)
        i += 1

    if len(experiments) == 0:
        print("no experiments found")
        exit()
        
    idx = int(input("enter experiment idx:"))
    if idx >= 0 and idx < i:
        print(f"experiment {i}: " + experiments[idx])
        return experiments[idx]
    else:
        print("invalid input")
        exit()

def load_particle_params(experiment_dir : str):
    particle_params_path = os.path.join(experiment_dir, "ParticleParams.json")
    with open(particle_params_path) as particle_params_file:
        return json.load(particle_params_file)
    exit("can't open particle params file")

def load_time_layers(data_path : str) -> list:
    files = []

    if not os.path.exists(data_path):
        return files
    
    directory = os.fsencode(data_path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        files.append(os.path.join(data_path, filename))

    return files

def print_time_layers(filenames : list, prefix : str = ""):
    for i, filename in enumerate(filenames):
        print(f"[{i}] {prefix}_{filename}")

def select_time_layer(filenames : list):
    idx = int(input("enter file idx:"))
    if idx >= 0 and idx < len(filenames):
        print(f"file {idx}: " + filenames[idx])
        return filenames[idx]
    else:
        return None

def input_time_layers(experiment_dir):
    layers_origin = "data"
    layers = load_time_layers(os.path.join(experiment_dir, layers_origin))
    if len(layers) == 0: 
        layers_origin = "dump"
        layers = load_time_layers(os.path.join(experiment_dir, layers_origin))

    if len(layers) == 0:
        print("no data files")
        exit()

    return layers

def load_dump(filename : str):
    dump = pd.read_csv(filename, usecols=["x", "y", "z", "itype"])
    return dump

def plot_fluid(dump):
    fluid = dump[dump["itype"] == 2]
    fluid_pos = np.array(fluid[["x", "y", "z"]])
    fluid_scatter = visuals.Markers(
        pos=fluid_pos,
        size=5,
        antialias=0,
        face_color=(0, 0, 1, .5),
        spherical=True
    )
    return fluid_scatter

def plot_boundary(dump):
    boundary = dump[dump["itype"] == -2]
    boundary_pos = np.array(boundary[["x", "y", "z"]])
    boundary_scatter = visuals.Markers(
        pos=boundary_pos,
        size=5,
        antialias=0,
        face_color=(1, 1, 1, .5),
        spherical=True
    )
    return boundary_scatter

#///////////////////////////////////////////////////////////////////////////////////////////////////

experiment_dir = input_experiment_dir()
particle_params = load_particle_params(experiment_dir)
time_layers = input_time_layers(experiment_dir)
current_layer = 0

def reload_dump(new_layer_id : int):
    global dump
    global current_layer
    global fluid_widget
    global boundary_widget
    global view
    global time_layers
    
    if new_layer_id == current_layer:
        return
    
    print(f"layer {new_layer_id}")
    current_layer = new_layer_id
    dump = load_dump(time_layers[current_layer])
    
    fluid_widget.parent = None
    boundary_widget.parent = None

    fluid_widget = plot_fluid(dump)
    boundary_widget = plot_boundary(dump)

    view.add(fluid_widget)
    view.add(boundary_widget)

def next_layer():
    reload_dump(min(current_layer + 1, len(time_layers) - 1))

def prev_layer():
    reload_dump(max(current_layer - 1, 0))

def record_screenshot():
    base_path = os.path.join(experiment_dir, "screenshots")
    os.makedirs(base_path, exist_ok=True)
    io.write_png(os.path.join(base_path, f"{current_layer:03}.png"), canvas.render())

def record_video():
    base_path = os.path.join(experiment_dir, "videos")
    raw_path = os.path.join(base_path, "raw")

    os.makedirs(raw_path, exist_ok=True)
    
    for layer_id, _ in enumerate(time_layers):
        reload_dump(layer_id)
        io.write_png(os.path.join(raw_path, f"{layer_id:03}.png"), canvas.render())
    
    with open(os.path.join(base_path, "render.bat"), "w") as file:
        file.write("@echo off\n")
        framerate = len(time_layers) // 10 + 1
        path = os.path.join(os.path.abspath(raw_path), "%%3d.png")
        file.write(f"ffmpeg -framerate {framerate} -i {path} -c:v libx264 -pix_fmt yuv420p video.mp4")

def toggle_fullscreen():
    canvas.fullscreen = not canvas.fullscreen

keys = {
    "v": record_video,
    "c": record_screenshot,
    "right": next_layer,
    "left": prev_layer,
    "f11": toggle_fullscreen,
    "escape": "close"
}

# canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
canvas = scene.SceneCanvas(keys=keys, size=(800, 600), show=True)

# Set up a viewbox to display the cube with interactive arcball
view = canvas.central_widget.add_view()
view.bgcolor = "black"#'#efefef'
# view.camera = 'arcball'
view.camera = 'turntable'
view.camera.up = "+y"
view.padding = 100

dump = load_dump(time_layers[0])
print("dump loaded")
print("nfluid:", particle_params["nfluid"])
print("nvirt:", particle_params["nvirt"])
print("ntotal:", particle_params["ntotal"])

fluid_widget = plot_fluid(dump)
boundary_widget = plot_boundary(dump)
view.add(fluid_widget)
view.add(boundary_widget)


if __name__ == '__main__' and sys.flags.interactive == 0:
    canvas.app.run()
