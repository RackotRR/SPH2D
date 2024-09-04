import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

experiment_dir = input("Enter experiment_dir:")

with open(os.path.join(experiment_dir, "ParticleParams.json")) as particle_params_file:
    particle_params = json.load(particle_params_file)

if particle_params["dim"] != 3:
    print("expected dim = 3")
    exit()

ntotal = particle_params["ntotal"]
print("ntotal:", ntotal)

dump = pd.read_csv(os.path.join(experiment_dir, "dump", "0.csv"), usecols=["x", "y", "z"])

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

ax.scatter(dump["x"], dump["y"], dump["z"])
plt.show()