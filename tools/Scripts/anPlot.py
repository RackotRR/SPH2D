import matplotlib.pyplot as plt
import pandas as pd
import os

target_path = input("Enter target dir path: ")

files = os.listdir(target_path)
layer : pd.DataFrame

axes = plt.axes()

i = 0
for filename in sorted(files):
    layer = pd.read_csv(os.path.join(target_path, filename), names=["x", "y"], delimiter=",")
    axes.plot(layer["x"], layer["y"], label=filename)
    i += 1

axes.grid()
axes.legend()
plt.show()

