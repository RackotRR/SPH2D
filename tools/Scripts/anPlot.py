import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

x_label = ""
y_label = ""
target_path = ""

if len(sys.argv) > 1:
    target_path = sys.argv[1]
else:
    target_path = input("Enter target dir path: ")
if len(sys.argv) == 4:
    x_label = sys.argv[2]
    y_label = sys.argv[3]


files = os.listdir(target_path)
layer : pd.DataFrame
axes = plt.axes()

for filename in sorted(files):
    layer = pd.read_csv(os.path.join(target_path, filename), names=["x", "y"], delimiter=",")
    axes.plot(layer["x"], layer["y"], label=filename)
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    
axes.grid()
axes.legend()
plt.show()

