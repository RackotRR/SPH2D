import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

N_PERIODS = 3
PERIOD = 1.

def gen_sin_wave():
    COUNT = 1000
    x = np.linspace(0, N_PERIODS * PERIOD, COUNT)
    y = np.sin(x)

    diff = np.random.normal(size=COUNT, scale=0.1)
    return x, y + diff

peaks = []

def smooth(start : float, end : float, x, y):
    X = x[(x >= start) & (x <= end)]
    Y = y[(x >= start) & (x <= end)]
    parabolic_params = np.polyfit(X, Y, 2)
    print(parabolic_params)

    parabola = np.poly1d(parabolic_params)
    parabola_x = np.linspace(start, end, 1000)
    parabola_y = parabola(parabola_x)
    plt.plot(parabola_x, parabola_y)

    peak_x = -0.5 * parabolic_params[1] / parabolic_params[0]
    print("peak_x:", peak_x)
    print("peak_y:", parabola(peak_x))

    plt.scatter([peak_x], [parabola(peak_x)])
    peaks.append(peak_x)

   
df = pd.read_csv("../experiment_name/analysis/space_at_0.csv", names=["x", "y"])

# x, y = gen_sin_wave()
x = df["x"]
y = df["y"]
plt.plot(x, y)

OFFSET = 1.3
WIDTH = 0.3
for k in range(N_PERIODS):
    s = OFFSET + k * PERIOD
    smooth(s, s + WIDTH * (1 + k * 0), x, y)

count = 0
sum_half_L = 0
for prev, next in zip(peaks[:-1], peaks[1:]):
    half_L = next - prev
    sum_half_L += half_L
    count += 1
    print(f"L_{count}:", half_L)
print("Avg L:", sum_half_L / count)


plt.show()