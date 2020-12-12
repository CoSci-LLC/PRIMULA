import numpy as np

p1 = "build/prob_failure.asc"
p2 = "resources/prob_failure_ref.asc"

a = np.loadtxt(p1, skiprows=6)
b = np.loadtxt(p2, skiprows=6)

print(((a - b) ** 2).mean())