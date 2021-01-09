import numpy as np

p1 = "build/out_malonno/prob_failure.asc"
p2 = "resources/prob_failure_ref.asc"

a = np.loadtxt(p1, skiprows=6)
b = np.loadtxt(p2, skiprows=6)

# a = a[b!= -9999]
# b = b[b!= -9999]

print(((a - b) ** 2).mean())