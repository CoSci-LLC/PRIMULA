import numpy as np

p1 = "build/Study_Test_UNIPV/prob_failure.asc"
p2 = "../PRIMULATests/Tests/Study_Test_UNIPV/reference_probls.asc"

a = np.loadtxt(p1, skiprows=6)
b = np.loadtxt(p2, skiprows=6)

# a = a[b!= -9999]
# b = b[b!= -9999]

print(((a - b) ** 2).mean())