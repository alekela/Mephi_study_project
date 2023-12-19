import numpy as np
import matplotlib.pyplot as plt


error_values = list(map(np.log, [2.259583333333333e-05, 2.252349004975124e-05, 2.2518492553723106e-05, 2.2517993785310715e-05, 2.6919503900000035e-05]))
step_values = list(map(np.log, [0.05, 0.005, 0.0005, 0.00005, 0.000005]))

plt.plot(step_values, error_values)

plt.show()
