import numpy as np

eps = np.half(1)
while np.half(1) + eps > np.half(1):
    eps /= np.half(2)
print("Machine epsilon, half:", eps * 2)
print("Theoretical machine epsilon, half:", np.finfo(np.half).eps, '\n')


eps2 = np.single(1)
while np.single(1) + eps2 > np.single(1):
    eps2 /= np.single(2)
print("Machine epsilon, single:", eps2 * 2)
print("Theoretical machine epsilon, single:", np.finfo(np.single).eps, '\n')


eps3 = np.double(1)
while np.double(1) + eps3 > np.double(1):
    eps3 /= np.double(2)
print("Machine epsilon, double:", eps3 * 2)
print("Theoretical machine epsilon, double:", np.finfo(np.double).eps, '\n')

print("--------------------------------------------------------------\n")

dwarf = np.half(1)
pre_dwarf = np.half(1)
while dwarf > np.half(0):
    pre_dwarf = dwarf
    dwarf /= np.half(2)
print("Machine dwarf, half:", pre_dwarf)
print("Theoretical machine dwarf, half:", np.finfo(np.half).tiny, '\n')

dwarf2 = np.single(1)
pre_dwarf2 = np.single(1)
while dwarf2 > np.single(0):
    pre_dwarf2 = dwarf2
    dwarf2 /= np.single(2)
print("Machine dwarf, single:", pre_dwarf2)
print("Theoretical machine dwarf, single:", np.finfo(np.single).tiny, '\n')

dwarf3 = np.double(1)
pre_dwarf3 = np.double(1)
while dwarf3 > np.double(0):
    pre_dwarf3 = dwarf3
    dwarf3 /= np.double(2)
print("Machine dwarf, double:", pre_dwarf3)
print("Theoretical machine dwarf, double:", np.finfo(np.double).tiny, '\n')

print("--------------------------------------------------------------\n")

maximum = np.half(1)
pre_maximum = np.half(1)
while maximum != np.half(np.inf):
    pre_maximum = maximum
    maximum *= np.half(2)
print("Machine maximum, half:", pre_maximum)
print("Theoretical machine maximum, half:", np.finfo(np.half).max, '\n')

maximum2 = np.single(1)
pre_maximum2 = np.single(1)
while maximum2 != np.single(np.inf):
    pre_maximum2 = maximum2
    maximum2 *= np.single(2)
print("Machine maximum, single:", pre_maximum2)
print("Theoretical machine maximum, single:", np.finfo(np.single).max, '\n')

maximum3 = np.double(1)
pre_maximum3 = np.double(1)
while maximum3 != np.double(np.inf):
    pre_maximum3 = maximum3
    maximum3 *= np.double(2)
print("Machine maximum, double:", pre_maximum3)
print("Theoretical machine maximum, double:", np.finfo(np.double).max, '\n')
