import matplotlib.pyplot as plt


with open("work2.csv") as f:
    data = f.readlines()

tmp = data.pop(0).split(',')
plt.xlabel(tmp[0])
plt.ylabel(tmp[1])
x = []
y = []
for i in data:
    tmp = list(map(float, i.split(',')))
    x.append(tmp[0])
    y.append(tmp[1])

plt.grid()
plt.plot(x, y)
plt.show()

