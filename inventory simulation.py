import numpy as np
import random
import matplotlib.pyplot as plt

n = 10
s = 0.5
S = 2
demand = []
replenish = []
x = [0]
y = [-s]
lambdas = np.array([1,2])
p = np.array([0.5,0.5])
for i in range(n):
    demand.append(random.uniform(0,1))
    if x[-1] < s:
        y.append(S - s)
        replenish.append(S - x[-1])
        x.append(max(S - demand[-1],0))
    else:
        y.append(x[-1] - s)
        replenish.append(0)
        x.append(max(x[-1] - demand[-1],0))
plt.plot(x)
plt.plot(y)
plt.plot(replenish)
plt.legend(['inventory','excess','replenish'])
plt.show()

