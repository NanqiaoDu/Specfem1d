import numpy as np 
import matplotlib.pyplot as plt

d=np.loadtxt('out.txt')
xg=d[:,0]
u=d[:,1]

plt.plot(xg,u)
plt.show()