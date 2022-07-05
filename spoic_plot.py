import numpy as np
import matplotlib.pyplot as plt

shift = 0
per = 16
data = np.loadtxt('RUN/spoic.out')

# plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D((2.0 + data[:,0]*np.cos(data[:,1]))*np.cos(data[:,2]),
#           (2.0 + data[:,0]*np.cos(data[:,1]))*np.sin(data[:,2]),
#              data[:,0]*np.sin(data[:,1]))

plt.figure()
plt.plot(data[shift::per,0]*np.cos(data[shift::per,1]),
         data[shift::per,0]*np.sin(data[shift::per,1]), ',')
