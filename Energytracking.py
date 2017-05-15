import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')





data = np.genfromtxt("Energy.txt")
# plt.plot(data[:,1],data[:,0], label="TOTAL E")
# plt.plot(data[:,1],data[:,2], label="KE")
plt.plot(data[:,1],data[:,3], label="PE")
# plt.yscale ('log')
plt.legend()
plt.show()

# data2 = np.genfromtxt("PEsheet.txt")
#
# ax.scatter(data2[:,1],data2[:,2],data2[:,0])
# plt.show()
