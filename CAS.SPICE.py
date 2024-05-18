import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import glob

spice.kclear()
mk = glob.glob('/home/cemaxim/cosp_1000/extras/mk/*(copy).tm')
spice.furnsh(mk)
step = 4000
utc = ['Jan 1, 2005','Sep 14, 2005']

etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])


print("ET One: {}, ET Two: {}".format(etOne, etTwo))
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]
positions, lightimes = spice.spkpos('Cassini', times, 'J2000', 'NONE', 'SATURN BARYCENTER')

spice.kclear()

positions = positions.T
fig = plt.figure(figsize=(9,9))
ax = fig.add_subplot(111,projection='3d')
ax.plot(positions[0], positions[1],positions[2])
ax.plot([0],[0],[0], 'o', markersize = 10, color = 'gold')
"""
ax.axes.set_xlim3d(left=-1e6, right=1e6) 
ax.axes.set_ylim3d(bottom=-1e6, top=1e6) 
ax.axes.set_zlim3d(bottom=-1e6, top=1e6) 
"""
plt.show
