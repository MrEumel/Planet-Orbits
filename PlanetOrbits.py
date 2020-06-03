import math
import numpy as np
import matplotlib.pyplot as plt
import time
#from astroquery.jplhorizons import Horizons

dt = 0.001
steps = 10000

G = 6.67 * 10 ** 20              #these do nothing for now, because G is set to 1 to simplify and I only factored in relative mass earth to sun 1:333000
m_x = 5.972 * 10 ** 24
m_sun = 1.989 * 10 ** 30

#Starting conditions

#obj = Horizons(id=1, location="@sun", epochs=Time("2017-01-01").jd, id_type='id').vectors()
#x_start = obj['x']
#vx_start = obj['v']

x_start = np.array([0.5, 0])
vx_start = np.array([0, 1.5])

xsun_start = np.array([0, 0])
vsun_start = np.array([0, 0])

#Starting arrays

x = np.zeros((2, steps))
vx = np.zeros((2, steps))
ax = np.zeros((2, steps))

r = np.zeros((1, steps))

xsun = np.zeros((2, steps))
vsun = np.zeros((2, steps))
asun = np.zeros((2, steps))

#Inserting starting conditions

x[:, 0] = x_start
xsun[:, 0] = xsun_start

#for i in range(0, steps):            #needed when sun at fixed position
#    xsun[:, i] = xsun_start

r[0, 0] = math.sqrt((x[0, 0]-xsun[0, 0])**2 + (x[1, 0]-xsun[1, 0])**2)
r3 = r[0, 0]**3

ax[:, 0] = - (x[:, 0]-xsun[:, 0])/r3
vx[:, 0] = vx_start + dt*ax[:, 0]

asun[:, 0] = - (xsun[:, 0]-x[:, 0])/(r3 * 333000)      #factor because m_earth to m_sun is 1:333000
vsun[:, 0] = vsun_start + dt*asun[:, 0]

#Calculating time development

for i in range(1, steps):
    x[:, i] = x[:, i-1] + dt*vx[:, i-1]
    xsun[:, i] = xsun[:, i - 1] + dt * vsun[:, i - 1]

    r[0, i] = math.sqrt((x[0, i]-xsun[0, i]) ** 2 + (x[1, i]-xsun[1, i]) ** 2)
    r3 = r[0, i]**3

    ax[:, i] = - (x[:, i] - xsun[:, i]) / r3
    vx[:, i] = vx[:, i - 1] + dt * ax[:, i - 1]
    asun[:, i] = - (xsun[:, i] - x[:, i]) / (r3 * 333000)          #factor because m_earth to m_sun is 1:333000
    vsun[:, i] = vsun[:, i - 1] + dt * asun[:, i - 1]

#Plotting trajectory

plt.plot(x[0], x[1])
plt.plot(xsun[0], xsun[1], 'ro')
plt.grid()
#plt.xlim(-1, 1)
#plt.ylim(-1, 1)
plt.show()

#print(time.process_time())


