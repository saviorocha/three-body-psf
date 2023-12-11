import math
import matplotlib.pyplot as plt
from scipy.constants import gravitational_constant as G
import numpy as np
from scipy.integrate import odeint
from matplotlib import animation

Tt = 365 * 360 * 24

mT = 5.97e24
mS = 1.98e30
mL = 7.36e22

xT_0 = 0
yT_0 = 149.6e9
xS_0 = 0
yS_0 = 0
xL_0 = 0
yL_0 = yT_0 + 384400000

vxT_0 = 30000  # vel média da terra ~30km/s
vyT_0 = 0
vxS_0 = 0    
vyS_0 = 0
vxL_0 = vxT_0 + 1000  # vel média da lua ~1km/s
vyL_0 = 0

def dSdt(S, t):
    xT, yT, xS, yS, xL, yL, vxT, vyT, vxS, vyS, vxL, vyL = S

    r12 = np.sqrt((xT-xS)**2 + (yT-yS)**2)  # terra-sol
    r13 = np.sqrt((xT-xL)**2 + (yT-yL)**2)  # terra-lua
    r23 = np.sqrt((xL-xS)**2 + (yL-yS)**2)  # lua-sol

    return [
        vxT,
        vyT,
        vxS,
        vyS,
        vxL,
        vyL,
        -G * (mS/r12**3 * (xT-xS) + mL/r13**3 * (xT-xL)),  # aceleração da terra
        -G * (mS/r12**3 * (yT-yS) + mL/r13**3 * (yT-yL)),
        -G * (mT/r12**3 * (xS-xT) + mL/r23**3 * (xS-xL)),  # aceleração do sol
        -G * (mT/r12**3 * (yS-yT) + mL/r23**3 * (yS-yL)),
        -G * (mT/r13**3 * (xL-xT) + mS/r23**3 * (xL-xS)),  # aceleração da lua
        -G * (mT/r13**3 * (yL-yT) + mS/r23**3 * (yL-yS))
    ]


t = np.linspace(0, 200000000, 1000)
sol = odeint(dSdt, y0=[xT_0, yT_0, xS_0, yS_0, xL_0, yL_0,
                       vxT_0, vyT_0, vxS_0, vyS_0, vxL_0, vyL_0], t=t)

t = sol.T
x1 = sol.T[0]
y1 = sol.T[1]
x2 = sol.T[2]
y2 = sol.T[3]
x3 = sol.T[4]
x3 = x3 + (10 * ( x3 - x1 ))
y3 = sol.T[5]
y3 = y3 + (10 * ( y3 - y1 ))


def animate(i):
    ln1.set_data([x1[i]], [y1[i]])
    ln2.set_data([x2[i]], [y2[i]])
    ln3.set_data([x3[i]], [y3[i]])


fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.grid()
ln1, = plt.plot([], [], 'bo', lw=3, markersize=6)
ln2, = plt.plot([], [], 'yo', lw=3, markersize=6)
ln3, = plt.plot([], [], 'ro', lw=3, markersize=6)
ax.set_ylim(-(yT_0 + 10**11), (yT_0 + 10**11))
ax.set_xlim(-(yT_0 + 10**11), (yT_0 + 10**11))
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
plt.title('Sistema Terra-Lua-Sol')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()
# ani.save('terraluasol.gif',writer='pillow',fps=30)
plt.close()

plt.plot(x1)
plt.plot(y1)
plt.legend(['x (m)', 'y (m)'], loc='upper left')

plt.xlabel('t (s)')
plt.title('Movimento da Terra')
plt.show()

