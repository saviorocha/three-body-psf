import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import gravitational_constant as Go
from scipy.integrate import odeint
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter

Tt = 365 * 360 * 24

G = Go*10000

mT = 5.97e24
mS = 1.98e30
mL = 7.36e22

xT_0 = 1e-10
yT_0 = 149.6e6
xS_0 = 1e-10
yS_0 = 1e-10
xL_0 = 1e-10
yL_0 = yT_0 + 384400

vxT_0 = 30 * Tt
# vxT_0 = (2 * math.pi * yT_0) / Tt * 10**2
vyT_0 = 1e-10
vxS_0 = 1e-10
vyS_0 = 1e-10
vxL_0 = vxT_0 + Tt 
# vxL_0 = vxT_0 + (2 * math.pi * yL_0) / Tt * 10**2
vyL_0 = 1e-10

def dSdt(S, t):
    xT, yT, xS, yS, xL, yL, vxT, vyT, vxS, vyS, vxL, vyL = S

    #vyT = vyT*10**2
    #vyS = vyS*10**2
    #vyL = vyL*10**2

    # vxT = vyT*10**2
    # vxS = vyS*10**2
    # vxL = vyL*10**2

    r12 = np.sqrt((xS-xT)**2 + (yS-yT)**2) # terra-sol
    r13 = np.sqrt((xL-xT)**2 + (yL-yT)**2) # terra-lua
    r23 = np.sqrt((xS-xL)**2 + (yS-yL)**2) # lua-sol

    return [
        vxT,
        vyT,
        vxS,
        vyS,
        vxL,
        vyL,
        G * (mS/r12**3 * (xS-xT) + mL/r13**3 * (xL-xT)), # aceleração da
        G * (mS/r12**3 * (yS-yT) + mL/r13**3 * (yL-yT)),
        G * (mT/r12**3 * (xT-xS) + mL/r23**3 * (xL-xS)),
        G * (mT/r12**3 * (yT-yS) + mL/r23**3 * (yL-yS)),
        G * (mT/r13**3 * (xT-xL) + mS/r23**3 * (xS-xL)),
        G * (mT/r13**3 * (yT-yL) + mS/r23**3 * (yS-yL))
    ]

t = np.linspace(0, 20, 1000)

sol = odeint(dSdt, y0=[xT_0, yT_0, xS_0, yS_0, xL_0, yL_0,
                       vxT_0, vyT_0, vxS_0, vyS_0, vxL_0, vyL_0], t=t)

t = sol.T
x1 = sol.T[0]
y1 = sol.T[1]
x2 = sol.T[2]
y2 = sol.T[3]
x3 = sol.T[4]
y3 = sol.T[5]

# print(x1)
# print(y1)
plt.plot(y1)
plt.show()

tt = 1/np.sqrt(6.67e-11 * 1.99e30 / (1.5e11)**3 ) # seconds
tt = tt / (60*60 * 24* 365.25) * np.diff(t)[0] # per time step (in years)

def animate(i):
    ln1.set_data([x1[i], x2[i], x3[i]], [y1[i], y2[i], y3[i]])
    #text.set_text('Time = {:.1f} Years'.format(i*Tt))
fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.grid()
ln1, = plt.plot([], [], 'ro', lw=3, markersize=6)
#text = plt.text(0, 1.75, 'asdasd', fontsize=20, backgroundcolor='white', ha='center')
ax.set_ylim(-(yT_0 + 10**8), (yT_0 + 10**8))
ax.set_xlim(-(yT_0 + 10**8), (yT_0 + 10**8))
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
plt.show()
ani.save('terraluasol.gif',writer='pillow',fps=30)

