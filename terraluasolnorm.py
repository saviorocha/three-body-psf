import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import gravitational_constant as G
from scipy.integrate import solve_ivp
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter

mT = 1
mS = 1
mL = 1
v1 = 0.39295
v2 = 0.09758
# mS = 333000
# mL = 0.012328308

xT_0 = 0
yT_0 = 1
xS_0 = 0
yS_0 = 0
xL_0 = 0
yL_0 = yT_0 + 0.1

vxT_0 = v1
vyT_0 = 0
vxS_0 = 0
vyS_0 = 0
vxL_0 = vxT_0 + 0.1
# vxL_0 = vxT_0 + (2 * math.pi * yL_0) / Tt * 10**2
vyL_0 = 0

def dSdt(t, S):
    xT, yT, xS, yS, xL, yL, vxT, vyT, vxS, vyS, vxL, vyL = S

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
        (mS/r12**3 * (xS-xT) + mL/r13**3 * (xL-xT)),
        (mS/r12**3 * (yS-yT) + mL/r13**3 * (yL-yT)),
        (mT/r12**3 * (xT-xS) + mL/r23**3 * (xL-xS)),
        (mT/r12**3 * (yT-yS) + mL/r23**3 * (yL-yS)),
        (mT/r13**3 * (xT-xL) + mS/r23**3 * (xS-xL)),
        (mT/r13**3 * (yT-yL) + mS/r23**3 * (yS-yL))
    ]

t = np.linspace(0, 20, 1000)


sol = solve_ivp(dSdt, (0,20), y0=[xT_0, yT_0, xS_0, yS_0, xL_0, yL_0,
                       vxT_0, vyT_0, vxS_0, vyS_0, vxL_0, vyL_0],
                method = 'DOP853', t_eval=t, rtol=1e-10, atol=1e-13)


t = sol.t
x1 = sol.y[0]
y1 = sol.y[1]
x2 = sol.y[2]
y2 = sol.y[3]
x3 = sol.y[4]
y3 = sol.y[5]

# print(x1)
# print(y1)

tt = 1/np.sqrt(6.67e-11 * 1.99e30 / (1.5e11)**3 ) # seconds
tt = tt / (60*60 * 24* 365.25) * np.diff(t)[0] # per time step (in years)

def animate(i):
    ln1.set_data([x1[i], x2[i], x3[i]], [y1[i], y2[i], y3[i]])
    # text.set_text('Time = {:.1f} Years'.format(i*Tt))
fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.grid()
ln1, = plt.plot([], [], 'ro', lw=3, markersize=6)
# text = plt.text(0, 1.75, 'asdasd', fontsize=20, backgroundcolor='white', ha='center')
ax.set_ylim(-5, 5)
ax.set_xlim(-5, 5)
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
plt.show()
ani.save('plan.gif',writer='pillow',fps=30)

