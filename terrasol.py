import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.constants import gravitational_constant as G
from matplotlib import animation
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter

Tt = 365 * 360 * 24

m1 = 5.97e24
m2 = 1.98e30
x1_0 = 149.6e6
y1_0 = 0
x2_0 = 0
y2_0 = 0
vx1_0 = 0
vy1_0 = 298000
vx2_0 = 0
vy2_0 = 0

def dSdt(S, t):
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = S
    r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    return [ vx1,
            vy1,
            vx2,
            vy2,
            G * (m2/r12**3 * (x2-x1)),
            G * (m2/r12**3 * (y2-y1)),
            G * (m1/r12**3 * (x1-x2)),
            G * (m1/r12**3 * (y1-y2))
        ]

t = np.linspace(0,1,10000)

sol = odeint(dSdt, y0=[x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0],
             t=t)

x1 = sol.T[0]
y1 = sol.T[1]
x2 = sol.T[2]
y2 = sol.T[3]

print(x1)
print("--------------------------")
print(y1)

def animate(i):
    ln1.set_data([x1[i], x2[i]], [y1[i], y2[i]])
    # text.set_text('Time = {:.2f} Years'.format(i*tt))
    
fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.grid()
ln1, = plt.plot([], [], 'ro--', lw=3, markersize=8)
# text = plt.text(0.7, 0.7, '')
ax.set_ylim(-(y1_0 + 10**8), (y1_0 + 10**8))
ax.set_xlim(-(y1_0 + 10**8), (y1_0 + 10**8))
ani = animation.FuncAnimation(fig, animate, frames=200, interval=50)
plt.show()
# ani.save('plan.gif',writer='pillow',fps=30)
