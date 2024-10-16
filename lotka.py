# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 17:07:56 2024

@author: Josh
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import animation

y0 = [10, 2]
alt_y0 = [5, 3]
alt2_y0 = [13, 1]
alt3_y0 = [4, 2.75]

t = np.linspace(0,50, num = 1000)

a = 1.1
b = 0.4
c = 0.4
d = 0.1

parameters = [a, b, c, d]
def yDot(variables, t, parameters):
    
    x, y = variables
    a, b, c, d = parameters
    
    dxdt = a*x - b*x*y
    dydt = -c*y + d*x*y
    
    return([dxdt, dydt])    


y = odeint(yDot, y0, t, args = (parameters, ) )
alt_y = odeint(yDot, alt_y0, t, args = (parameters, ) )
alt2_y = odeint(yDot, alt2_y0, t, args = (parameters, ) )
alt3_y = odeint(yDot, alt3_y0, t, args = (parameters, ) )

y_slice0 = y[:,0]
y_slice1 = y[:,1]

alt_y_slice0 = alt_y[:,0]
alt_y_slice1 = alt_y[:,1]

alt2_y_slice0 = alt2_y[:,0]
alt2_y_slice1 = alt2_y[:,1]

alt3_y_slice0 = alt3_y[:,0]
alt3_y_slice1 = alt3_y[:,1]



##################################################################################

x = np.linspace(0, 20, 20)
y = np.linspace(0, 20, 20)
X, Y = np.meshgrid(x, y)

U, V = yDot([X,Y], 0, parameters)

N = np.sqrt(U**2 + V**2)
U2, V2 = U/N, V/N

fig, ax = plt.subplots()


ax.quiver(X, Y, U2, V2, N, cmap=plt.cm.plasma)

alt3_traj, = ax.plot([], [], c = 'm', label = '[4,2.75]')
alt_traj, = ax.plot([], [], c = 'g', label = '[5,3]')
traj, = ax.plot([], [], c = 'r', label = '[10,2]')
alt2_traj, = ax.plot([], [], c = 'b', label = '[13,1]')

current_traj, = ax.plot([], [], 'ro')
alt_current_traj, = ax.plot([], [], 'go')
alt2_current_traj, = ax.plot([], [], 'bo')
alt3_current_traj, = ax.plot([], [], 'mo')
ax.set(xlim=[-1, 20], ylim=[-1, 20], xlabel='Prey', ylabel='Predators')
ax.legend(title = '[Prey, Predator] at t=0')



def updateTraj(frame):
    traj.set_data(y_slice0[0:frame], y_slice1[0:frame])
    current_traj.set_data((y_slice0[frame - 1], ), (y_slice1[frame - 1], ))
    
    alt_traj.set_data(alt_y_slice0[0:frame], alt_y_slice1[0:frame])
    alt_current_traj.set_data((alt_y_slice0[frame - 1], ), (alt_y_slice1[frame - 1], ))
    
    alt2_traj.set_data(alt2_y_slice0[0:frame], alt2_y_slice1[0:frame])
    alt2_current_traj.set_data((alt2_y_slice0[frame - 1], ), (alt2_y_slice1[frame - 1], ))
    
    alt3_traj.set_data(alt3_y_slice0[0:frame], alt3_y_slice1[0:frame])
    alt3_current_traj.set_data((alt3_y_slice0[frame - 1], ), (alt3_y_slice1[frame - 1], ))
    
    
    return(traj, current_traj, alt_traj, alt_current_traj, alt2_traj, alt2_current_traj, alt3_traj, alt3_current_traj)

ani2 = animation.FuncAnimation(fig, func=updateTraj, frames=len(y_slice0), interval=20, blit=True)
plt.show()

###################################################################################

fig, ax = plt.subplots()
line1, = ax.plot([], [], c='r', linestyle='dotted', label='Prey')  # Unpack the line object correctly
line2, = ax.plot([], [], c='r', label='Predator')
ax.set(xlim=[0, 50], ylim=[0, 15], xlabel='Time', ylabel='Population')
ax.legend()


def update(frame):
    # For each frame, update the data stored on each line.
    line1.set_data(t[:frame], y_slice0[:frame])
    line2.set_data(t[:frame], y_slice1[:frame])
    return (line1, line2)

# Create animation
ani = animation.FuncAnimation(fig, func=update, frames=len(t), interval=20, blit=True)
plt.show()

#################################################################################
