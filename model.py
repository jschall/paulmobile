import numpy as np
from math import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def sqrt_ctrl(err, p, amax):
    lin_dist = amax/p**2

    if err > lin_dist:
        return sqrt(2.*amax*(err-lin_dist/2.))
    elif err < -lin_dist:
        return -sqrt(2*amax*(-err-lin_dist/2.))
    else:
        return err*p

def dyn(y, t):
    g = 9.8065
    Mw = .1
    Iw = Mw/2.
    r = 0.1
    Mp = 2.5
    Ip = 2.*0.4**2+.5*(2-0.4)**2
    l = 0.4

    beta = 2*Mw+2*Iw+Mp
    alpha = Ip*beta+2*Mp*l**2*(Mw+Iw/r**2)

    A23 = Mp**2*g*e**2/alpha
    A43 = Mp*g*l*beta/alpha
    B2 = (Ip+Mp*l**2-Mp*l*r)/(r*alpha)
    B4 = (Mp*l-r*beta)/(r*alpha)

    A = np.matrix([
        [0,1,0,0],
        [0,0,A23,0],
        [0,0,0,1],
        [0,0,A43,0]])


    B = np.matrix([
        [0],
        [B2],
        [0],
        [B4]])

    desPos = 5.*sin(t*0.4)
    desVel = sqrt_ctrl(desPos-y[0], .275, radians(10)*A23)# + 2.*cos(t*0.4)
    if desVel > 10.:
        desVel=10.
    elif desVel < -10.:
        desVel=-10.
    #desVel = 10.
    desAcc = (desVel-y[1])*.9# + -0.8*sin(t*0.4)
    if desAcc > radians(10)*A23:
        desAcc = radians(10)*A23
    elif desAcc < -radians(10)*A23:
        desAcc = -radians(10)*A23
    #desAcc = 0.
    desPhi = desAcc/A23
    desPhiDot = (desPhi-y[2])*3.
    torque = (desPhiDot-y[3])*20.
    if torque > 1.5:
        torque = 1.5
    elif torque < -1.5:
        torque = -1.5

    return np.squeeze(np.asarray(A*np.matrix(y).T+B*torque))

state = np.array([0.,1.,0.0,0.])

from visual import *
t = 0.
dt = 1./60.
body = cylinder(pos=(0,0.1,-0.2),axis=(0,0,.4),radius=.1, color=color.orange)
pendulum = cylinder(pos=(0,0,0), axis=(0,2,0), radius=.02, color=color.orange)
despos = sphere(pos = (0,0,0), radius=.05)
cylinder(pos=(0,0,0), axis=(0,-.1,0),radius=10., material=materials.wood)
while True:
    rate(60)
    state += dyn(state, t)*dt
    body.pos = (state[0],0.1,-0.2)
    pendulum.pos = (state[0],0.1,0)
    pendulum.axis = (2.*sin(state[2]),2.*cos(state[2]),0)
    despos.pos = (5.*sin(t*0.4),1,0)
    t += dt

x = np.linspace(0,20,1000)

y = odeint(dyn,(0.,0.,radians(0),0.),x)

l1 = plt.plot(x,y[0:,0])
l2 = plt.plot(x,y[0:,1])
l3 = plt.plot(x,y[0:,2]*degrees(1))
l4 = plt.plot(x,y[0:,3]*degrees(1))
plt.legend(["$x$", "$\dot x$", "$\phi$", "$\dot \phi$"])
plt.show()
