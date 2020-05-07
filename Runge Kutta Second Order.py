import numpy as np
#Author: Jailson39
"Initially, these algorithm implements a integrator for the problem of the"
"Harmonic Oscilator. In addition, it can be updated to any other problem that"
"involves de calculation of a numerical approximation for the analytical "
"solution of a first order ODE."

#RK2:

m = 1
k = 8
xt0 = 2
yt0 = 2
vxt0 = 1
vyt0 = 1
axt0 = -k/m*xt0
ayt0 = -k/m*yt0
passos = 950
tempomax = 10
deltat = tempomax / passos
vx = np.ndarray(passos + 1)
vy = np.ndarray(passos + 1)
x = np.ndarray(passos + 1)
y = np.ndarray(passos + 1)
ax = np.ndarray(passos + 1)
ay = np.ndarray(passos + 1)
vx[0] = vxt0
vy[0] = vyt0
x[0]  = xt0
y[0]  = yt0
ax[0] = axt0
ay[0] = ayt0

for i in range(0,passos):
    k1Xx = vx[i]*deltat                 #papel de k1x pra x
    k1Xy = vy[i]*deltat                 #papel de k1x pra y

    k1vx = -(k/m)*x[i]*deltat           #x(t+dt/2) = x(t) + k1x/2   #papel de k1v pra x
    k1vy = -(k/m)*y[i]*deltat           #v(t+dt/2) = v(t) + k1v/2   #papel de k1v pra y

    k2Xx = (vx[i] + k1vx/2)*deltat       #k2x = v(t+dt/2)dt          #papel de k2x pra x
    k2Xy = (vy[i] + k1vy/2)*deltat       #k2v = -(k/m)*x(t+dt/2)dt   #papel de k2x pra y

    k2vx = -(k/m)*(x[i] + k1Xx/2)*deltat
    k2vy = -(k/m)*(y[i] + k1Xy/2)*deltat

    ax[i+1] = -k*x[i]/m
    ay[i+1] = -k*y[i]/m
    vx[i+1] = vx[i] + k2vx
    vy[i+1] = vy[i] + k2vy
    x[i+1] = x[i] + k2Xx
    y[i+1] = y[i] + k2Xy
