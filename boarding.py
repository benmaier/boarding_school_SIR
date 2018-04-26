import numpy as np
from scipy.integrate import ode

import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica Neue'

import matplotlib.pyplot as pl

# set equation of motion for SIR dynaics
def dxdt(t,y,eta,rho):

    S = y[0]
    I = y[1]
    R = y[2]

    dy = np.zeros(3)
    dy[0] = -eta*S*I
    dy[1] = +eta*S*I - rho*I
    dy[2] = +rho*I

    return dy

# load data
data = np.loadtxt('british_boarding_school.txt')
t_b = data[:,0] - 1
I_b = data[:,1]

# fit parameters from 
eta = 1.66 # unit [1/d]
rho = 1/2.2 # unit [1/d]

# initial values
N = 763 # number of pupils
I_0 = 3.0 / N # initially infected
S_0 = 1 - I_0 # iinitially susceptible
R_0 = 0
t_0 = 0 # intitial time

# initial y-vector
y0 = np.array([S_0,I_0,R_0])

# initialize integrator
r = ode(dxdt)

# Runge-Kutta with step size control
r.set_integrator('dopri5')

# set initial values
r.set_initial_value(y0,t_0)

# set transmission rate and recovery rate
r.set_f_params(eta,rho)

# max value of time and points in time to integrate to
t_max = t_b.max()
N_spacing_in_t = 100

# create vector of time points you want to evaluate
t = np.linspace(t_0,t_max,N_spacing_in_t)

# create vector of positions for those times
result = np.zeros((3,len(t)))


# loop through all demanded time points
for it, t_ in enumerate(t):

    # get result of ODE integration
    y = r.integrate(t_)

    # write result to result vector
    result[:,it] = y

# plot result
color = ['#777777','#A95648','#88885A', ]

fig = pl.figure(figsize= (4,3))
pl.xlabel('time $t$ [days]')
pl.ylabel('number of individuals')

pl.plot(t_b,I_b,marker='o',ls='None',ms= 7,zorder = 100,mfc='w',mec='k')

# save first to make mpl adjust the figure, then save again
fig.tight_layout()
fig.savefig('boarding_school_{0}.pdf'.format(0))


pl.ylim([-10,N+10])
pl.xlim([-0.5,t_max+0.5])
fig.tight_layout()
fig.savefig('boarding_school_{0}.pdf'.format(0))

for iplot, EPI in enumerate([1,0,2]):
    pl.plot(t,result[EPI,:]*N,'-',lw=3,c=color[EPI])
    pl.ylim([-10,N+10])
    pl.xlim([-0.5,t_max+0.5])
    fig.tight_layout()
    fig.savefig('boarding_school_{0}.pdf'.format(iplot+1))

pl.show()
