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

# fit parameters from 
eta = 1.66 # unit [1/d]
rho = 1/2.2 # unit [1/d]

# initial values
N = 763 # number of pupils
I_0 = 0.001 / N # initially infected
S_0 = 1 - I_0 # iinitially susceptible
R_0 = 0
t_0 = 0 # intitial time

# initial y-vector
y0 = np.array([S_0,I_0,R_0])

R0 = np.append(0,np.linspace(1,5,30))

print(R0)

R_inf = np.zeros_like(R0)

for ir0, r0 in enumerate(R0):

    # initialize integrator
    r = ode(dxdt)
    
    # Runge-Kutta with step size control
    r.set_integrator('dopri5')
    
    # set initial values
    r.set_initial_value(y0,t_0)

    # set transmission rate and recovery rate
    eta = r0*rho
    print(eta)
    r.set_f_params(eta,rho)
    
    # max value of time and points in time to integrate to
    t_max = 10000
    N_spacing_in_t = 3
    
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

    R_inf[ir0] = result[2,-1]


# plot result
color = ['#777777','#A95648','#88885A', ]

fig = pl.figure(figsize= (4,3))
pl.xlabel(r'number $R_0=\eta/\rho$ of secondary infections')
pl.ylabel('number of recovered')

fig.tight_layout()

pl.plot(R0,R_inf,color = 'k',lw=2)
#pl.plot(R0[R0>1],1-1.0/R0[R0>1])

fig.savefig('R0.pdf')

pl.show()
