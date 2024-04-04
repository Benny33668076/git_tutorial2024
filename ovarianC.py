import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the parameters and initial conditions
gamma = 1.25; # growth rate of proliferating cells
alpha = 0.6685; # transition rate from Q to P
beta = 0.05; # transition rate from P to Q
mu = 0.477; # transition rate from Q to D
phi = 0.01; # decay rate of dead cells
K = 1000; # carrying capacity

P0 = 0.18; # initial number of proliferating cells
Q0 = 0; # initial number of quiescent cells
D0 = 0; # initial number of dead cells


# Define the system of ODEs 
def system(y, t):
    P, Q, D = y

    dPdt = gamma * P * (1 - P / K) - beta * P + alpha * Q  
    dQdt = beta * P - (alpha + mu) * Q  
    dDdt = mu * Q - phi * D  

    return [dPdt, dQdt, dDdt]

# Time span (in days)
tspan = np.linspace(0,100)

# Initial conditions vector 
y0 =[P0,Q0,D0]

# Solve ODEs using scipy's odeint 
y=odeint(system,y0,tspan)

# Plotting the solution 
plt.plot(tspan,y[:,0],'r-',label='P')
plt.plot(tspan,y[:,1],'g-',label='Q')
plt.plot(tspan,y[:,2],'b-',label='D')

plt.xlabel('Time (days)')
plt.ylabel('P, Q and D')
plt.legend(loc='best')
plt.title('Dynamics of ovarian cancer')
plt.show()


# Plotting P(t)
plt.figure(figsize=(8, 5))
plt.plot(tspan, y[:,0], 'r-', label='P')
plt.xlabel('Time (days)')
plt.ylabel('Number of P Cells')
plt.title('Dynamics of P Cells in Ovarian Cancer')
plt.legend(loc='best')
plt.show()

# Plotting Q(t)
plt.figure(figsize=(8, 5))
plt.plot(tspan, y[:,1], 'g-', label='Q')
plt.xlabel('Time (days)')
plt.ylabel('Number of Q Cells')
plt.title('Dynamics of Q Cells in Ovarian Cancer')
plt.legend(loc='best')
plt.show()

# Plotting D(t)
plt.figure(figsize=(8, 5))
plt.plot(tspan, y[:,2], 'b-', label='D')
plt.xlabel('Time (days)')
plt.ylabel('Number of D Cells')
plt.title('Dynamics')
plt.title('Dynamics of Q Cells in Ovarian Cancer')
plt.legend(loc='best')
plt.show()