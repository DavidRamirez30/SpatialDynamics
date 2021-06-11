import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import toeplitz       

### Model equations
def GoL(y, t, b, beta, gamma, dC, n1, n2, F, G):
    dy = y * (1 - y) * (1 + beta * ( ((F**n1)*(y**n1)) / (1 + ((F**n1)*(y**n1))) )) * ( (1 + ((b**n2)*(G**n2)*(y**n2)) ) / (1 + gamma * ((b**n2)*(G**n2)*(y**n2)) )) - dC*y;
    return dy, t
    
### Create interaction matrices
def formMatrix(n, a, q, d):
    
    # define cooperation (positive interaction) input matrix
    F = np.zeros((n,n))
    
    # define competition (negative interaction) input matrix
    G = np.zeros((n,n))    
    
    # assessing the elements of these matrices
    for i in range(n):
        
        # ith row of F defines the positive input to the ith cluster
        # from others
        F[i,:] = a[i,:]/d[:,i]
        
        # ith row of G defines the negative input to the ith cluster
        # from others
        G[i,:] = q[i,:]/d[:,i]
        
        # interaction of every cluster with itself (diagonal terms)
        # are set to 0 (zero)
        F[i,i] = 0
        G[i,i] = 0
        
    return F, G

### Create model

# enter part of number if the job was splitted, default value is 1
Part = 1

# total initial number of clusters
n = 16

# define spatial configuration
# square
lis1 = []
lis2 = []
for i in range(int(np.sqrt(n))):
    for k in range(int(np.sqrt(n))):
        lis1.append(i)
        lis2.append(k)
        
pos = np.array([lis1, lis2])



# calculate distances between each clusters
# it generates a n-by-n matrix whose ith column represents distances 
# to the ith cluster from neighbors
dst=[]
for x in range(len(pos[0])):
    p1 = np.array((pos[0][x],pos[1][x]))
    for y in range(len(pos[1])):
        p2 = np.array((pos[0][y],pos[1][y]))
        val = np.linalg.norm(p1-p2)
        dst.append(val)

dst = np.array([dst])
dst = np.reshape(dst, (16, 16))

# scale the nearest neighbor distance
dmin = 20

# for example, try for different values of dmin 
# (nearest neighbor distance)
dminvec = np.array([dmin])

# define kinetic model parameters
b = 14.1; # ratio of competition/cooperation length scales 
beta = 100; # coefficient of max. promotion due to cooperation
gamma = 100; # coefficient of max. suppression due to competition
dC = 0.2; # death rate
n1 = 3; # Hill coefficient for cooperation (positive interaction) function
n2 = 3; # Hill coefficient for competition (negative interaction) function

# This constant scales the maximum initial value from there 
C0m = 1e-3;

# define simulation time interval
tend = 200
tspan = np.linspace(0, tend)


# positive (cooperative) interactions
# intracluster interaction coefficient
# (multiplies the spatial distances, if desired)
ma = 2
sa = 1.0*ma

# here also the symmetry in the interactions can be set
# a = toeplitz(np.random.normal(ma, sa, (n,1)))
a = np.random.normal(ma, sa, (n, n))

# negative (competitive) interactions
mq = 2
sq = 1.0*mq
# q = toeplitz(np.random.normal(mq, sq, (n,1)))
q = np.random.normal(mq, sq, (n, n))

xx = 0

# define the initial condition 
# a randomly generated initial condition of
# occupied or unoccupied for each vertex
mid = np.random.randint(2, size=16)
C = mid.tolist()
C.insert(0, 1)
C = np.array([C])
C = C.reshape(17,1)

m = (C.shape)[1]

# create an outout matrix
output = np.zeros((m*len(dminvec),1*n+2))

# non-negative ODE?

for zz in range(m):
    y0 = C[1:,]
    y0 = C0m * y0
    tv = y0.tolist()
    
    yyy = .001, .001, .001, 0, 0, 0, .001, .001, .001, .001, 0, .001, 0, .001, .001, .001
    
    for bb in range(len(dminvec)):
        b = b
        beta = beta
        gamma = gamma
        dC = dC
        dmin = dminvec
        d = dmin*dst
        
        F, G = formMatrix(n, a, q, d)
        
        output[xx][0] = C[0][zz]
        output[xx][1] = dmin
        
        soln = odeint(GoL, yyy, tspan, args=(b, beta, gamma, dC, n1, n2, F, G,))

fig1 = plt.figure(num=1, clear=True)
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(pos[0], pos[1], "ko")
ax1.set(xlabel=r'x position ($\times d_{min}$ = ' + str(dmin) + ')',
        ylabel=r'y position ($\times d_{min}$ = ' + str(dmin) + ')',
        title='Square lattice')
ax1.grid(True)
ax1.axis('equal')
fig1.tight_layout()









