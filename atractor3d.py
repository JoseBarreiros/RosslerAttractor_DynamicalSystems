#######################
#Jose Barreiros, 2016 #
#Complex Systems      #
#Cornell University   #
#######################

#plot atractor 3D and 2D with Rugge Kutta
import math, random
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D



# Step Size
h = .1
#Initial conditions for Rossler system
x = 0
y = 0
z = 0

# Parameters for the Rossler System
a = .41 #.3
b = 2   #2
c = 4   #4



# Functions that define the system
def f(x,y,z):
	global a,b,c
	dxdt = -y-z
	return dxdt

def g(x,y,z):
	global a,b,c
	dydt = x + a * y
	return dydt

def e(x,y,z):
	global a,b,c
	dzdt = b + z * (x - c)
	return dzdt

# randomly perturb the initial conditions to create variable time series
x = x + random.random() / 2.0
y = y + random.random() / 2.0
z = z + random.random() / 2.0

dataX0 = []
dataY0 = []
dataZ0 = []
yList = []
xList = []
zList = []
tList = []
lamdaList = []
lyapunovList = []

t = 1

xList.append(x)
yList.append(y)
zList.append(z)
tList.append(z)
# Use the 4th order Runge-Kutta method
def rk4o(x, y, z):
	global h
	k1x = h*f(x, y, z)
	k1y = h*g(x, y, z)
	k1z = h*e(x, y, z)

	k2x = h*f(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)
	k2y = h*g(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)
	k2z = h*e(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)

	k3x = h*f(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)
	k3y = h*g(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)
	k3z = h*e(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)

	k4x = h*f(x + k3x, y + k3y, z + k3z)
	k4y = h*g(x + k3x, y + k3y, z + k3z)
	k4z = h*e(x + k3x, y + k3y, z + k3z)

	x = x + k1x/6.0 + k2x/3.0 + k3x/3.0 + k4x/6.0
	y = y + k1y/6.0 + k2y/3.0 + k3y/3.0 + k4y/6.0
	z = z + k1z/6.0 + k2z/3.0 + k3z/3.0 + k4z/6.0

	return [x,y,z]

t = 1
changeInTime = h
startLE = True
perturb=0
while changeInTime < 2000: # Perform 20000 / h iterations

	[x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1])

	xList.append(x)
	yList.append(y)
	zList.append(z)
	tList.append(t) 
	if 200 < changeInTime: # Remove the transient after 200 / h iterations
		if startLE:
			cx = xList[t-1] + perturb
			cy = yList[t-1]
			cz = zList[t-1]
			startLE = False

	t = t + 1
	#tList.append(t) 
	#t = t + 1
	changeInTime += h


print(len(tList))
print(len(xList))
#plt.plot(xList, yList, '-', linewidth=0.1)  #2D
#f, axarr = plt.subplots(2, sharex=True)

fig = plt.figure(figsize=plt.figaspect(.4),facecolor='white')
fig.suptitle('Rossler attractor\n a=%s'%(str(a)),fontsize=18)
ax = fig.add_subplot(1,2,1)
ax.plot(xList, yList, '-', linewidth=0.1)
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)

ax = fig.add_subplot(1,2,2,projection='3d')
ax.plot(xList, yList, zList, '-', linewidth=0.1) 
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)


plt.show()


