#######################
#Jose Barreiros, 2016 #
#Complex Systems      #
#Cornell University   #
#######################

#plot return map with Rugge Kutta
import math, operator, random
import numpy as np
import matplotlib.pyplot as plt

dataX0 = []
dataY0 = []
dataZ0 = []
yList = []
xList = []
zList = []
y_v=[]
z_v=[]

count = 1
t = 1
##
x = 0
y = 0
z = 0

xList.append(x)
yList.append(y)
zList.append(z)

h = .01

a = .398
b = 2
c = 4

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

X0 = []
Y0 = []
Z0 = []
previousLocalMax = x
localMax = x
t = 1
changeInTime = h
avgChange = []
LastPoint = changeInTime
while changeInTime < 20000 and len(dataX0) < 1000:

	[x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1])

	xList.append(x)
	yList.append(y)
	zList.append(z)
	if 201 < changeInTime:
                # Print y and z points when x is between 0.3 and 0.4
		if x <= 0.4 and x >= 0.3:
			print(y , "," , z)
			y_v.append(y)
			z_v.append(z)

		if x < xList[t-1] and xList[t-2] < xList[t-1]:
			previousLocalMax = localMax
			localMax = xList[t-1]
			dataX0.append(localMax)
			dataY0.append(previousLocalMax)
			# Calculate the change between this and last point
			avgChange.append(changeInTime - LastPoint)
			LastPoint = changeInTime

	t = t + 1
	changeInTime += h
#plt.plot(dataY0, dataX0, 'bo', y_v, z_v, 'mo')


ax=plt.figure(facecolor='white')
plt.plot(dataY0, dataX0, 'b.', markersize=2)
ax.suptitle('Return map - Rossler attractor\n a=%s b=%s c=%s'%(str(a),str(b),str(c)),fontsize=18)
plt.xlabel('X[n]', fontsize=14)
plt.ylabel('X[n+1]', fontsize=14)
plt.show()

#plt.plot(y_v, z_v, 'mo')
#plt.show()



