#f1=gcurve()
#fT=gcurve(color=color.red)
#fU=gcurve(color=color.blue)
from vpython import * 
L1=1
L2=1

pivot=sphere(pos=vector(0,L1,0), radius=0.05)

m1=sphere(pos=pivot.pos-vector(0,L1,0),radius=0.05,color=color.red)
stick1=cylinder(pos=pivot.pos,axis=m1.pos-pivot.pos,radius=0.015,color=color.yellow)

m2=sphere(pos=m1.pos-vector(0,L2,0),radius=0.05,color=color.red)
stick2=cylinder(pos=m1.pos, axis=m2.pos-m1.pos, radius=0.015,color=color.yellow)

M1=1
M2=1

g=9.8
theta1=75*pi/180
theta2=180*pi/180
theta1dot=0
theta2dot=0.8

m1.pos=pivot.pos+vector(L1*sin(theta1),-L1*cos(theta1),0)
m2.pos=m1.pos+vector(L2*sin(theta2),-L2*cos(theta2),0)
stick1.axis=m1.pos-pivot.pos
stick2.pos=m1.pos
stick2.axis=m2.pos-m1.pos

a=(M1+M2)*L1**2
b=M2*L1*L2
c=-(M1+M2)*g*L1*theta1
d=M2*L1*L2
k=M2*L2**2
f=-M2*g*L2*theta2

t=0
dt=0.0001
attach_trail(m2, retain=200)
while t<20:
    rate(10000)
    a=(M1+M2)*L1**2
    b=M2*L1*L2
    c=-(M1+M2)*g*L1*theta1
    d=M2*L1*L2
    k=M2*L2**2
    f=-M2*g*L2*theta2
    
   # theta2ddot=(c/b-(a*f)/(b*d))/(1-a*k/(b*d))
    #theta1ddot=(f-k*theta2ddot)/d
    theta1ddot=(-g*(2*M1+M2)*sin(theta1)-M2*g*sin(theta1-2*theta2)-2*sin(theta1-theta2)*M2*(L2*theta2dot**2+L1*cos(theta1-theta2)*theta1dot**2))/(L1*(2*M1+M2-M2*cos(2*theta1-2*theta2)))
    theta2ddot=(2*sin(theta1-theta2)*((M1+M2)*L1*theta1dot**2+g*(M1+M2)*cos(theta1)+L2*M2*cos(theta1-theta2)*theta2dot**2))/(L2*(2*M1+M2-M2*cos(2*theta1-2*theta2)))
    theta2dot=theta2dot+theta2ddot*dt
    theta1dot=theta1dot+theta1ddot*dt
    theta1=theta1+theta1dot*dt
    theta2=theta2+theta2dot*dt
    m1.pos=pivot.pos+vector(L1*sin(theta1),-L1*cos(theta1),0)
    m2.pos=m1.pos+vector(L2*sin(theta2),-L2*cos(theta2),0)
    stick1.axis=m1.pos-pivot.pos
    stick2.pos=m1.pos
    stick2.axis=m2.pos-m1.pos
    t=t+dt
   # T=.5*M1*L1**2*theta1dot**2+.5*M2*(L1**2*theta1**2+2*L1*L2*cos(theta1)*cos(theta2)*theta1dot*theta2dot+L2**2*theta2dot**2+2*L1*L2*sin(theta1)*sin(theta2)*theta1dot*theta2dot)
    #U=-M1*L1*g*cos(theta1)-M2*g*(L1*cos(theta1)+L2*cos(theta2))
    x1=L1*sin(theta1)
    y1=-L1*cos(theta1)
    x1dot=L1*theta1dot*cos(theta1)
    y1dot=L1*theta1dot*sin(theta1)
    x2=x1+L2*sin(theta2)
    y2=y1-L2*cos(theta2)
    x2dot=L1*theta1dot*cos(theta1)+L2*theta2dot*cos(theta2)
    y2dot=L1*theta1dot*sin(theta1)+L2*theta2dot*sin(theta2)
    T=.5*M1*(x1dot**2+y1dot**2)+.5*M2*(x2dot**2+y2dot**2)
    U=M1*g*y1+M2*g*y2
    
    
    E=T+U
#   f1.plot(t,E)
#   fT.plot(t,T)
#   fU.plot(t,U)