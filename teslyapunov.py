from mimetypes import init
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque

## On pose les constantes de notre problème#
g = 9.80
m1 = 1 
m2 = 1
l1 = 1
l2 = 1
tmax = 180
k =0.02 #0.0002
l = l1+l2
qi1 = np.pi/2
qi2 = np.pi/2
pi1 = 0
pi2 = 0
N = int(tmax/k)  #int(tmax/k)
fps = 1/k
history_len = 300  # how many trajectory points to display


ini = np.array([[qi1, qi2, pi1, pi2]])
ini_modif = np.array([[qi1+0.001, qi2+0.001, pi1+0.001,pi2+0.001]])


def M(q1, q2,p1,p2):
    
    dq1 = p1
    dq2 = p2
    dp1 = (-g*(2*m1+m2)*np.sin(q1)-m2*g*np.sin(q1-2*q2)-2*np.sin(q1-q2)*m2*(p2**2*l2+p1**2*l1*np.cos(q1-q2)))/(l1*(2*m1+m2-m2*np.cos(2*q1-2*q2)))
    dp2 = (2*np.sin(q1-q2)*(p1**2*l1*(m1+m2)+g*(m1+m2)*np.cos(q1)+p2**2*l2*m2*np.cos(q1-q2)))/(l2*(2*m1+m2-m2*np.cos(2*(q1-q2))))
    return dq1,dq2,dp1,dp2



def RK4(q1, q2, p1, p2):

    ancien = np.array([q1, q2, p1, p2])
    k1 = np.array(M(q1, q2, p1, p2))
    k2 = np.array(M(q1 + k*k1[0]/2, q2 + k*k1[1]/2, p1 + k*k1[2]/2, p2 + k*k1[3]/2))
    k3 = np.array(M(q1 + k*k2[0]/2, q2 + k*k2[1]/2, p1 + k*k2[2]/2, p2 + k*k2[3]/2))
    k4 = np.array(M(q1 + k*k3[0], q2 + k*k3[1], p1 + k*k3[2], p2 + k*k3[3]))

    new = ancien + (k/6) * (k1 + 2*k2 + 2*k3 + k4)
    return new

sol = np.empty((0,4))
sol = np.concatenate((sol, ini))
sol_modified = np.empty((0,4))
sol_modified = np.concatenate((sol_modified, ini_modif))



#On boucle Runge-Kutta N fois  sur chacun des points
def calcul(sol):
    t = np.linspace(0, tmax, N+1)
    for i in range(int(N+1)):
        a = RK4(sol[i,0], sol[i,1], sol[i,2], sol[i,3])
        a = np.array([a])
        
        sol = np.concatenate((sol, a))

    return sol, t

# on calcul tout

solution, t = calcul(sol)

solution2, t = calcul(sol_modified)


#position des pendules #1
x1 = l1*np.sin(solution[:,0])
#print(x1.shape)
y1 = -l1*np.cos(solution[:,0])
x2 = x1 + l2*np.sin(solution[:,1])
y2 = y1 - l2*np.cos(solution[:,1])
y2s = 2 + y1 - l2*np.cos(solution[:,1])
#2
x11 = l1*np.sin(solution2[:,0])
y11 = -l1*np.cos(solution2[:,0])
x22 = x11 + l2*np.sin(solution2[:,1])
y22 = y11 - l2*np.cos(solution2[:,1])
y2ss = 2 + y11 - l2*np.cos(solution2[:,1])

### définition de l'énergie

Ec1, Ec2 = np.zeros(N+1), np.zeros(N+1)
Ep1, Ep2 = np.zeros(N+1), np.zeros(N+1)
Em1, Em2 = np.zeros(N+1), np.zeros(N+1)
Em = np.zeros(N+1)

def energie():
    Ec1[0] = (1./2)*m1*(l1**2)*(solution[0,0]**2)
    Ec2[0]=(1./2)*m2*(l1**2*solution[0,2]**2+l2**2*solution[0,3]**2+2*l1*l2*solution[0,2]*solution[0,3]*np.cos(solution[0,0]-
solution[0,1]))
   
    Ep1[0]=-m1*g*l1*np.cos(solution[0,0])
    Ep2[0]=-m2*g*(l1*np.cos(solution[0,0])+l2*np.cos(solution[0,1]))
    
    Em2[0]=Ec2[0]+Ep2[0]
    Em1[0]=Ec1[0]+Ep1[0]
    
    for i in range(N):
        Ec1[i+1]=(1./2)*m1*(l1**2)*(solution[i,0]**2)
        Ec2[i+1]=(1./2)*m2*(l1**2*solution[i,2]**2+l2**2*solution[i,3]**2+2*l1*l2*solution[i,2]*solution[i,3]*np.cos(solution[i,0]-
solution[i,1]))
        
 
        Ep1[i+1]=-m1*g*l1*np.cos(solution[i,1])
        Ep2[i+1]=-m2*g*(l1*np.cos(solution[i,1])+l2*np.cos(solution[i,2]))
        
        Em1[i+1]=Ec1[i+1]+Ep1[i+1]
        Em2[i+1]=Ec2[i+1]+Ep2[i+1]

    Em = Em1 + Em2
    return Em

####On crée le graphe pour l'animation 

pos1 = x1 + y1
pos2 = x11 + y11
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-l, l), ylim=(-l, 1.))
ax.set_aspect('equal')
ax.grid()
line, = ax.plot([], [], 'o-', lw=2)
line2, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], 'r-', lw=1, ms=0.2)
trace2, = ax.plot([], [], 'y-', lw=1, ms=0.2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)
history_xx, history_yy = deque(maxlen=history_len), deque(maxlen=history_len)

#On définit l'animation
#####################################################################################################################
######### activer pour avoir 2 pendules /!\/!\ PROBLEME LE 2EME PENDULE EST 2 FOIS PLUS LENT ET TREMBLE
#####################################################################################################################
def animate(i):

    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]
    thisxx = [0, x11[i], x22[i]]
    thisyy = [0, y11[i], y22[i]]

    if i == 0 : 
        history_x.clear()
        history_y.clear()
        history_xx.clear()
        history_yy.clear()
    history_x.appendleft(thisx[2])
    history_y.appendleft(thisy[2])
    history_xx.appendleft(thisxx[2])
    history_yy.appendleft(thisyy[2])
    line.set_data(thisx, thisy)
    line2.set_data(thisxx, thisyy)
    trace.set_data(history_x, history_y)
    trace2.set_data(history_xx, history_yy)
    time_text.set_text(time_template % (i*k))
    return line, time_text, trace, line2, trace2
def animate2(i):

    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]
   

    if i == 0 : 
        history_x.clear()
        history_y.clear()
       
    history_x.appendleft(thisx[2])
    history_y.appendleft(thisy[2])
   
    line.set_data(thisx, thisy)
    
    trace.set_data(history_x, history_y)
    
    time_text.set_text(time_template % (i*k))
    return line, time_text, trace

# ani = animation.FuncAnimation(
#     fig, animate, N-1,  interval = k*1000, blit = True) ## a priori interval = fps = 1/k
# ani.save('double_pendule.mp4')
# plt.show()
# plt.close()

################Calcul exposant de Lyapunov
def Jacobienne(q1,q2,p1,p2): ##matrice Jacobienne de notre système d'EDO
    deno1 = l1*(2*m1+m2-m2*np.cos(2*(q1-q2)))
    deno1q1 = l1*(2*m2*np.sin(2*(q1-q2)))
    deno1q2 = -l1*(2*m2*np.sin(2*(q1-q2)))
    J11, J12, J13, J14 = 0, 0, 1, 0
    J21, J22, J23, J24 = 0, 0, 0, 1
    J1 = np.array([J11, J12, J13, J14])
    J2 = np.array([J21, J22, J23, J24])

    w1q1 = (-g*(2*m1+m2)*np.cos(q1)-m2*g*np.cos(q1-2*q2)-2*np.cos(q1-q2)*m2*(p2**2*l2+p1**2*l1*np.cos(q1-q2))-2*np.sin(q1-q2)*m2*(p2**2*l2-p1**2*l1*np.sin(q1-q2)))
    w1 = (-g*(2*m1+m2)*np.sin(q1)-m2*g*np.sin(q1-2*q2)-2*np.sin(q1-q2)*m2*(p2**2*l2+p1**2*l1*np.cos(q1-q2)))/(l1*(2*m1+m2-m2*np.cos(2*q1-2*q2)))
    w1q2 = 2*m2*np.cos(q1-2*q2)+2*np.cos(q1-q2)*m2*(p2**2*l2+p1**2*l1*np.cos(q1-q2))-2*np.sin(q1-q2)*m2*(p2**2*l2+p1**2*l1*np.sin(q1-q2))
    J31, J32, J33, J34 = (w1q1*deno1-w1*deno1q1)/(deno1**2),(w1q2*deno1-w1*deno1q2)/(deno1**2),(-4*m2*l1*p1*np.sin(q1-q2)*np.cos(q1-q2))/deno1,(-4*m2*l2*p2*np.sin(q1-q2))

    J3 = np.array([J31, J32, J33, J34])
    if J3.shape != (4,):
        J3 = np.reshape(J3, (4))


    deno2 = l2*(2*m1+m2-m2*np.cos(2*(q1-q2)))
    deno2q1 = 2*l2*m2*np.sin(2*(q1-q2))
    deno2q2 = - deno2q1
    w2 = (2*np.sin(q1-q2)*(p1**2*l1*(m1+m2)+g*(m1+m2)*np.cos(q1)+p2**2*l2*m2*np.cos(q1-q2)))/(l2*(2*m1+m2-m2*np.cos(2*(q1-q2))))
    w2q1 = 2*np.cos(q1-q2)*(p1**2*l1*(m1+m2)+g*(m1+m2)*np.cos(q1)+p2**2*l2*m2*np.cos(q1-q2))+2*np.sin(q1-q2)*(-g*(m1+m2)*np.sin(q1)-p2**2*l2*m2*np.sin(q1-q2))
    w2q2 = -2*np.cos(q1-q2)*(p1**2*l1*(m1+m2)+g*(m1+m2)*np.cos(q1)+p2**2*l2*m2*np.cos(q1-q2))+2*np.sin(q1-q2)*(p2**2*l2*m2*np.sin(q1-q2))
    J41 = (w2q1*deno2-w2*deno2q1)/(deno2**2)
    J42 = ((w2q2*deno2-w2*deno2q2)/(deno2**2))
    J43 = (4*np.sin(q1-q2)*p1*l1*(m1+m2))/(deno2)
    J44 = (4*np.sin(q1-q2)*p2*l2*m2*np.cos(q1-q2))/(deno2)
    J4 = np.array([J41, J42, J43, J44])
    if J4.shape != (4,):
        J4 = np.reshape(J4, (4))
    J = np.zeros((4,4))
    J[0,:] = J1
    J[1,:] = J2
    J[2,:] = J3
    J[3,:] = J4
    return J
    ######test jacobienne
# print('jacobienne', Jacobienne(1,2,3,4))
# def coor_delta(q1,q2,p1,p2,q1m,q2m,p1m,p2m):
#     d1 = q1-q1m
#     d2 = q2-q2m
#     d3 = p1-p1m
#     d4 = p2-p2m
#     return d1, d2, d3, d4
# dq1, dq2, dp1, dp2 = coor_delta(1,2,3,4,1+1e-8,2+1e-8,3+1e-8,4+1e-8)
q1, q2, p1, p2 = 1, 2, 3, 4
q1m, q2m, p1m, p2m = q1+1e-8, q2+1e-8, p1+1e-8, p2+1e-8
d1, d2, d3, d4 = q1-q1m, q2-q2m, p1-p1m, p2-p2m
def delta(dq1,dq2,dp1,dp2):
    d = np.zeros((4,1))

    # #d[0,0], d[1,0], d[2,0], d[3,0] = coor_delta(q1,q2,p1,p2,q1m,q2m,p1m,p2m)####Deuxieme tenta
    d[0,0] = dq1##premier 
    d[1,0] = dq2##premier
    d[2,0] = dp1##premier
    d[3,0] = dp2##premier
    #print('delta',d)
    return d

######Test delta
# blbl = delta(d1,d2,d3,d4)
# print('delta',blbl)
######
def norme_delta(q1,q2,p1,p2,q1m,q2m,p1m,p2m): ###On veut definir delta comme etant l'ecart entre les deux solutions
    norme = np.sqrt((q1-q1m)**2+(q2-q2m)**2+(p1-p1m)**2+(p2-p2m)**2)
    return norme

def Jacobdelta(q1,q2,p1,p2,dq1,dq2,dp1,dp2):
    j = np.dot(Jacobienne(q1,q2,p1,p2), delta(dq1,dq2,dp1,dp2))

    return j
#Jacobdelta(q1,q2,p1,p2,d1,d2,d3,d4)
##########################
###On doit definir de nouvelles coordonnees ∂x = x'-x
# def coor_delta(q1,q2,p1,p2,q1m,q2m,p1m,p2m):
#     d1 = q1-q1m
#     d2 = q2-q2m
#     d3 = p1-p1m
#     d4 = p2-p2m
#     return d1, d2, d3, d4
# dq1, dq2, dp1, dp2 = coor_delta(1,2,3,4,1+1e-8,2+1e-8,3+1e-8,4+1e-8)
def RK(q1,q2,p1,p2, dq1, dq2, dp1, dp2):

    ancien = np.array([q1, q2, p1, p2])
    k1 = np.array(Jacobdelta(q1,q2,p1,p2,dq1,dq2,dp1,dp2))
    k2 = np.array(Jacobdelta(q1,q2,p1,p2,dq1 + k*k1[0]/2, dq2 + k*k1[1]/2, dp1 + k*k1[2]/2, dp2 + k*k1[3]/2))
    k3 = np.array(Jacobdelta(q1,q2,p1,p2,dq1 + k*k2[0]/2, dq2 + k*k2[1]/2, dp1 + k*k2[2]/2, dp2 + k*k2[3]/2))
    k4 = np.array(Jacobdelta(q1,q2,p1,p2,dq1 + k*k3[0], dq2 + k*k3[1], dp1 + k*k3[2], dp2 + k*k3[3]))

    new = ancien + (k/6) * (k1 + 2*k2 + 2*k3 + k4)
    print(new)
    return new
RK(q1,q2,p1,p2, d1,d2,d3,d4)
