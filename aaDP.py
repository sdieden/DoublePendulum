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
tmax = 60
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
ini_modif = np.array([[qi1+0.001, qi2+0.001, pi1,pi2]])


def M(q1, q2,p1,p2):
    
    一 = p1
    二 = p2
    三 = (-g*(2*m1+m2)*np.sin(q1)-m2*g*np.sin(q1-2*q2)-2*np.sin(q1-q2)*m2*(p2**2*l2+p1**2*l1*np.cos(q1-q2)))/(l1*(2*m1+m2-m2*np.cos(2*q1-2*q2)))
    # 啊 = -g*(2*m1*m2)*np.sin(月)
    # print('tot', 三)
    # print('1',啊)
    # print('2',-m2*g*np.sin(月-2*朋) )
    # print('3',-2*np.sin(月-朋)*m2*(林**2*l2+木**2*l1*np.cos(月-朋)))
    # print('4',1/(l1*(2*m1+m2-m2*np.cos(2*月-2*朋))))
    四 = (2*np.sin(q1-q2)*(p1**2*l1*(m1+m2)+g*(m1+m2)*np.cos(q1)+p2**2*l2*m2*np.cos(q1-q2)))/(l2*(2*m1+m2-m2*np.cos(2*(q1-q2))))
    return 一,二,三,四



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
sol_modified = np.concatenate((sol, ini_modif))



#On boucle Runge-Kutta N fois  sur chacun des points
def calcul(sol):
    t = np.linspace(0, tmax, N+2)
    for i in range(int(N+1)):
        print(i)
        #sol[i+1,0], sol[i+1,1], sol[i+1,2], sol[i+1,3] = RK(sol[i,0], sol[i,1], sol[i,2], sol[i,3])
        a = RK4(sol[i,0], sol[i,1], sol[i,2], sol[i,3])
        a = np.array([a])
        
        sol = np.concatenate((sol, a))

    return sol, t

# on calcul tout

solution, t = calcul(sol)

#solution2, t = calcul(sol_modified)


#position des pendules #1
x1 = l1*np.sin(solution[:,0])
#print(x1.shape)
y1 = -l1*np.cos(solution[:,0])
x2 = x1 + l2*np.sin(solution[:,1])
y2 = y1 - l2*np.cos(solution[:,1])
y2s = 2 + y1 - l2*np.cos(solution[:,1])
#2
# x11 = l1*np.sin(solution2[:,0])
# y11 = -l1*np.cos(solution2[:,0])
# x22 = x11 + l2*np.sin(solution2[:,1])
# y22 = y11 - l2*np.cos(solution2[:,1])
# y2ss = 2 + y11 - l2*np.cos(solution2[:,1])

### définition de l'énergie

Ec1, Ec2 = np.zeros(N+2), np.zeros(N+2)
Ep1, Ep2 = np.zeros(N+2), np.zeros(N+2)
Em1, Em2 = np.zeros(N+2), np.zeros(N+2)
Em = np.zeros(N+2)
print(Em.shape)

def energie():
    print(Ec1)
    print()
    Ec1[0] = (1./2)*m1*(l1**2)*(solution[0,0]**2)
    Ec2[0]=(1./2)*m2*(l1**2*solution[0,2]**2+l2**2*solution[0,3]**2+2*l1*l2*solution[0,2]*solution[0,3]*np.cos(solution[0,0]-
solution[0,1]))
   
    Ep1[0]=-m1*g*l1*np.cos(solution[0,0])
    Ep2[0]=-m2*g*(l1*np.cos(solution[0,0])+l2*np.cos(solution[0,1]))
    
    Em2[0]=Ec2[0]+Ep2[0]
    Em1[0]=Ec1[0]+Ep1[0]
    
    for i in range(N+1):
        Ec1[i+1]=(1./2)*m1*(l1**2)*(solution[i,0]**2)
        Ec2[i+1]=(1./2)*m2*(l1**2*solution[i,2]**2+l2**2*solution[i,3]**2+2*l1*l2*solution[i,2]*solution[i,3]*np.cos(solution[i,0]-
solution[i,1]))
        
 
        Ep1[i+1]=-m1*g*l1*np.cos(solution[i,1])
        Ep2[i+1]=-m2*g*(l1*np.cos(solution[i,1])+l2*np.cos(solution[i,2]))
        
        Em1[i+1]=Ec1[i+1]+Ep1[i+1]
        Em2[i+1]=Ec2[i+1]+Ep2[i+1]

    Em = Em1 + Em2
    print(Em.shape)
    return Em

Em = energie()
####On crée le graphe pour l'animation 

pos1 = x1 + y1
# pos2 = x11 + y11
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
#print(history_x, history_y)
#On définit l'animation
#####################################################################################################################
######### activer pour avoir 2 pendules /!\/!\ PROBLEME LE 2EME PENDULE EST 2 FOIS PLUS LENT ET TREMBLE
#####################################################################################################################
# def animate(i):

#     thisx = [0, x1[i], x2[i]]
#     thisy = [0, y1[i], y2[i]]
#     thisxx = [0, x11[i], x22[i]]
#     thisyy = [0, y11[i], y22[i]]

#     if i == 0 : 
#         history_x.clear()
#         history_y.clear()
#         history_xx.clear()
#         history_yy.clear()
#     history_x.appendleft(thisx[2])
#     history_y.appendleft(thisy[2])
#     history_xx.appendleft(thisxx[2])
#     history_yy.appendleft(thisyy[2])
#     line.set_data(thisx, thisy)
#     line2.set_data(thisxx, thisyy)
#     trace.set_data(history_x, history_y)
#     trace2.set_data(history_xx, history_yy)
#     time_text.set_text(time_template % (i*k))
#     return line, time_text, trace, line2, trace2
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

ani = animation.FuncAnimation(
    fig, animate2, N-1,  interval = k*1000, blit = True) ## a priori interval = fps = 1/k
ani.save('double_pendule.mp4')
plt.show()
plt.close()

fig, ax2=plt.subplots()
ax2.plot(t, Em)
plt.title('Evolution energetique')
plt.xlabel('Temps (seconde)')
plt.ylabel('Energies (joule)')
plt.show()

#####################################
#2 CALCUL DU PLUS GRAND COEFFICIENT DE LYAPUNOV
#####################################


