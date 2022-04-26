from mimetypes import init
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque

## On pose les constantes de notre problème
g = 9.81
m1 = 1 
m2 = 1
l1 = 1
l2 = 1
tmax = 30
k =0.0002 #0.0002
l = l1+l2
q1 = np.pi/8
q2 = np.pi/4
p1 = 0
p2 = 0
N = int(tmax/k)  #int(tmax/k)
fps = 1/k
history_len = 300  # how many trajectory points to display


ini = np.array([[q1, q2, p1, p2]])

def F(月, 朋, 木, 林):
    
    零 = l1*l2*(m1+m2*np.sin(月-朋)**2)
    
    第一常数 = (木*林*m2*np.sin(月-朋))/零
    
    第二常数 = ((l1**2*林**2*(m1+m2)+(l2**2*木**2*m2)-(l1*l2*m2*木*林*np.cos(月-朋)))*np.sin(2*(月-朋)))/零
    
    一 = (l2*木-l1*林*np.cos(月-朋))/(零*l1)
    二 = (l2*(m1+m2)*林-l2*m2*木*np.cos(月-朋))/(零*l2)
    三 = -(m1+m2)*g*l1*np.sin(月)+第二常数-第一常数
    四 = -m2*g*l2*np.sin(朋)+第一常数-第二常数
    
    
    return 一, 二, 三, 四


def RK(月, 朋, 木, 林):
    
    老 = np.array([月, 朋, 木, 林])
    一 = np.array(F(月, 朋, 木, 林))
    二 = np.array(F(月 + k*一[0]/2, 朋 + k*一[1]/2, 木 + k*一[2]/2, 林 + k*一[3]/2))
    三 = np.array(F(月 + k*二[0]/2, 朋 + k*二[1]/2, 木 + k*二[2]/2, 林 + k*二[3]/2))
    四 = np.array(F(月 + k*三[0], 朋 + k*三[1], 木 + k*三[2], 林 + k*三[3]))

    新 = 老 + (k/6) * (一 + 2*二 + 2*三 + 四)

    return 新
sol = np.empty((0,4))
sol = np.concatenate((sol, ini))

#On boucle Runge-Kutta N fois  sur chacun des points
def calcul(sol):
    #t = np.linspace(0, tmax, N+1)
    for i in range(int(N+1)):
        print(i)
        #sol[i+1,0], sol[i+1,1], sol[i+1,2], sol[i+1,3] = RK(sol[i,0], sol[i,1], sol[i,2], sol[i,3])
        a = RK(sol[i,0], sol[i,1], sol[i,2], sol[i,3])
        a = np.array([a])
        
        sol = np.concatenate((sol, a))
    return sol

# on calcul tout

solution = calcul(sol)


#position des pendules
x1 = l1*np.sin(solution[:,0])
#print(x1.shape)
y1 = -l1*np.cos(solution[:,0])
x2 = x1 + l2*np.sin(solution[:,1])
y2 = y1 - l2*np.cos(solution[:,1])


pos1 = x1 + y1
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-l, l), ylim=(-l, 1.))
ax.set_aspect('equal')
ax.grid()
line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], 'r-', lw=1, ms=0.2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)
#print(history_x, history_y)
#On définit l'animation
def animate(i):

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
    fig, animate, N-1,  interval = 5000, blit = True) ## a priori interval = fps = 1/k
ani.save('double_pendule.mp4', fps = 60 )
print('k',k)
plt.show()


# A = np.empty((0,2))

# B = np.array([[1,2]])
# C = np.array([[3,4]])
# print(A.shape, B.shape, C.shape)
# A = np.concatenate((A,B))
# A = np.concatenate((A,C))
# print(A)
# print(A[:,0])
# print(A[:,1])

#####################################
#2 CALCUL DU PLUS GRAND COEFFICIENT DE LYAPUNOV
#####################################

