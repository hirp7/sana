# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 18:34:23 2022

@author: hiro7
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad
from matplotlib.animation import FuncAnimation


pi,sin,cos,exp,sqrt = [np.pi,np.sin,np.cos,np.exp,np.sqrt]




def main():
    
    #f = lambda x:1 if 0.4<=x<=0.6 else 0
    #Amp = 0.1
    f1 = lambda x:4/5/pi*x - 1/5 if pi/4<=x<=pi/2 else -4/5/pi*x + 3/5 if pi/2 <=x<=3*pi/4 else 0  
    #g = lambda x:0.5 if 0.4<=x<=0.6  else 0    
    g1 = lambda x: 0
    analytical1 = lambda x,t: 1.6/pi**2*((2-sqrt(2))*cos(t)*sin(x) - 1/9*(2+sqrt(2))*cos(3*t)*sin(3*x) +\
        1/25*(2+sqrt(2))*cos(5*t)*sin(5*t))
    string = PDE_1D(0,pi,100,0,1,100,f,g,1,5,1,analytical1)
    
    fig,ax = plt.subplots()
    ax.plot(string.x,string.f_x)
    
    string.animation("test")
    
    
    return 0

class PDE_1D(object):
    def __init__(self,x0,x1,nx,t0,t1,nt,f,g,c,modes,analytical_flag = None,analysis = None):
        self.x = np.linspace(x0,x1,nx)
        self.t = np.linspace(t0,t1,nt)
        self.f_x = np.array([f(i) for i in self.x])
        self.g_x = np.array([g(i) for i in self.x])
        self.sol_all = np.array([[u(i_x,i_t,f,g,L,1,1) for i_x in self.x] for i_t in self.t]) #[time frame*x*modes]
        self.sol = self.sol_all.sum(axis=2)
        self.analytical_flag = analytical_flag
        if analytical_flag != None:
            self.analysis = np.array([[analysis(i,j) for i in self.x] for j in self.t])
        else:
            self.analysis = None
        
    #def plot(self,mode):
        
    def plot(self,mode):
        return 0
        
    def animation(self,name):        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        anim = FuncAnimation(fig, lambda x:frame_update(x,ax,self.x,self.sol,self.analytical_flag,self.analysis), frames=np.arange(self.t.shape[0]), interval=200)
        
        anim.save(name+".gif", writer="imagemagick")
        plt.close()
        return 0


#https://www.yutaka-note.com/entry/matplotlib_artist_anim#Matplotlib%E3%81%A7%E3%82%A2%E3%83%8B%E3%83%A1%E3%83%BC%E3%82%B7%E3%83%A7%E3%83%B3%E3%82%92%E4%BD%9C%E6%88%90%E3%81%99%E3%82%8B%EF%BC%92%E3%81%A4%E3%81%AE%E6%96%B9%E6%B3%95
def u(x,t,f,g,L,c,n): 
    """
    f(x),g(x) is a function
    c is a constant c**2 = T/roh T:power, roh: mass per length of the string
    """
    b1 = lambda n: 2/L*quad(lambda x:f(x)*sin(n*pi*x/L),0,L)[0]
    b2 = lambda n: 2/(c*n*pi)*quad(lambda x:g(x)*sin(n*pi*x/L),0,L)[0]
    #lambdan = c*n*pi/L
    u_sum = np.array([(b1(i)*cos(c*i*pi/L*t) + \
                       b2(i)*sin(c*i*pi/L*t))*sin(i*pi*x/L)\
                      for i in (np.arange(n)+1) ])
    return u_sum

def frame_update(frame,ax,x,y,analytical_flag = None, analysis = None):
    ax.cla() # ax をクリア
    ax.plot(x,y[frame,:],'red')
    #ax.plot(x,test_ana[frame],'blue',linestyle = '--')
    ax.set_ylim(-y.max(),y.max())
    #ax.text(0.5,0.5,f'{t[frame]}')

    if analytical_flag != None:
        ax.plot(x,analysis[frame],'blue',linestyle = '--')
        ax.legend(['numerical','analytical'],loc = 'upper right')

    


if __name__ == '__main__':    
    main()




L = pi
#f = lambda x:1 if 0.4<=x<=0.6 else 0
Amp = 0.1
f = lambda x:4/5/pi*x - 1/5 if pi/4<=x<=pi/2 else -4/5/pi*x + 3/5 if pi/2 <=x<=3*pi/4 else 0  
#g = lambda x:0.5 if 0.4<=x<=0.6  else 0    
g = lambda x: 0
x = np.linspace(0,L,100)

fig,ax = plt.subplots(2,1)
ax[0].plot(x,[f(i) for i in x])    
ax[1].plot(x,[g(i) for i in x])    

t = np.linspace(0,5,20)    
test_frames = np.array([[u(i_x,i_t,f,g,L,1,1) for i_x in x] for i_t in t])
test = test_frames.sum(axis = 2)

fig = plt.figure()
ax = fig.add_subplot(111)

def frame_update(frame,ax,analytical = None):
    ax.cla() # ax をクリア
    ax.plot(x,test[frame,:],'red')
    #ax.plot(x,test_ana[frame],'blue',linestyle = '--')
    ax.set_ylim(-test.max(),test.max())
    if analytical != None:
        ax.plot(x,analytical[frame],'blue',linestyle = '--')
        ax.legend(['numerical','analytical'],loc = 'upper right')
        

anim = FuncAnimation(fig, frame_update, frames=np.arange(t.shape[0]), interval=200)

anim.save("string.gif", writer="imagemagick")
plt.close()
    
"""
for d in frames:
   artists = func(d, *fargs)
   fig.canvas.draw_idle()
   fig.canvas.start_event_loop(interval)
"""
"""
thermal simulation
diffusion equation
"""
def u(x,t,f,L,c,modes): 
    """
    f(x),g(x) is a function
    c is a constant c**2 = T/roh T:power, roh: mass per length of the string
    """
    b1 = lambda n: 2/L*quad(lambda x:f(x)*sin(n*pi*x/L),0,L)[0]
    #b2 = lambda n: 2/(c*n*pi)*quad(lambda x:g(x)*sin(n*pi*x/L),0,L)[0]
    #lambdan = c*n*pi/L
    u_sum = np.array([b1(i)*sin(i*pi/L*x)*exp(-(c*i*pi/L)**2*t)\
                      for i in (np.arange(modes)+1) ])
    return u_sum


L = 80
#f = lambda x:1 if 0.4<=x<=0.6 else 0
f = lambda x:100*sin(pi*x/80)
x = np.linspace(0,L,100)
c = sqrt(1.158)
fig,ax = plt.subplots()
ax.plot(x,[f(i) for i in x])    

t = np.linspace(0,700,20)    
test_frames = np.array([[u(i_x,i_t,f,L,c,3) for i_x in x] for i_t in t])
test = test_frames.sum(axis = 2)

u_ana = lambda x,t:100*sin(pi*x/80)*exp(-0.001785*t)
test_ana = np.array([[u_ana(i_x,i_t) for i_x in x] for i_t in t])





fig = plt.figure()
ax = fig.add_subplot(111)

def update(frame):
    ax.cla() # ax をクリア
    ax.plot(x,test[frame,:],'red')
    ax.plot(x,test_ana[frame],'blue',linestyle = '--')
    ax.set_ylim(0,100)
    ax.legend(['numerical','analytical'],loc = 'upper right')
    ax.text(0.5,0.5,f'{t[frame]}')

anim = FuncAnimation(fig, update, frames=np.arange(t.shape[0]), interval=200)

anim.save("c03.gif", writer="imagemagick")
plt.close()
    
def u_2D(x,y,t,f,L,c,modes): 
    """
    f(x),g(x) is a function
    c is a constant c**2 = T/roh T:power, roh: mass per length of the string
    """
    b1 = lambda n: 2/L*quad(lambda x:f(x)*sin(n*pi*x/L),0,L)[0]
    #b2 = lambda n: 2/(c*n*pi)*quad(lambda x:g(x)*sin(n*pi*x/L),0,L)[0]
    #lambdan = c*n*pi/L
    u_sum = np.array([b1(i)*sin(i*pi/L*x)*exp(-(c*i*pi/L)**2*t)\
                      for i in (np.arange(modes)+1) ])
    return u_sum
