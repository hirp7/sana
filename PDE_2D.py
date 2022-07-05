# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 21:13:39 2022

@author: hiro7
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad
from matplotlib.animation import FuncAnimation


pi,sin,cos,exp,sqrt = [np.pi,np.sin,np.cos,np.exp,np.sqrt]

def main():
    a = 1
    b = 1
    c = 1
    d = a+b+c
    print(f'{d}')
    return 0

if __name__ == '__main__':
    
    main()
    
    
def u_1D(x,t,f,L,c,modes): 
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
