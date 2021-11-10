
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


x_0=0.01
v_0=0.
t_0=0.
S=0.75
p=2
alpha=5.
t_max=40.
dt=1.
t=0.

def a11():
    return 1.

def a12(dt):
    return -dt/2.

def a21(xn,vn,dt):
    return -dt/2. * (-2.0*alpha*xn*vn - 1.)

def a22(xn,dt):
    return 1.0 - dt/2.0*alpha*(1.0 - xn**2)

def g(x, v, alpha):
    return (alpha*(1.- x**2)*v - x)


def F(xn, vn, xn_prev, vn_prev, dt, alpha):
    return xn - xn_prev - dt/2. * (vn_prev + vn)


def G(xn, vn, xn_prev, vn_prev, dt, alpha):
    return vn - vn_prev - dt/2. * (g(xn_prev, vn_prev, alpha) + g(xn, vn, alpha))



def rk2(xn, vn, dt, alpha):
    k_1x = vn
    k_1v = alpha*(1.-xn**2)*vn - xn

    k_2x = vn + dt*k_1v
    k_2v = alpha * (1.-(xn + dt*k_1x)**2) * (vn + dt*k_1v) - (xn + dt*k_1x)
    
    xn_next = xn + dt/2.*(k_1x + k_2x)
    vn_next = vn + dt/2.*(k_1v + k_2v)

    return xn_next, vn_next


def trapezy(xn_prev, vn_prev, dt, alpha, delta=1e-10):
    xn= xn_prev
    vn= vn_prev
    while True:
        a_11 = a11()
        a_12 = a12(dt)
        a_21 = a21(xn,vn,dt)
        a_22 = a22(xn,dt)

        dx = (-F(xn, vn, xn_prev, vn_prev, dt, alpha)*a_22 + G(xn, vn, xn_prev, vn_prev, dt, alpha)*a_12) / (a_11*a_22 - a_12*a_21)
        dv = (-G(xn, vn, xn_prev, vn_prev, dt, alpha)*a_11 + F(xn, vn, xn_prev, vn_prev, dt, alpha)*a_21) / (a_11*a_22 - a_12*a_21)

        xn += dx
        vn += dv

        if abs(dx) < delta and abs(dv) < delta:
            break
    
    xn_next = xn_prev + dt/2.0*(vn_prev +vn)
    vn_next = vn_prev + dt/2.0*(g(xn_prev, vn_prev, alpha) + g(xn, vn, alpha))

    return xn_next, vn_next



def calculateE( p1, p2):
        return (p1 -p2)/(2**p-1.)


def solve(schemat_numeryczny,TOL):

        dt=1.
        t=0.
        xn=x_0
        vn=v_0

        t_tab = np.array([dt])
        dt_tab = np.array([dt])
        x_tab = np.array([xn])
        v_tab = np.array([vn])

        
        i=0
        while True:
            xn1_2, vn1_2 = schemat_numeryczny(x_tab[i], v_tab[i],dt,alpha)
            xn2_2, vn2_2 = schemat_numeryczny(xn1_2, vn1_2, dt,alpha)

            xn2_1, vn2_1 = schemat_numeryczny(x_tab[i], v_tab[i],2.0*dt,alpha)

            E_x = calculateE(xn2_2, xn2_1) 
            E_v = calculateE(vn2_2, vn2_1) 

            if max(abs(E_x), abs(E_v)) < TOL:
                t_tab = np.append(t_tab, t_tab[i] + 2*dt)
                x_tab = np.append(x_tab, xn2_2)
                v_tab = np.append(v_tab, vn2_2)
                dt_tab = np.append(dt_tab, dt)
                i+=1

            dt = ((S * TOL / max(abs(E_x), abs(E_v)))**( 1.0/(p + 1.0))) * dt

            if t_tab[i-1] >= t_max:
                break
        return x_tab, v_tab, dt_tab, t_tab

    
def plot_results(t_tab1,t_tab2,dt_tab1,dt_tab2,x_tab1, x_tab2, v_tab1, v_tab2):
    plt.figure(0)
    plt.plot(t_tab1, x_tab1, label="TOL = 1e-2")
    plt.plot(t_tab2, x_tab2, label="TOL = 1e-5")
    plt.legend()
    plt.show()
    plt.figure(1)
    plt.plot(t_tab1, v_tab1, label="TOL = 1e-2")
    plt.plot(t_tab2, v_tab2, label="TOL = 1e-5")
    plt.legend()
    plt.show()
    plt.figure(2)
    plt.plot(t_tab1, dt_tab1, label="TOL = 1e-2")
    plt.plot(t_tab2, dt_tab2, label="TOL = 1e-5")
    plt.legend()
    plt.show()
    plt.figure(3)
    plt.plot(x_tab1, v_tab1, label="TOL = 1e-2")
    plt.plot(x_tab2, v_tab2, label="TOL = 1e-5")
    plt.legend()
    plt.show()       

if __name__ == "__main__":
    TOL=1e-2

    x_tab1, v_tab1, dt_tab1, t_tab1=solve(rk2,TOL)
    x_tab2, v_tab2, dt_tab2, t_tab2=solve(trapezy,TOL)

    TOL=1e-5

    x_tab3, v_tab3, dt_tab3, t_tab3=solve(rk2,TOL)
    x_tab4, v_tab4, dt_tab4, t_tab4=solve(trapezy,TOL)

    plot_results(t_tab1,t_tab3,dt_tab1,dt_tab3,x_tab1,x_tab3,v_tab1,v_tab3)
    plot_results(t_tab2,t_tab4,dt_tab2,dt_tab4,x_tab2,x_tab4,v_tab2,v_tab4)