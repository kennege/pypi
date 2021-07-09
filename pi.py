# https://mathscholar.org/2019/02/simple-proofs-archimedes-calculation-of-pi/
# https://www.davidhbailey.com/dhbpapers/pi-formulas.pdf

import math
import numpy as np

###### Archimedes Method (1)
def archimedes():
    k = 0 
    A = 2*(3**0.5) 
    B = 3 

    err = 1e-20
    while A-B > err:
        A = 2*A*B/(A+B)
        B = math.sqrt(A*B)
        k += 1

    pi = (A+B)/2
    return pi

####### Brent-Salamin algorithm (2)
def brent_salamin(K=6):
    k = 0
    A = 1
    B = 1/math.sqrt(2)
    S = 1/2

    while k<K:
        Ak = (A+B)/2
        Bk = (A*B)**0.5
        Ck = Ak**2 - Bk**2
        Sk = S - 2**(k+1)*Ck
        Pk = (2*Ak**2)/Sk

        A = Ak
        B = Bk
        S = Sk
        P = Pk
        k+=1

    pi = P
    return pi

####### Borwein cubic algorithm (3)
def borwein_cubic(K=6):
    k = 0
    A = 1/3
    S = ((3**0.5) - 1)/2

    while k<K:
        R = 3/(1+2*(1-S**3)**(1/3))
        S = (R-1)/2
        A = R**2*A - 3**k*(R**2-1)
        k+=1
    
    pi = 1/A
    return pi

####### Borwein quartic algorithm (4)
def borwein_quartic(K=3):
    k = 0 
    A = 6 - 4*(2**0.5)
    B = (2**0.5) - 1

    while k<K:
        B = (1-(1-B**4)**(1/4))/(1+(1-B**4)**(1/4))
        A = A*(1+B)**4 - 2**(2*k+3)*B*(1+B+B**2)
        k+=1

    pi = 1/A
    return pi

####### Borwein nonic algorithm (5)
def borwein_nonic(K=3):
    k = 0
    A = 1/3
    R = ((3**0.5)-1)/2
    S = (1-R**3)**(1/3)

    while k<K:
        T = 1+2*R
        U = (9*R*(1+R+R**2))**(1/3)
        V = T**2 + T*U+U**2
        W = (27*(1+S+S**2))/V
        A = W*A+3**(2*k-1)*(1-W)
        S = ((1-R)**3)/((T+2*U)*V)
        R = (1-S**3)**(1/3)
        k+=1
    
    pi = 1/A
    return pi

######## Leibniz (6)
def leibniz(K=1e6):
    k = 0
    A = 0
    while k<K:
        A += (((-1)**k)/(2*k+1))
        k+=1
    pi = 4*A
    return pi

######## Euler (7)
def euler(self,K=1e6):
    k = 1
    A = 0
    while k<K:
        A += (1/(k**2))    
        k+=1
    pi = (6*A)**0.5

    return pi

######## Chudnovsky (8)
def chudnovsky(self,K=2):
    k = 0
    A = 0
    while k<K:
        num = ((-1)**k)*(math.factorial(6*k))*(13591409+(545140134*k))
        den = (math.factorial(3*k))*(math.factorial(k)**3)*(640320**(3*k+(3/2)))
        A += (num/den)
        k+=1
    pi = 1/(12*A)
    return pi


######## intergration (9)
def integrate(K=1e6):
    k=0
    x = -1
    A = 0
    dx = 2/K
    while k<K:
        x+= dx
        A += (math.sqrt(1-x**2))*dx
        k+=1
    pi = 2*A
    return pi

######## Viete 1593 (10)
def viete(K=100):
    k = 1
    num = 0
    A=1
    while k<K:
        num = math.sqrt(2 + num)
        A *= num/2
        k+=1
    pi = 2/A
    return pi

####### John Wallis 1655 (11)
def wallis(K=1e6):
    k = 1
    A = 1
    while k<K:
        A *= ((2*k)/(2*k-1))*((2*k)/(2*k+1))
        k+=1
    pi = 2*A
    return pi

###### Gregory-Leibniz (12)
def g_leibniz(K=1e6):
    k = 1
    A = 0
    while k<K:
        A += ((-1)**(k+1) * 4)/(2*k-1)
        k+=1
    pi = A
    return pi

###### Nilakantha (13)
def nilakantha(K=1000):
    k = 1
    A = 3
    while k<K:
        A += (-1)**(k+1) * 4/(((2*k)*(2*k+1)*(2*k+2)))
        k+=1
    pi = A
    return pi 

###### Riemann zeta(2) (14)
def riemann_2(K=1e6):
    k=1
    A = 0
    while k<K:
        A+=1/(k**2)
        k+=1
    pi = (6*A)**0.5
    return pi

###### Riemann zeta(4) (15)
def riemann_4(K=1e6):
    k=1
    A = 0
    while k<K:
        A+=1/(k**4)
        k+=1
    pi = (90*A)**0.25
    return pi   

###### Plouffe 2006 (16)
def plouffe(K=60):
    k=1
    A = 0
    B = 0
    C = 0
    while k<K:
        A+=1/(k*(np.exp(k*math.pi))-1)
        B+=1/(k*(np.exp(2*k*math.pi))-1)
        C+=1/(k*(np.exp(4*k*math.pi))-1)
        k+=1
    return 72*A - 96*B + 24*C

###### Double Integral (17)
def double_integral(K=5):
    # need to use trapezoidal rule
    x = -K
    t = -K
    dx = 2
    dt = 2
    A = 0
    while t<dt*K-1:
        x=t
        while x<dx*K-1:        
            if x != 0:
                A += math.exp(-t**2-(1/(2*x**2))+x*t)
            x+=dx
        t+=dt
    return A
        
###### continued fraction (18) 
def continued(K=100):
    A = 0
    k=K
    while k>0:
        num = (2*k-1)**2
        A = num / (6 + A)
        k-=1
    pi = 3 + A
    return pi

###### sin improper integral (19)
def integrate_sin(K=1e6):
    k=0
    x = -K
    A = 0
    dx = 2
    while k<K:
        if x == 0:
            x=1e-20
        A += (math.sin(x)/x)*dx
        x+=dx
        k+=1
    pi = A
    return pi    

if __name__ == "__main__":
    
    print(math.pi)
    print(integrate_sin())

