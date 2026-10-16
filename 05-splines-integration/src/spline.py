import numpy as np
from load_foil import load_foil
import matplotlib.pyplot as plt
import sys

def spline(x, y, n, yp1, ypn) :
    p = qn = sig = un = 0.0
    u = np.zeros(n)
    y2 = np.zeros(n)
    
    if yp1 > 0.99*10**30 :
        y2[0] = 0.0
        u[0] = 0.0
    else :
        y2[0] = -0.5
        u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1)
    for i in range(1, n-1) :
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = sig*y2[i-1]+2.0
        y2[i] = (sig-1.0)/p
        u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
    if ypn > 0.99*10**30 :
        qn = un = 0.0
    else :
        qn = 0.5
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0)
    for k in range(n-2, -1, -1) :
        y2[k] = y2[k]*y2[k+1]+u[k]
    return y2


def splint(xa, ya, y2a, n, x) :
    klo = 0
    khi = n-1
    k = 0
    
    while khi-klo > 1 :
        k = int((khi + klo) / 2)
        if xa[k] > x :
            khi = k
        else :
            klo = k
            
    h = xa[khi] - xa[klo]
    if h == 0 :
        print("Bad xa input to routine splint")
        return
    
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b-b)*y2a[khi])*(h*h)/6
    
    return y

def test_spline_with_sin() :
    x = np.linspace(-2*np.pi, 2*np.pi, 1000)
    y = np.sin(x)
    n = len(x)
    yp1 = np.cos(x[0])
    ypn = np.cos(x[n-1])
    fx = spline(x, y, n, yp1, ypn)

    plt.plot(x, fx, label='Dérivée seconde de sin avec spline', color='black')
    plt.plot(x, -np.sin(x), '--', label='Dérivée seconde de sin calculée manuellement : -sin(x)', color='red')
    plt.xlabel("x")
    plt.ylabel("Dérivée seconde de sin(x)")
    plt.legend()
    plt.title("Test spline sur sin(x)")
    plt.grid(True)
    plt.show()

def test_spline_with_x2() :
    x = np.linspace(-10, 10, 1000)
    y = x**2
    n = len(x)
    yp1 = 2*x[0]
    ypn = 2*x[n-1]
    fx = spline(x, y, n, yp1, ypn)
    fx_naturel = spline(x, y, n, 1e30, 1e30)

    plt.plot(x, fx, label='Dérivée seconde de x2 avec spline', color='black')
    #plt.plot(x, fx_naturel, label='Dérivée seconde de x2 avec spline naturel', color='blue')
    f = [2 for i in range(n)]
    plt.plot(x, f, '--', label='Dérivée seconde de x2 calculée manuellement : 2', color='red')
    plt.xlabel("x")
    plt.ylabel("Dérivée seconde de x2")
    plt.legend()
    plt.title("Test spline sur x2")
    plt.grid(True)
    plt.show()

def test_spline_with_x3() :
    x = np.linspace(-10, 10, 1000)
    y = x**3
    n = len(x)
    yp1 = 3*x[0]**2
    ypn = 3*x[n-1]**2
    fx = spline(x, y, n, yp1, ypn)
    fx_naturel = spline(x, y, n, 1e30, 1e30)

    plt.plot(x, fx, label='Dérivée seconde de x3 avec spline', color='black')
    #plt.plot(x, fx_naturel, label='Dérivée seconde de x3 avec spline naturel', color='blue')
    plt.plot(x, 6*x, '--', label='Dérivée seconde de x3 calculée manuellement : 6x', color='red')
    plt.xlabel("x")
    plt.ylabel("Dérivée seconde de x3")
    plt.legend()
    plt.title("Test spline sur x3")
    plt.grid(True)
    plt.show()
    
def test_airfoil() :
    
    if len(sys.argv) == 1:
        print( "Please provide a input file")
        exit(1)
    
    (dim,ex,ey,ix,iy) = load_foil( sys.argv[1] )
    
    ne = int(dim[0])
    ye = spline(ex, ey, ne, ey[0], ey[ne-1])
    def fe(x) :
        return splint(ex, ey, ye, ne, x)
    fex = [ fe(i) for i in ex ]
    
    ni = int(dim[1])
    yi = spline(ix, iy, ni, iy[0], iy[ni-1])
    def fi(x) :
        return splint(ix, iy, yi, ni, x)
    fix = [ fi(i) for i in ix ]

    plt.plot(ex, fex)
    plt.plot(ix, fix)
    plt.ylim(-0.3,0.3)
    plt.xlim(0, 1)
    plt.xticks(np.arange(0,1,0.1))
    plt.grid()
    plt.show()
    

if __name__ == "__main__" :
    #test_spline_with_sin()
    #test_spline_with_x2()
    test_spline_with_x3()
    test_airfoil()
