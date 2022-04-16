import matplotlib.pyplot as plt
import numpy as np
from math import asin, log, sqrt
import bisect
from sage.plot.line import line
from sage.plot.point import point
from sage.calculus.integration import numerical_integral as integral_numerical
from sage.libs.pari import pari
from sage.rings.fast_arith import prime_range

def fast_aplist(p):
    """
    compute a list of the a_p value of all elliptic curves over F_p
    """
    ap_list = np.array([])
    for a in range(0, int(2 * sqrt(p))):
        val = int(4*p - a**2)
        hurwitz_number = int(pari(int(4*p - a**2)).qfbhclassno()) #need to change sign because of implementaion of hurwitz class number        
        ret_list = np.zeros(hurwitz_number)
        ret_list += float(a/float(2 * sqrt(p)))
        ap_list = np.append(ap_list, ret_list)
        if a != 0:
            ap_list = np.append(ap_list, -ret_list)
    return ap_list

def histogram(num_bins, p):
    """
    plot histogram of v, divided in num_bins bins.
    p is only given to be added to title of histogram
    """
    v = fast_aplist(p)
    n, bins, patches = plt.hist(v, num_bins, density=True)
    angle = np.linspace(0, np.pi, 150) 
    radius = 1
    x = radius * np.cos(angle) 
    y = radius * np.sin(angle) * 2 / np.pi #stretching so that area = 1
    plt.plot(x, y)
    plt.xlabel('$a_p/2\\sqrt{p}$')
    plt.ylabel('Frequency')
    plt.title('Vertical Sato-Tate-Distribution for p = {}'.format(p))
    plt.xlim(-1,1)
    plt.grid(True)
    plt.show()


def liney(y, xmin,xmax):
    return line([(xmin,y),(xmax,y)], rgbcolor=(1,0,0))

def Xab(a,b):
    bb = (asin(b)/2.0 + b*sqrt(1.0-b**2.0)/2.0)
    aa = (asin(a)/2.0 + a*sqrt(1.0-a**2.0)/2.0)
    def X(T):
        return (asin(T)/2 + T*sqrt(1-T**2)/2 - aa)/(bb - aa)
    return X

def sorted_aplist(p):
    v = fast_aplist(p)
    v.sort()
    return v

def Ypab(p, a=-1, b=1):
    v = sorted_aplist(p)
    denom = bisect.bisect_right(v, float(b)) - bisect.bisect_left(v, float(a))
    try:
        normalize = float(1)/denom
    except:
        def Y(T):
            return 1.0
        return Y
    start_pos = bisect.bisect_left(v, float(a))

    def Y(T):
        # find position that T would go in if it were inserted
        # in the sorted list v.
        n = bisect.bisect_right(v, float(T)) - start_pos
        return n * normalize
    return Y

global _delta

def Delta(p, a, b, max_points=300, verb=False, norm='l2'):
    """
    Delta_{a}**{b} function:
    INPUT: C - cutoff
    a,b - evaluate over the interval (a,b)
    max_points - number of points used in numerical integral
    """
    global _delta
    key = (p,a,b,max_points)
    #do not cache while changing the function!
    try:
        return _delta[key]
    except NameError:
        _delta = {}
    except KeyError:
        pass
    X = Xab(a,b)
    Y = Ypab(p,a,b)
    def h(T):
        if norm=='l2':
            return (X(T) - Y(T))**2.0
        if norm=='l1':
            return abs(X(T) - Y(T))
    
    val, err = integral_numerical(h, a, b, max_points=max_points, algorithm='qag', rule=1, eps_abs=1e-10, eps_rel=1e-10)
    
    _delta[key] = (val, err)
    if verb:
        print("p: ", p, " val: ", val)
    if val < 0:
        val = -val
        print(val)
    return val, err

def plot_Delta(pmax, step_size=10, max_points=100, a=-1, b=1):
    v = [(p,Delta(p, a, b, max_points=max_points)[0]) for p in prime_range(pmax)[4:][::step_size]] #skip 2,3,5,7 and then take only every step_size-th prime
    return line(v,rgbcolor=(0,0,0), ymin=0, ymax=0.1)


def theta(p, a=-1, b=1, max_points=300, verb=False):
    val, err = Delta(p, a, b, max_points=max_points, verb=verb)
    if verb: print(val, err)
    return -log(val)/log(p), val, err

def theta_interval(p, a=-1, b=1, max_points=300):
    val, err = Delta(p, a, b, max_points=max_points)
    return -log(val-abs(err))/log(p), -log(val+abs(err))/log(p)

def compute_theta(p, step_size=None, a=-1, b=1, max_points=300, verb=False):
    if verb:
        print('verbose')
    a,b = (float(a), float(b))
    if not step_size:
        step_size = max(1, int(p/(20 * log(p)))) # grows with pi(p)

    def f(p):
        if verb: print(p)
        z = theta(p, a, b, max_points=max_points, verb=verb)
        if verb: print(z)
        return z[0]

    #calls theta (= log etc. angewandt auf Delta)
    return [(x,f(x)) for x in prime_range(p)[25:][::step_size]]#skip the first 25 primes and then take only every step_size-th prime

def compute_theta_interval(p, step_size=None, a=-1, b=1, max_points=300, verb=False):
    if verb:
        print('verbose')
    a,b = (float(a), float(b))
    if not step_size:
        step_size = max(1, int(p/(20 * log(p)))) # grows with pi(p)
    vmin = []; vmax = []
    for C in prime_range(p)[25:][::step_size]:#skip the first 25 primes and then take only every step_size-th prime
        zmin,zmax = theta_interval(C, a, b, max_points=max_points)
        vmin.append((C, zmin))
        vmax.append((C, zmax))
        if verb: 
            print(C, zmin, zmax)
    return vmin, vmax

def plot_theta_interval(p, clr=(0,0,0), *args, **kwds):
    vmin, vmax = compute_theta_interval(p, *args, **kwds)
    v = compute_theta(p, *args, **kwds)
    grey = (0.7,0.7,0.7)
    p = line(vmin,rgbcolor=grey, ymin=0,ymax=1.5) + line(vmax,rgbcolor=grey) + point(v,rgbcolor=clr) + line(v,rgbcolor=clr) + liney(1,0, p)
    p.axes_labels(['$p$', '$\\theta$'])

def plot_theta(p, clr=(0,0,0), *args, **kwds):
    v = compute_theta(p, *args, **kwds)
    return point(v,rgbcolor=clr, ymin=0, ymax=2) + line(v, rgbcolor=clr) + liney(1,0,p)
