import matplotlib.pyplot as plt
import numpy as np
from math import asin, log, sqrt
import bisect

def fast_aplist(p):
    """
    compute a list of the a_p value of all elliptic curves over F_p
    """
    ap_list = np.array([])
    for a in range(0, int(2 * sqrt(p))):
        hurwitz_number = int(pari(int(4*p - a^2)).qfbhclassno()) #need to change sign because of implementaion of hurwitz class number
        list = np.zeros(hurwitz_number)
        list += a/float(2 * sqrt(p))
        ap_list = np.append(ap_list, list)
        if a != 0:
            ap_list = np.append(ap_list, -list)
    return ap_list

def histogram(v, num_bins, p):
    """
    plot histogram of v, divided in num_bins bins.
    p is only given to be added to title of histogram
    """
    n, bins, patches = plt.hist(v, num_bins, density=True)
    angle = np.linspace(0, np.pi, 150) 
    radius = 1
    x = radius * np.cos(angle) 
    y = radius * np.sin(angle) * 2 / np.pi #stretching so that area = 1
    plt.plot(x, y)
    plt.xlabel('$a_p/2\\sqrt{p}$')
    plt.ylabel('Frequency')
    plt.title('Vertical Sato-Tate-Distribution for p = {}'.format(p))
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    plt.xlim(-1,1)
    #plt.ylim(0, 1)
    plt.grid(True)
    plt.show()


def liney(y, xmin,xmax):
    return line([(xmin,y),(xmax,y)], rgbcolor=(1,0,0))

def Xab(a,b):
    bb = (asin(b)/2.0 + b*sqrt(1.0-b^2.0)/2.0)
    aa = (asin(a)/2.0 + a*sqrt(1.0-a^2.0)/2.0)
    def X(T):
        return (asin(T)/2 + T*sqrt(1-T^2)/2 - aa)/(bb - aa)
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
    Delta_{a}^{b} function:
    INPUT: C - cutoff
    a,b - evaluate over the interval (a,b)
    max_points - number of points used in numerical integral
    """
    global _delta
    key = (p,a,b,max_points)
    #do not cache while changing the function!
    #try:
    #    return _delta[key]
    #except NameError:
    #    _delta = {}
    #except KeyError:
    #    pass
    X = Xab(a,b)
    Y = Ypab(p,a,b)
    def h(T):
        if norm=='l2':
            return (X(T) - Y(T))^2.0
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
    v = [(p,Delta(p, a, b, max_points=max_points)[0]) for p in prime_range(10,pmax)[::step_size]]
    return line(v,rgbcolor=(0,0,0), ymin=0, ymax=0.1)