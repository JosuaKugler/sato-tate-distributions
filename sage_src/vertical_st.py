#Functions:
#
#-fast_aplist(p)
#   compute a list of the a_p value of all elliptic curves over F_p
#-sorted_aplist(p)
#   sort ap_list(p)
#-histogram(p)
#   plot histogram of fast_aplist(p), divided in num_bins bins
#-Xab (from python_utils)
#   return a function X(T)
#    where X(T) is the area under the arc of the semicircle from a to T
#-Ypab
#   return a function Y(T)
#   where Y(T) is the area under the 'histogram' from a to T
#-Delta
#compute Delta_a^b(p)
#   where max_points is the number of points used in numerical integration
#   and norm is either 'l2' or 'l1'
#-plot_Delta
#   plot Delta_a^b(p) from 0 to pmax (use every stepsize-th prime)
#-theta
#   compute theta from Delta
#-theta_error_bound
#   compute error bound from numerical integration
#-theta_range
#   compute list of (x, theta(x)) pairs where x is a prime <= p
#   use every stepsize-th prime for x
#   where max_points is the number of points used in numerical integration
#   if step_size is not given the function estimates a proper size depending on p
#-theta_error_bound_range
#   compute lists of (x, theta_max(x)) or (x, theta_min(x)) pairs where x is a prime <= p
#   use every stepsize-th prime for x
#   where max_points is the number of points used in numerical integration
#   and theta_max/theta_min are upper/lower error bounds for the result of the numerical integration
#   if step_size is not given the function estimates a proper size depending on p
#-plot_theta
#plot theta(x) for x prime from 0 to p
#   if show_error_bound is True the upper and lower error_bound is shown

import matplotlib.pyplot as plt
import numpy as np
from math import log, sqrt
import bisect
from sage.plot.line import line
from sage.plot.point import point
from sage.calculus.integration import numerical_integral as integral_numerical
from sage.libs.pari import pari
from sage.rings.fast_arith import prime_range
from utils import Xab

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

def sorted_aplist(p):
    """
    sort ap_list(p)
    """
    v = fast_aplist(p)
    v.sort()
    return v

def histogram(num_bins, p):
    """
    plot histogram of fast_aplist(p), divided in num_bins bins
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

def Ypab(p, a=-1, b=1):
    """
    return a function Y(T)
    where Y(T) is the area under the 'histogram' from a to T
    """
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

def Delta(p, a, b, max_points=300, norm='l2'):
    """
    compute Delta_a^b(p)
    where max_points is the number of points used in numerical integration
    and norm is either 'l2' or 'l1'
    """
    global _delta
    key = (p,a,b,max_points)
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
    return val, err

def plot_Delta(pmax, step_size=10, max_points=100, a=-1, b=1):
    """
    plot Delta_a^b(p) from 0 to pmax (use every stepsize-th prime)
    """
    v = [(p,Delta(p, a, b, max_points=max_points)[0]) for p in prime_range(pmax)[4:][::step_size]] #skip 2,3,5,7 and then take only every step_size-th prime
    p = line(v,rgbcolor=(0,0,0), ymin=0, ymax=0.1)
    p.axes_labels(['$p$', '$\\Delta$'])
    return p

def theta(p, a=-1, b=1, max_points=300):
    """
    compute theta from Delta
    """
    val, err = Delta(p, a, b, max_points=max_points)
    return -log(val)/log(p), val, err

def theta_error_bound(p, a=-1, b=1, max_points=300):
    """
    compute error bound from numerical integration
    """
    val, err = Delta(p, a, b, max_points=max_points)
    return -log(val-abs(err))/log(p), -log(val+abs(err))/log(p)

def theta_range(p, step_size=None, a=-1, b=1, max_points=300):
    """
    compute list of (x, theta(x)) pairs where x is a prime <= p
    use every stepsize-th prime for x
    where max_points is the number of points used in numerical integration
    if step_size is not given the function estimates a proper size depending on p
    """
    a,b = (float(a), float(b))
    if not step_size:
        step_size = max(1, int(p/(20 * log(p)))) # grows with pi(p)

    def f(p):
        z = theta(p, a, b, max_points=max_points)
        return z[0]

    #calls theta (= log etc. angewandt auf Delta)
    return [(x,f(x)) for x in prime_range(p)[25:][::step_size]]#skip the first 25 primes and then take only every step_size-th prime

def theta_error_bound_range(p, step_size=None, a=-1, b=1, max_points=300):
    """
    compute lists of (x, theta_max(x)) or (x, theta_min(x)) pairs where x is a prime <= p
    use every stepsize-th prime for x
    where max_points is the number of points used in numerical integration
    and theta_max/theta_min are upper/lower error bounds for the result of the numerical integration
    if step_size is not given the function estimates a proper size depending on p
    """
    a,b = (float(a), float(b))
    if not step_size:
        step_size = max(1, int(p/(20 * log(p)))) # grows with pi(p)
    vmin = []; vmax = []
    for C in prime_range(p)[25:][::step_size]:#skip the first 25 primes and then take only every step_size-th prime
        zmin,zmax = theta_error_bound(C, a, b, max_points=max_points)
        vmin.append((C, zmin))
        vmax.append((C, zmax))
    return vmin, vmax

def plot_theta(p, clr=(0,0,0), show_error_bound=True, *args, **kwds):
    """
    plot theta(x) for x prime from 0 to p
    if show_error_bound is True the upper and lower error_bound is shown
    """
    if show_error_bound:
        vmin, vmax = theta_error_bound_range(p, *args, **kwds)
    v = theta_range(p, *args, **kwds)
    grey = (0.7,0.7,0.7)
    theta_line = line(v,rgbcolor=clr,ymin=0,ymax=1.5, legend_label='$\\theta(p)$')
    theta_points = point(v,rgbcolor=clr)
    if show_error_bound:
        error_bound = line(vmin,rgbcolor=grey, legend_label='numerical integration error bound') + line(vmax,rgbcolor=grey)
    reference_line = line([(0,1),(p,1)], rgbcolor=(1,0,0), legend_label='$\\theta=1.0$')
    
    p = theta_line + theta_points + reference_line
    if show_error_bound:
        p = p + error_bound

    p.axes_labels(['$p$', '$\\theta$'])
    p.legend(True)

    return p