from math import log, sqrt
import matplotlib.pyplot as plt
import numpy as np
import bisect
from sage.plot.line import line
from sage.plot.point import point
from sage.calculus.integration import numerical_integral as integral_numerical
from sage.rings.fast_arith import prime_range
from python_utils import Xab

class SatoTate:
    def __init__(self, E):
        self._E = E
        self._normalized_aplist = []

    def anlist(self, n):
        return self._E.anlist(n)

    def normalized_aplist_with_caching(self, n):
        primes = prime_range(n)
        pi = len(primes)
        if len(self._normalized_aplist) >= pi:
            return self._normalized_aplist[0:pi]
        anlist = self.anlist(n)
        v = [float(anlist[p])/float(2 * sqrt(p)) for p in primes]
        self._normalized_aplist = v
        return v.copy()
    
    def normalized_aplist(self, n):
        anlist = self.anlist(n)
        v = [float(anlist[p])/float(2 * sqrt(p)) for p in prime_range(n)]
        return v

    def sorted_aplist(self, n):
        #way faster, for 10^6: 4.1s instead of 60s
        v = self.normalized_aplist_with_caching(n)
        v.sort()
        return v
            
    def YCab(self, Cmax, a=-1, b=1):
        v = self.sorted_aplist(Cmax)
        
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
    
    
    def Delta(self, C, a, b, max_points=300):
        """
        Delta_{a}^{b} function:
        INPUT: C - cutoff
        a,b - evaluate over the interval (a,b)
        max_points - number of points used in numerical integral
        """
        key = (C,a,b,max_points)
        try:
            return self._delta[key]
        except AttributeError:
            self._delta = {}
        except KeyError:
            pass
        X = Xab(a,b)
        Y = self.YCab(C,a,b)
        def h(T):
            return (X(T) - Y(T))**2.0
        
        val, err = integral_numerical(h, a, b, max_points=max_points, algorithm='qag', rule=1, eps_abs=1e-10, eps_rel=1e-10)
        
        self._delta[key] = (val, err)
        return val, err
    
    def plot_Delta(self, Cmax, plot_points=400, max_points=100, a=-1, b=1):
        v = [(x,self.Delta(x, a, b, max_points=max_points)[0]) for x in range(0, Cmax, int(Cmax/plot_points))]
        p = line(v,rgbcolor=(0,0,0), ymin=0, ymax=0.1)
        p.axes_labels(['$C$', '$\\Delta$'])
        return p
    
    def theta(self, C, a=-1, b=1, max_points=300):
        val, err = self.Delta(C, a, b, max_points=max_points)
        return -log(val)/log(C), val, err
    
    def theta_interval(self, C, a=-1, b=1, max_points=300):
        val, err = self.Delta(C, a, b, max_points=max_points)
        return -log(val-abs(err))/log(C), -log(val+abs(err))/log(C)
    
    def compute_theta(self, Cmax, plot_points=30, a=-1, b=1, max_points=300, verbose=False):
        a,b = (float(a), float(b))
        
        def f(C):
            z = self.theta(C, a, b, max_points=max_points)
            if verbose: print(C, z)
            return z[0]

        #calls theta (= log etc. angewandt auf Delta)
        return [(x,f(x)) for x in range(100, Cmax, int(Cmax/plot_points))]
    
    def compute_theta_interval(self, Cmax, plot_points=30, a=-1, b=1, max_points=300, verbose=False):
        a,b = (float(a), float(b))
        vmin = []; vmax = []
        for C in range(100, Cmax, int(Cmax/plot_points)):
            zmin,zmax = self.theta_interval(C, a, b, max_points=max_points)
            vmin.append((C, zmin))
            vmax.append((C, zmax))
            if verbose: 
                print(C, zmin, zmax)
        return vmin, vmax
    
    def plot_theta(self, Cmax, clr=(0,0,0), show_error_bound=True, *args, **kwds):
        self.sorted_aplist(Cmax)
        if show_error_bound:
            vmin, vmax = self.compute_theta_interval(Cmax, *args, **kwds)
        v = self.compute_theta(Cmax, *args, **kwds)
        grey = (0.7,0.7,0.7)
        theta_line = line(v,rgbcolor=clr, ymin=0,ymax=1.2, legend_label='$\\theta(C)$') 
        theta_points = point(v,rgbcolor=clr)
        if show_error_bound:
            error_bound = line(vmin,rgbcolor=grey, legend_label='numerical integration error bound') + line(vmax,rgbcolor=grey)
        reference_line = line([(0,1),(Cmax,1)], rgbcolor=(1,0,0), legend_label='$\\theta=1.0$')
        
        p = theta_line + theta_points + reference_line
        if show_error_bound:
            p = p + error_bound

        p.axes_labels(['$C$ ', '$\\theta$'])
        p.legend(True)
        
        return p
        
    def histogram(self, Cmax, num_bins):
        '''
        modified with matplotlib plotting
        '''
        v = self.normalized_aplist(Cmax)
        print("num_bins:", num_bins)
        n, bins, patches = plt.hist(v, num_bins, density=True)
        angle = np.linspace(0, np.pi, 150)
        radius = 1
        x = radius * np.cos(angle) 
        y = radius * np.sin(angle) * 2 / np.pi #stretching so that area = 1
        plt.plot(x, y)
        plt.xlabel('$a_p/\\sqrt{p}$')
        plt.ylabel('Frequency')
        plt.title('Sato-Tate-Distribution for E = {}'.format(self._E))
        plt.xlim(-1,1)
        plt.grid(True)
        plt.show()
