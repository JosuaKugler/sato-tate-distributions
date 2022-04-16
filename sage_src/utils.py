from math import asin, sqrt

def Xab(a,b):
    """
    return a function X(T)
    where X(T) is the area under the arc of the semicircle from a to T
    """
    bb = (asin(b)/2.0 + b*sqrt(1.0-b**2.0)/2.0)
    aa = (asin(a)/2.0 + a*sqrt(1.0-a**2.0)/2.0)
    def X(T):
        return (asin(T)/2.0 + T*sqrt(1.0-T**2.0)/2.0 - aa)/(bb - aa)
    return X