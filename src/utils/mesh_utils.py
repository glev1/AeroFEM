import numpy as np

def r_disc(N, b, r):
    """Discretization with parameter r

    Splits span domain [-b/2, b/2] into N points using a
    geometric progression that refines tip region.

    Args:
        N (int): Number of points
        b (float): Span length
        r (float): Ratio between tip and root elements
        length

    Returns:
        array: y points positions

    """
    n = int((N-1)/2)
    if (r==1):
        y = np.linspace(-b/2, b/2, N)
    else:
        q = np.exp((np.log(r))/(n-1))
        L = (b/2)*((1-q)/(1-q**(n)))
        y2 = np.array([])

        for i in range(0,n+1):
            y1 = L*(1-q**(i))/(1-q)
            y2 = np.append(y1,y2)

        y3 = np.negative(y2)
        y4 = np.concatenate((y3,y2))
        y = np.unique(y4)
    return y
