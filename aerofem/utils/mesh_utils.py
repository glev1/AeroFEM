import numpy as np
import cmath
from typing import List, Union


def r_disc(N: int,
           b: float,
           r: float) -> np.array:
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

class Node():

    """
    Creates a Node object.
    """

    def __init__(self,
                label: int,
                coords: List[float],
                chords: List[float],
                cl_alpha: List[float],
                alpha_l0: List[float],
                theta: List[float]):
        # Initialize variables        
        self.label = label
        self.coords = coords
        self.chords = chords
        self.cl_alpha = cl_alpha
        self.alpha_l0 = alpha_l0
        self.theta = theta

class LL2():

    """
    Creates a LL2 element.
    """

    def __init__(self,
                nodes: List[Union[int, List[float]]],
                ngauss: int = 4):
        self.nodes, self.ngauss = nodes, ngauss

        [self.qsi, self.w] = np.polynomial.legendre.leggauss(self.ngauss)
        self.h = self.nodes[1].coords[0] - self.nodes[0].coords[0]
        self.ye = (self.nodes[1].coords[0] + self.nodes[0].coords[0])/2


    def shape_functions(self,
                        qsi: float) -> np.array:

        """Linear shape functions

        Compute shape functions N1 and N2 in some position 
        qsi of the bi-unitary domain.

        Args:
            qsi(float): Parametrised position

        Returns:
            array: Shape function values
        """
        N1 = lambda qsi: (1.0 - qsi) / 2.0
        N2 = lambda qsi: (1.0 + qsi) / 2.0
        Nmat = np.array([[N1(qsi), N2(qsi)]])
        return Nmat

    def compute_matrix_M(self, Uinf):
        n1 = self.nodes[0].label
        n2 = self.nodes[1].label

        i = np.array([n1, n1, n2, n2])
        j = np.array([n1, n2, n1, n2])

        w = self.w
        qsi = self.qsi

        Ke = np.zeros((2,2))
        coords = np.array([[self.nodes[0].coords[0]],[self.nodes[1].coords[0]]])
        chords = np.array([[self.nodes[0].chords[0]],[self.nodes[1].chords[0]]])
        cl_alphas = np.array([[self.nodes[0].cl_alpha[0]],[self.nodes[1].cl_alpha[0]]])


        for k in range(self.ngauss):
            Nmat = self.shape_functions(qsi[k])
            y = np.matmul(Nmat, coords)
            y = y[0][0]
            chord = np.matmul(Nmat, chords)
            cl_alpha = np.matmul(Nmat,cl_alphas)

        Ke += (self.h/6.0) *  np.array([[2,1],[1,2]]) * (2.0/(cl_alpha * chord * Uinf))

        values = np.array([Ke[0][0], Ke[0][1], Ke[1][0], Ke[1][1]])

        return i, j, values


    def compute_matrix_K(self, element, Uinf):
        n1i = self.nodes[0].label
        n2i = self.nodes[1].label
        n1j = element.nodes[0].label
        n2j = element.nodes[1].label

        i = np.array([n1i, n1i, n2i, n2i])
        j = np.array([n1j, n2j, n1j, n2j])
        Ke = np.zeros((2,2))

        ye = self.ye
        h_e = self.h

        yk = element.ye
        h_k = element.h
        log = lambda x: cmath.log(x)

        if (n1i == n1j): # CASO 1
            Ke[0][0] = 1.0/2.0
            Ke[0][1] = -1.0/2.0
            Ke[1][0] = -1.0/2.0
            Ke[1][1] = 1.0/2.0

        elif (n2i == n1j): # CASO 2A
            Ke[0][0] = (h_e*h_k - h_e**2*np.log(h_e) + h_k**2*np.log(h_k) + h_e**2*np.log(h_e + h_k) - h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)
            Ke[0][1] = -(h_e*h_k - h_e**2*np.log(h_e) + h_k**2*np.log(h_k) + h_e**2*np.log(h_e + h_k) - h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)
            Ke[1][0] = (-(h_e*h_k) - h_e**2*np.log(h_e) - h_k*(2*h_e + h_k)*np.log(h_k) + h_e**2*np.log(h_e + h_k) + 2*h_e*h_k*np.log(h_e + h_k) + h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)
            Ke[1][1] = (h_e*h_k + h_e**2*np.log(h_e) + h_k*(2*h_e + h_k)*np.log(h_k) - h_e**2*np.log(h_e + h_k) - 2*h_e*h_k*np.log(h_e + h_k) - h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)

        elif (n1i == n2j): # CASO 2B
            Ke[0][0] = (h_e*h_k + h_e**2*np.log(h_e) + h_k*(2*h_e + h_k)*np.log(h_k) - h_e**2*np.log(h_e + h_k) - 2*h_e*h_k*np.log(h_e + h_k) - h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)
            Ke[0][1] = (-(h_e*h_k) - h_e**2*np.log(h_e) - h_k*(2*h_e + h_k)*np.log(h_k) + h_e**2*np.log(h_e + h_k) + 2*h_e*h_k*np.log(h_e + h_k) + h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)
            Ke[1][0] = -(h_e*h_k - h_e**2*np.log(h_e) + h_k**2*np.log(h_k) + h_e**2*np.log(h_e + h_k) - h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)
            Ke[1][1] = (h_e*h_k - h_e**2*np.log(h_e) + h_k**2*np.log(h_k) + h_e**2*np.log(h_e + h_k) - h_k**2*np.log(h_e + h_k))/(2.*h_e*h_k)

        else:
            d = (ye - yk) * 2.0
            k11 = (4*h_e*h_k - (d**2 - 3*h_e**2 + 2*d*(h_e - h_k) - 2*h_e*h_k + h_k**2)*log(d - h_e - h_k) + (d + h_e - h_k)**2*log(d + h_e - h_k) +  (d**2 - 3*h_e**2 + 2*h_e*h_k + h_k**2 + 2*d*(h_e + h_k))*log(d - h_e + h_k) - (d + h_e + h_k)**2*log(d + h_e + h_k))/(8.*h_e*h_k)
            k12 = (-4*h_e*h_k + (d - h_e - h_k)*(d + 3*h_e - h_k)*log(d - h_e - h_k) - (d + h_e - h_k)**2*log(d + h_e - h_k) - (d - h_e + h_k)*(d + 3*h_e + h_k)*log(d - h_e + h_k) + (d + h_e + h_k)**2*log(d + h_e + h_k))/(8.*h_e*h_k)
            k21 = (-4*h_e*h_k + (-d + h_e + h_k)**2*log(d - h_e - h_k) - (d - 3*h_e - h_k)*(d + h_e - h_k)*log(d + h_e - h_k) - (d - h_e + h_k)**2*log(d - h_e + h_k) +  (d - 3*h_e + h_k)*(d + h_e + h_k)*log(d + h_e + h_k))/(8.*h_e*h_k)
            k22 = (4*h_e*h_k - (-d + h_e + h_k)**2*log(d - h_e - h_k) + (d - 3*h_e - h_k)*(d + h_e - h_k)*log(d + h_e - h_k) + (d - h_e + h_k)**2*log(d - h_e + h_k) -  (d - 3*h_e + h_k)*(d + h_e + h_k)*log(d + h_e + h_k))/(8.*h_e*h_k)

            Ke[0][0] = k11.real
            Ke[0][1] = k12.real
            Ke[1][0] = k21.real
            Ke[1][1] = k22.real

        values = (1.0/(4.0*np.pi*Uinf)) * np.array([Ke[0][0], Ke[0][1], Ke[1][0], Ke[1][1]])

        return i, j, values

    def compute_vector_f(self, AoA, Uinf):

        n1 = self.nodes[0].label
        n2 = self.nodes[1].label
        i = np.array([n1,n2])

        w = self.w
        qsi = self.qsi

        alphas_l0 = np.array([[self.nodes[0].alpha_l0[0]],[self.nodes[1].alpha_l0[0]]])
        thetas = np.array([[self.nodes[0].theta[0]],[self.nodes[1].theta[0]]])
        fe = np.zeros((2,1))
        for k in range(self.ngauss):
            Nmat = self.shape_functions(qsi[k])
            alpha_l0 = np.matmul(Nmat,alphas_l0)
            theta = np.matmul(Nmat,thetas)

        fe += (self.h/2.0) * np.array([[1],[1]]) * (AoA - alpha_l0 + theta)

        values = np.array([fe[0],fe[1]])
        return i, values

class Mesh():
    def __init__(self):
        return

    def create_mesh(self, mesh_type, Nelem, elem_type,r):
        self.elem_type = elem_type
        self.mesh_type = mesh_type
        self.Nelem = Nelem
        if (mesh_type == 'r'):
            if (elem_type == 'LL2'):
                self.Nnodes = Nelem+1
                self.coords = np.zeros((self.Nnodes, 2))
                self.coords[:,0] = r_disc(self.Nnodes,self.span,r)
                #self.coords[:,1] = np.zeros(self.Nnodes)
                self.coords[:,1] = 0.0*self.coords[:,0]**2
                self.create_nodes()

                for k in range(Nelem):
                    new_element = LL2((self.nodes[k], self.nodes[k+1]))
                    self.elements.append(new_element)
