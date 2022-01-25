import numpy as np
from aerofem.utils.aero_utils import Wing
from aerofem.utils.data_utils import save_object
from aerofem.utils.mesh_utils import Mesh, Node, LL2

class LLFourier():
    def __init__(self):
        return


class LLGalerkin(Wing, Mesh, LL2):
    def __init__(self):
        self.nodes = []
        self.elements = []
    
    def set_flycond(self, Uinf, AoA):
        self.Uinf = Uinf
        self.AoA = AoA*np.pi/180

    def create_nodes(self):
        for k in range(self.Nnodes):
            chord = np.interp(self.coords[k], self.stations, self.chords)
            cl_alpha = np.interp(self.coords[k], self.stations, self.cl_alpha)
            alpha_l0 = np.interp(self.coords[k], self.stations, self.alpha_l0)
            theta = np.interp(self.coords[k], self.stations, self.theta)
            new_node = Node(k, self.coords[k], chord, cl_alpha, alpha_l0, theta)
            self.nodes.append(new_node)
        return

    def assembly(self):
        NGDL = 1

        K1 = np.zeros((self.Nnodes*NGDL, self.Nnodes*NGDL))
        K2 = np.zeros((self.Nnodes*NGDL, self.Nnodes*NGDL))
        f = np.zeros(self.Nnodes*NGDL)

        for element1 in self.elements:
            [i1, j1, K1_elem] = element1.compute_matrix_M(self.Uinf)
            [i3, f_elem] = element1.compute_vector_f(self.AoA)

            for k in range(len(i1)):
                K1[i1[k], j1[k]] += K1_elem[k]

            for k in range(len(i3)):
                f[i3[k]] += f_elem[k]

            for element2 in self.elements:
                [i2, j2, K2_elem] = element1.compute_matrix_K(element2, self.Uinf)
                for k in range(len(i2)):
                    K2[i2[k], j2[k]] += K2_elem[k]

        self.K1, self.K2 = K1, K2
        self.K = K1 + K2


        self.K[:,0] = 0
        self.K[0,:] = 0
        self.K[0,0] = 1
        self.K[:,-1] = 0
        self.K[-1,:] = 0
        self.K[-1,-1] = 1

        f[0] = 0
        f[-1] = 0

        self.f = f
        return

    def solve(self):
        invK = np.linalg.inv(self.K)
        self.circ = np.matmul(invK, self.f)
        return self.circ

    def compute_coeff(self):
        gamma_int = 0.0
        for element in self.elements:
            gamma_int += (self.circ[element.nodes[1].label]+ self.circ[element.nodes[0].label]) * (element.h/2.0)

        self.CL = 2.0 * gamma_int/(self.area*self.Uinf)
        self.alpha_i = np.matmul(self.K2,self.circ)
        self.CDi = (2/(self.Uinf*self.area))*np.matmul(self.circ,self.alpha_i)


        print('CL_Galerkin:', self.CL)
        print('CDi_Galerkin:', self.CDi)

        return

    def save(self,
            filename: str):
        save_object(self, filename)

