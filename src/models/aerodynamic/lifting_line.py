import numpy as np

class LLFourier():
    def __init__(self):
        return


class LLGalerkin():
    def __init__(self):
        self.nodes = []
        self.elements = []

    def assembly(self):
        NGDL = 1

        K1 = np.zeros((self.Nnodes*NGDL, self.Nnodes*NGDL))
        K2 = np.zeros((self.Nnodes*NGDL, self.Nnodes*NGDL))
        f = np.zeros(self.Nnodes*NGDL)

        for element1 in self.elements:
            [i1, j1, K1_elem] = element1.K1_matrix(self.Uinf)
            [i3, f_elem] = element1.f_vector(self.AoA, self.Uinf)

            for k in range(len(i1)):
                K1[i1[k], j1[k]] += K1_elem[k]

            for k in range(len(i3)):
                f[i3[k]] += f_elem[k]

            for element2 in self.elements:
                [i2, j2, K2_elem] = element1.K2_matrix(element2, self.Uinf)
                for k in range(len(i2)):
                    K2[i2[k], j2[k]] += K2_elem[k]

        self.K1 = K1
        self.K2 = K2

        self.K = K1 + K2


        K[:,0] = 0
        K[0,:] = 0
        K[0,0] = 1
        K[:,-1] = 0
        K[-1,:] = 0
        K[-1,-1] = 1

        f[0] = 0
        f[-1] = 0

        self.f = f
        self.K = K
        return K, f, K1, K2

    def solve(self):
        invK = np.linalg.inv(self.K)
        self.circ = np.matmul(invK, self.f)
        return self.circ

    def compute_coeff(self)
        gamma_int = 0.0
        for element in self.elements:
            gamma_int += (self.circ[element.nodes[1].label]+ self.circ[element.nodes[0].label]) * (element.h/2.0)

        CL = 2.0 * gamma_int/(self.area*self.Uinf)
        alpha_i = np.matmul(K2,self.circ)
        CDi = (2/(self.Uinf*self.area))*np.matmul(self.circ,alpha_i)


        print('CL_Galerkin:', CL)
        print('CDi_Galerkin:', CDi)

        return CL, CDi
