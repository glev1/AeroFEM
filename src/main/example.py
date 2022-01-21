import numpy as np
from src.models.aerodynamic.lifting_line import LLGalerkin

span = 12

vec_Nelem = np.array([70])
vec_r = np.array([0.1])
N = vec_Nelem[0]+1
stations = -(span/2)*np.cos(np.linspace(0,np.pi,N))
chords = 1 * np.ones(np.size(stations))
area = span*chords[0]
cl_alpha = (2*np.pi)* np.ones(np.size(stations))
S = span*(chords[0] + chords[-1])
rho = 1.225

#FLAP
alpha_l0_sf = 0*np.pi/180
alpha_l0_cf = 0*np.pi/180
alpha_l0 = alpha_l0_sf * ((stations < -span/4)*1 + (stations > span/4)*1) + alpha_l0_cf * ((stations > -span/4)*1 + (stations < span/4)*1) - alpha_l0_cf

theta = 0.0 * np.ones(np.size(stations))


mesh_type = 'r'
elem_type = 'LL2'

problem = LLGalerkin()

problem.create_wing_from_sections(stations, chords, cl_alpha, alpha_l0, theta)
