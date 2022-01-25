if __name__=="__main__":

    import numpy as np
    from aerofem.models.aerodynamic.lifting_line import LLGalerkin
    from aerofem.utils.data_utils import get_param

    span = get_param('example', 'span')
    
    Nelem = 70
    r = 0.1
    N = Nelem+1
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

    problem.create_mesh(mesh_type, Nelem, elem_type, r)

    problem.set_flycond(70, 5)

    problem.assembly()

    problem.solve()

    problem.compute_coeff()