if __name__=="__main__":

    import pickle
    from aerofem.models.aerodynamic.lifting_line import LLGalerkin
    from aerofem.utils.data_utils import get_param, save_object, load_object

    problem = LLGalerkin()

    span = get_param('example_wing', 'span')
    stations = get_param('example_wing', 'stations')
    chords = get_param('example_wing', 'chords')
    cl_alphas = get_param('example_wing', 'cl_alphas')
    alpha_l0s = get_param('example_wing', 'alpha_l0s')
    thetas = get_param('example_wing', 'thetas')

    problem.create_wing_from_sections(stations, chords, cl_alphas, alpha_l0s, thetas)

    mesh_type = get_param('example_mesh', 'mesh_type')
    Nelem = get_param('example_mesh', 'Nelem')
    r = get_param('example_mesh', 'r')
    elem_type = get_param('example_mesh', 'elem_type')

    problem.create_mesh(mesh_type, Nelem, elem_type, r)

    problem.set_flycond(70, 5)

    problem.assembly()

    problem.solve()

    problem.compute_coeff()

    save_object(problem, 'example')

    problem_loaded = load_object('projects','example')