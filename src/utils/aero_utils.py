import numpy as np

class Wing():
    def __init__():
        return
        
    def create_wing_from_sections(self, stations, chords, cl_alpha, alpha_l0, theta,
                 span = False, aspect_ratio = False, area = False,
                 taper_ratio = False):

         self.stations = stations
         self.span = stations[-1] - stations[0]
         self.chords = chords
         self.area = np.trapz(chords,stations)
         self.cl_alpha = cl_alpha
         self.alpha_l0 = alpha_l0
         self.theta = theta
