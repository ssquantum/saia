"""Single Atom Image Analysis
Stefan Spence 15/04/19

a class to collect histogram statistics

"""
import numpy as np

class histo_handler:
    """Append histogram statistics to a list. The values are:
    user variable, loading probability, bg count, bg width, 
    signal count, signal width, separation, threshold, images processed
    """
    def __init__(self, atom_index=0, atom_symbol='Cs '):
        self.headers = np.array(['User variable', 'Number of images processed',
            'No atoms' , 'Single atom', 'Both atoms',
            'Loading probability', 'Error in loading probability',
            'Background peak count', 'Background peak width', 
            'Signal peak count', 'Signal peak width', 'Separation', 
            'Fidelity', 'Error in fidelity', 'Threshold'])
        self.vals     = np.array([]) # the variables are in the columns - [:,i]
        self.xvals    = [] # variables to plot on the x axis
        self.yvals    = [] # variables to plot on the y axis
        self.i        = atom_index  # indicates the index of this handler in the list
        self.X        = atom_symbol # the name of the atom that this handler deals with

    def load_from_log(self, fname):
        """load data from a log file, the first column is histogram number"""
        self.vals = np.loadtxt(fname, skiprows=3, delimiter=',')[:,1:]