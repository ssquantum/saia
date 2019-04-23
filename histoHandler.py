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
    def __init__(self):
        self.headers = ['User variable', 'Number of images processed', 
            'Loading probability',
            'Background peak count', 'Background peak width', 
            'Signal peak count', 'Signal peak width', 
            'Separation', 'Threshold']
        self.vals     = [] # the variables are in the columns - [:,i]
        self.xvals    = [] # variables to plot on the x axis
        self.yvals    = [] # variables to plot on the y axis