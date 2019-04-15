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
        self.headers = ['user variable', 'loading probability', 
        'background peak count', 'background peak width', 
        'signal peak count', 'signal peak width', 'separation of peaks', 
        'threshold', 'images processed']
        self.vals     = [] # the variables are in the columns - [i,:]
        