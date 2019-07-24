"""Single Atom Image Analysis
Stefan Spence 15/04/19

a class to collect histogram statistics

"""
import numpy as np

class histo_handler:
    """Append histogram statistics to a list. These are defined in
    a dictionary so that they can each be individually managed.
    """
    def __init__(self):
        # histogram statistics and variables for plotting:
        self.stats_dict = {'Hist #':np.array([], dtype=int),
        'Counts above : below threshold':np.array([], dtype=str),
        'User variable':np.array([]),
        'Number of images processed':np.array([], dtype=int), 
        'Loading probability':np.array([]), 
        'Error in Loading probability':np.array([]),
        'Background peak count':np.array([], dtype=int), 
        'Background peak Poissonian width':np.array([], dtype=int), 
        'Background peak width':np.array([], dtype=int), 
        'Background peak error':np.array([]), 
        'Signal peak count':np.array([], dtype=int), 
        'Signal peak Poissonian width':np.array([], dtype=int),
        'Signal peak width':np.array([], dtype=int), 
        'Signal peak error':np.array([]),
        'Separation':np.array([]), 
        'Fidelity':np.array([]), 
        'Error in Fidelity':np.array([]), 
        'Threshold':np.array([])}
        # variables that won't be saved for plotting:
        self.temp_vals = {key:0 for key in self.stats_dict.keys()}
        self.xvals    = [] # variables to plot on the x axis
        self.yvals    = [] # variables to plot on the y axis
        
    def load_from_log(self, fname):
        """load data from a log file. Expect the first 3 rows to be comments.
        The 3rd row gives the column headings."""
        header=''
        with open(fname, 'r') as f:
            rows = f.read().split('\n')
        rows = list(filter(None, rows)) # get rid of empty row, usually from \n at end of file

        # get headers
        try:
            header = rows[2]
        except IndexError:
            print('Load from log warning: Invalid log file. Data was not loaded.')
            return 0
        
        # remove comments, retain compatability with old column heading
        header = header.replace('#', '').replace('loading', 'Loading'
                                        ).replace('fidelity', 'Fidelity')
        # make list
        header = np.array(header.replace('Histogram', 'Hist #').split(', '))

        # get data
        if np.size(rows) < 4:
            return 0 # no data to be loaded
        data = np.array([rows[i+3].split(',') for i in range(len(rows)-3)])
        if np.size(data) < np.size(header):
            return 0 # insufficient to be loaded
        
        n = len(data[:,0]) # number of points on the plot
        for key in self.stats_dict.keys():
            index = np.where(header == key)[0]
            if np.size(index): # if the key is in the header
                self.stats_dict[key] = np.array(data[:,index], dtype=self.stats_dict[key].dtype).reshape(n)
            else: # load an empty array
                self.stats_dict[key] = np.zeros(n, dtype=self.stats_dict[key].dtype)
        
        return 1 # success