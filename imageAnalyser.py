"""Single Atom Image Analysis
Stefan Spence 26/02/19

 - watch the image_read_path directory for new images
 - save the new image with label into a dated subdirectory under image_storage_path
 - find the position of the atom in the image and integrate the signal
 - determine atom presence by comparison with a threshold count
 - plot a histogram of signal counts, which defines the threshold
 - save file references to easily re-analyse data

Assume that there are two peaks in the histogram which are separated by a 
region of zeros

Assuming that image files are ASCII

Use Qt to send the signal from the watchdog to a real-time plot
"""
import os
import sys
import numpy as np
import pyqtgraph as pg  # not as flexible as matplotlib but works a lot better with qt
import shutil
# change directory to this file's location
os.chdir(os.path.dirname(os.path.realpath(__file__))) 
import time
# some python packages use PyQt4, some use PyQt5...
try:
    from PyQt4.QtCore import QThread, pyqtSignal
    from PyQt4.QtGui import (QApplication, QPushButton, QWidget, QLabel, QAction,
            QGridLayout, QMainWindow, QMessageBox, QLineEdit, QIcon, QFileDialog,
            QDoubleValidator, QComboBox, QMenu, QActionGroup) 
except ModuleNotFoundError:
    from PyQt5.QtCore import QThread, pyqtSignal
    from PyQt5.QtGui import (QApplication, QPushButton, QWidget, QLabel, QAction,
            QGridLayout, QMainWindow, QMessageBox, QLineEdit, QIcon, QFileDialog,
            QDoubleValidator, QComboBox, QMenu, QActionGroup) 
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


def gauss(x, A, x0, sig):
    """Gaussian centred at x0 with amplitude A and standard deviation sig.
    When A = 1 this is the normal distribution. Note that Gaussian beam
    propagation uses a 1/e^2 width wx = 2*sig."""
    return A* np.exp(-(x-x0)**2 /2. /sig**2) 

          
def est_param(x, y):
    """Generator function to estimate the parameters for a Guassian fit. 
    Assume that the data has peaks of decreasing amplitude separated by regions of zero."""
    A = np.max(y)          # amplitude is peak height
    ind = np.argmax(y)     # index of peak
    x0 = x[ind]            # centre is position of peak height
    
    # estimate half max at the closest index above the peak where it drops below A/2
    try:
        sig = abs(x0 - x[ind + np.where(y[ind:] < A/2.)[0][0]]) /1.177  # s.d. is HWHM/sqrt(2ln2) 
        # assuming there's no y-background-offset of the Gaussian
        top = ind + np.where(y[ind:] == 0)[0][0]  # upper edge of the peak is where the counts return to zero
        
    except IndexError:
        sig = abs(x0 - x[ind-1])/1.177  # probably reached the end of the list of x values
        top = len(x) - 1
    
    return A, ind, x0, sig, top
    
####    ####    ####    ####
    
# set up an event handler that is also a QObject through inheritance of QThread
class system_event_handler(FileSystemEventHandler, QThread):
    """The event handler responds to file creation events and emits the path
    to the file as a signal"""
    event_path = pyqtSignal(str)
    
    def __init__(self, image_storage_path, dexter_sync_file_name, date):
        super().__init__()
        
        self.dfn = ""    # dexter file number
        self.last_event_path = ""   # last event processed 
        self.image_storage_path = image_storage_path  # directory where we copy images to
        self.dexter_sync_file_name = dexter_sync_file_name # path to file where dexter syncs its file #
        self.date = date # today's date
        # self.init_t = time.time()      # time of initiation: use to test how long it takes to realise an event is started
        self.event_t = 0           # time taken to process the last event
        self.end_t = time.time()   # time at end of event
        self.idle_t = 0            # time between events
    
    def on_any_event(self, event):
        """On a new image being written, save the file with a synced label into the 
        image storage dir"""
        t0 = time.time()
        self.idle_t = self.end_t - t0   # duration between end of last event and start of current event
        if event.event_type == 'created' or event.event_type == 'modified': # ignoring deletion or move events
            # print(event.event_type)
            # print(event.src_path)
            # print("This function appears to run twice somteimes....")
            # get Dexter file number  
            with open(self.dexter_sync_file_name, 'r') as sync_file:
                self.dfn = str(int(sync_file.read()))
            
            # copy file with labeling: [species]_[date]_[Dexter file #] ---- this will overwrite if file already exists
            new_file_name = self.image_storage_path+r'\Cs-133_'+self.date+'_'+self.dfn+'.'+event.src_path.split(".")[-1]
            try:
                shutil.copyfile(event.src_path, new_file_name)
            except PermissionError:
                print("WARNING: added a pause because python tried to access the file before the other program had let go")
                time.sleep(0.01)
                shutil.copyfile(event.src_path, new_file_name)
            
            # os.remove(event.src_path)  # delete the old file so that we can see a new created file event
            
            self.last_event_path = new_file_name  # update last event path
            self.event_path.emit(new_file_name)  # emit signal
            
            self.end_t = time.time()       # time at end of current event
            self.event_t = self.end_t - t0 # duration of event
        
####    ####    ####    ####        
    
# setup up a watcher to detect changes in the image read directory
class dir_watcher(QThread):
    """Watches a directory to detect changes in the files present"""
    def __init__(self):
        super().__init__()
        
        # load paths used from config.dat
        self.dirs_list = self.get_dirs()  # handy list contains them all
        self.image_storage_path, self.log_file_path, self.dexter_sync_file_name, self.image_read_path = self.dirs_list
        
        # create the watchdog object
        self.observer = Observer()
        
        # get the date to be used for file labeling
        date = time.strftime("%d %b %B %Y", time.localtime()).split(" ") # day short_month long_month year
        self.date = date[0] + date[1] + date[3]  # [day][month][year]
        self.image_storage_path += r'\%s\%s\%s'%(date[3],date[2],date[0])
        
        self.event_handler = system_event_handler(self.image_storage_path, 
                                self.dexter_sync_file_name, self.date)
        
        # create image storage directory by date if it doesn't already exist
        os.makedirs(self.image_storage_path, exist_ok=True) # requies version > 3.2
        
        # initiate observer
        self.observer.schedule(self.event_handler, self.image_read_path, recursive=False)
        self.observer.start()
    
    @staticmethod # static method can be accessed without making an instance of the class
    def get_dirs():
        """Load the paths used from the config.dat file or prompt user if 
        it can't be found"""
        # load config file for directories or prompt user if first time setup
        try:
            with open('./config.dat', 'r') as config_file:
                config_data = config_file.read().split("\n")
        except FileNotFoundError:
            print("config.dat file not found. This file is required for directory references.")
            with open(input('Please supply the absolute path to config.dat\t\t'), 'r') as config_file:
                config_data = config_file.read().split("\n")
                
        for row in config_data:
            if "image storage path" in row:
                image_storage_path = row.split('--')[-1] # where image files are saved
            elif "log file path" in row:
                log_file_path = row.split('--')[-1]      # where dat files of saved data and log files are saved
            elif "dexter sync file" in row:
                dexter_sync_file_name = row.split('--')[-1]   # where the txt from Dexter with the latest file # is saved
            elif "image read path" in row:
                image_read_path = row.split('--')[-1]    # where camera images are stored just after acquisition
                
        if os.path.split(dexter_sync_file_name)[0] == image_read_path:
            print("WARNING: The directory watcher acts on all file change events, so the Dexter sync file path and image read path must be different.")
        return [image_storage_path, log_file_path, dexter_sync_file_name, image_read_path]
        
    @staticmethod
    def print_dirs(image_storage_path, log_file_path, dexter_sync_file_name, image_read_path):
        """Return a string containing information on the paths used"""
        outstr = '// list of required directories for SAIA\n'
        outstr += 'image storage path\t--'+image_storage_path+'\n'
        outstr += 'log file path\t\t--'+log_file_path+'\n'
        outstr += 'dexter sync file\t\t--'+dexter_sync_file_name+'\n'
        outstr += 'image read path\t--'+image_read_path+'\n'
        return outstr
    
    def run(self):
        pass
        
    def save_config(self):
        with open('./config.dat', 'w+') as config_file:
            config_file.write(self.print_dirs(*self.dirs_list))
            
            
####    ####    ####    ####
        
# convert an image into its pixel counts to put into a histogram
class image_handler:
    """load an ROI image centred on the atom, integrate the counts,
    then compare to the threshold
    for speed, make an array of counts with length n. If the number of images
    analysed exceeds (n-10) of this length then append another n"""
    def __init__(self):
        self.max_count = 2**16          # max count expected from the image
        self.delim = ' '                # delimieter to use when opening files
        self.n = 10000                  # length of array for storing counts
        self.counts = np.zeros(self.n)  # integrated counts from atom
        self.mean_count = np.zeros(self.n) # list of mean counts in image - estimates background 
        self.std_count = np.zeros(self.n)  # list of standard deviation of counts in image
        self.xc_list = np.zeros(self.n) # horizontal positions of max pixel
        self.yc_list = np.zeros(self.n) # vertical positions of max pixel
        self.xc = 0                     # ROI centre x position 
        self.yc = 0                     # ROI centre y position
        self.roi_size = -1              # ROI length in pixels. default -1 takes the whole image
        self.pic_size = 512             # number of pixels in an image
        self.thresh = 1                 # initial threshold for atom detection
        self.atom = np.zeros(self.n)    # deduce presence of an atom by comparison with threshold
        # file label list length < integrated counts so that no data is lost when extra spaces have to be appended
        self.files = np.array([None]*(self.n)) # labels of files. 
        self.im_num = 0                 # number of images processed
        self.im_vals = np.array([])     # the data from the last image is accessible to an image_handler instance
        self.bin_array = []             # if bins for the histogram are supplied, plotting can be faster
        
    def set_pic_size(self, im_name):
        """Set the pic size by looking at the number of columns in a file"""
        im_vals = np.genfromtxt(im_name, delimiter=self.delim)
        self.pic_size = int(np.size(im_vals[0]) - 1) # the first column of ASCII image is row number
        return self.pic_size

    def reset_arrays(self):
        """Reset all of the histogram array data to zero"""
        self.files = np.array([None]*(self.n)) # labels of files. 
        self.counts = np.zeros(self.n)  # integrated counts from atom
        self.mean_count = np.zeros(self.n) # list of mean counts in image - estimates background 
        self.std_count = np.zeros(self.n)  # list of standard deviation of counts in image
        self.xc_list = np.zeros(self.n) # horizontal positions of max pixel
        self.yc_list = np.zeros(self.n) # vertical positions of max pixel
        self.atom = np.zeros(self.n)    # deduce presence of an atom by comparison with threshold
        self.im_num = 0                 # number of images processed
        
        
    def load_full_im(self, im_name):
        """return an array with the values of the image"""
        # np.array(Image.open(im_name)) # for bmp images
        # return np.genfromtxt(im_name, delimiter=self.delim)#[:,1:] # first column gives column number
        return np.loadtxt(im_name, delimiter=self.delim,
                              usecols=range(1,self.pic_size))
        
    def process(self, im_name):
        """Get the data from an image """
        try:
            self.add_count(im_name)
            
        except IndexError: # this is a bad exception - the error might be from a bad ROI rather than reaching the end of the arrays
            # filled the array of size n so add more elements
            if self.im_num % (self.n - 10) == 0 and self.im_num > self.n / 2:
                self.counts = np.append(self.counts, np.zeros(self.n))
                self.mean_count = np.append(self.counts, np.zeros(self.n)) 
                self.std_count = np.append(self.counts, np.zeros(self.n))
                self.xc_list = np.append(self.xc_list, np.zeros(self.n))
                self.yc_list = np.append(self.yc_list, np.zeros(self.n))
                self.atom = np.append(self.atom, np.zeros(self.n))
                self.files = np.append(self.files, np.array([None]*self.n))
            self.add_count(im_name)

    def add_int_count(self, im_name):
        """Fill in the next index of the counts by summing over the ROI region and then 
        getting a counts/pixel. 
        Fill in the next index of the file, xc, yc, mean, std arrays."""
        self.im_vals = np.genfromtxt(im_name, delimiter=self.delim)[self.yc-self.roi_size//2:
            self.yc+self.roi_size//2, self.xc-self.roi_size//2:self.xc+self.roi_size//2]
        # mean counts in the image (which should already be an ROI)
        self.counts[self.im_num] = np.sum(self.im_vals) / np.size(self.im_vals)
        
        # naming convention: [Species]_[date]_[Dexter file #]
        self.files[self.im_num] = im_name.split("_")[-1].split(".")[0]
        
        # find the position of the largest pixel
        # note that the statistics of the background are undermined by always taking the max.
        try:
            self.xc_list[self.im_num], self.yc_list[self.im_num] = np.where(self.im_vals == np.max(self.im_vals))
        except ValueError: # same maximum value found in more that one position
            xcs, ycs = np.where(self.im_vals == np.max(self.im_vals))
            self.xc_list[self.im_num], self.yc_list[self.im_num] = xcs[0], ycs[0]
    
        # background statistics: mean count and standard deviation across image
        self.mean_count[self.im_num] = np.mean(self.im_vals)
        self.std_count[self.im_num] = np.std(self.im_vals, ddof=1)
        
        self.im_num += 1
            
    def add_count(self, im_name):
        """Fill the next index of the counts with the max count in a pixel.
        Fill in the next index of the file, xc, yc, mean, std arrays."""
        # could speed up by removing this if statement by having two different functions:
        # then the bool toggle to use the ROI changes which function is used.
        if self.roi_size > 0:
            self.im_vals = self.load_full_im(im_name)[self.yc-self.roi_size//2:
            self.yc+self.roi_size//2, self.xc-self.roi_size//2:self.xc+self.roi_size//2]  # array of pixel values in ROI
        else:
            # note that ASCII file formats are variable... might have the last column is empty, comes out as NaN: [:,:-1]
            # might fail if trying to go really fast because the file hasn't been filled with data yet
            # if os.stat(im_name).st_size: # check size of file in bytes (0 if unwritten) - this check only takes 0.2ms
            try:
                self.im_vals = self.load_full_im(im_name)
            except IndexError:
                print("File was empty, waiting 0.01s and trying again")
                time.sleep(0.01)
                self.im_vals = self.load_full_im(im_name)
        
        # take the max count in the image (undermines the statistics of the background)
        self.counts[self.im_num] = np.max(self.im_vals)
        
        # naming convention: [Species]_[date]_[Dexter file #]
        self.files[self.im_num] = im_name.split("_")[-1].split(".")[0]
        
        # find the position of the largest pixel
        # note that the statistics of the background are undermined by always taking the max.
        try:
            self.xc_list[self.im_num], self.yc_list[self.im_num] = np.where(self.im_vals == np.max(self.im_vals))
        except ValueError: # same maximum value found in more that one position
            xcs, ycs = np.where(self.im_vals == np.max(self.im_vals))
            self.xc_list[self.im_num], self.yc_list[self.im_num] = xcs[0], ycs[0]
    
        # background statistics: mean count and standard deviation across image
        self.mean_count[self.im_num] = np.mean(self.im_vals)
        self.std_count[self.im_num] = np.std(self.im_vals, ddof=1)
        
        self.im_num += 1
        
        # probability 6e-7 of being outside of 5 sigma
        # threshold_estimate = np.mean(im_vals) + 5*np.std(im_vals, ddof=1)
            
    def hist_and_thresh(self):
        """Make a histogram of the photon counts and determine a threshold for 
        single atom presence."""
        if np.size(self.bin_array) > 0: 
            occ, bins = np.histogram(self.counts[:self.im_num], self.bin_array) # fixed bins. 
        else:
            occ, bins = np.histogram(self.counts[:self.im_num]) # no bins provided, do automatic binning
        
        # determine a threshold separating the background peak from the single atom peak
        # there might also be peaks from more than one atom, or they might be overlapping...
        # a decent estimate of parameters is important for fitting
        # guess the threshold as the midpoint of the zero region between peaks
        # or if that doesn't exist, then the point of maximum curvature
        # peak value, index of peak, position of peak, standard deviation
        peak1, ind1, x1, sd1, edge1 = est_param(bins[:-1], occ)
        try:
            peak2, ind2, x2, sd2, edge2 = est_param(bins[edge1:-1], occ[edge1:])
        except ValueError:
            peak2, ind2, x2, sd2, edge2 = est_param(bins[:2 * ind1 - edge1 - 1], occ[:2 * ind1 - edge1])

        self.thresh = 0.5*(x1 + x2)  # midpoint between the two peaks
        
        return bins, occ, self.thresh

    def histogram(self):
        """Make a histogram of the photon counts but don't update the threshold"""
        if np.size(self.bin_array) > 0: 
            occ, bins = np.histogram(self.counts[:self.im_num], self.bin_array) # fixed bins. 
        else:
            occ, bins = np.histogram(self.counts[:self.im_num]) # no bins provided, do automatic binning
        return bins, occ, self.thresh
        
    def set_roi(self, im_name='', dimensions=[]):
        """Set the ROI for the image either by finding the position of the max 
        in the file im_name, or by taking user supplied dimensions [xc, yc, 
        roi_size]. The default is to use supplied dimensions."""
        if np.size(dimensions) != 0:
            self.xc, self.yc, self.roi_size = list(map(int, dimensions))
            return 1
            
        elif len(im_name) != 0:
            # presume the supplied image has an atom in and take the max
            # pixel's position at the centre of the ROI
            # note that the statistics of the background are undermined by always taking the max.
            im_vals = self.load_full_im(im_name)
            try:
                self.xc, self.yc = np.where(im_vals == np.max(im_vals))
            except ValueError: # same maximum value found in more that one position
                xcs, ycs = np.where(im_vals == np.max(im_vals))
                self.xc, self.yc = xcs[0], ycs[0]
            return 1
            
        else:
            # print("set_roi usage: supply im_name to get xc, yc or supply dimensions [xc, yc, l]")
            return 0 
        
    def load_from_csv(self, file_name):
        """Load back in the counts data from a stored csv file, leavning space
        in the arrays to add new data as well"""
        data = np.genfromtxt(file_name, delimiter=',', dtype=str)
        fd = data[1:,1:].astype(float) # the numerical data

        # the first row is the header
        self.files = np.concatenate((self.files[:self.im_num], data[1:,0], np.array([None]*self.n)))
        self.counts = np.concatenate((self.counts[:self.im_num], fd[:,0], np.zeros(self.n)))
        self.atom = np.concatenate((self.atom[:self.im_num], fd[:,1], np.zeros(self.n)))
        self.xc_list = np.concatenate((self.xc_list[:self.im_num], fd[:,2], np.zeros(self.n)))
        self.yc_list = np.concatenate((self.yc_list[:self.im_num], fd[:,3], np.zeros(self.n)))
        self.mean_count = np.concatenate((self.mean_count[:self.im_num], fd[:,4], np.zeros(self.n)))
        self.std_count = np.concatenate((self.std_count[:self.im_num], fd[:,5], np.zeros(self.n)))
        self.im_num += np.size(data[1:,0]) # now we have filled this many extra columns.

        
    def save_state(self, save_file_name):
        """Save the processed data to csv"""
        self.hist_and_thresh() # update the threshold estimate
         
        # atom is present if the counts are above threshold
        self.atom[:self.im_num] = self.counts[:self.im_num] // self.thresh 
        
        out_arr = np.array((self.files[:self.im_num], self.counts[:self.im_num], 
            self.atom[:self.im_num], self.xc_list[:self.im_num], 
            self.yc_list[:self.im_num], self.mean_count[:self.im_num],
            self.std_count[:self.im_num])).T
                    
        np.savetxt(save_file_name, out_arr, fmt='%s', delimiter=',',
                header='File, Counts, Atom Detected (threshold=%s), X-pos (pix), Y-pos (pix), Mean Count, s.d.'
                %int(self.thresh))
            

####    ####    ####    ####

# main GUI window contains all the widgets                
class main_window(QMainWindow):
    """Use Qt to produce the window where the histogram plot is shown
    have a simple interface for the user to control the limits of the plot 
    and number of bins in the histogram
    
    This GUI was produced with help from http://zetcode.com/gui/pyqt5/"""
    def __init__(self):
        super().__init__()
        self.dir_watcher = None  # a button will initiate the dir watcher
        self.image_handler = image_handler() # class to process images
        self.initUI()   # make the widgets
        self.init_DW()  # ask the user if they want to start the dir watcher
        self.t0 = time.time()  # time of initiation
        self.int_time = 0      # time taken to process an image
        self.plot_time = 0     # time taken to plot the graph

    def initUI(self):
        # grid layout: central main plot, params right, dir watcher status at bottom
        self.centre_widget = QWidget()
        self.setCentralWidget(self.centre_widget)
        grid = QGridLayout()
        self.centre_widget.setLayout(grid)

        # make sure user input is float:
        double_validator = QDoubleValidator()

        # menubar at top keeps things tidy
        menubar = self.menuBar()
        # file menubar allows you to save/load data
        file_menu = menubar.addMenu('File')
        
        save_hist = QAction('Save histogram', self) # save current hist to csv
        save_hist.triggered.connect(self.save_hist_data)
        
        load_menu = QMenu('Load histogram data', self)  # drop down menu for loading hist
        load_dir = QAction('From Files', self) # from image files
        load_dir.triggered.connect(self.load_from_files)
        load_menu.addAction(load_dir)
        load_csv = QAction('From csv', self) # from csv of hist data
        load_csv.triggered.connect(self.load_from_csv)
        load_menu.addAction(load_csv)

        load_im = QAction('Load Image', self) # display a loaded image
        load_im.triggered.connect(self.load_image)
        
        file_menu.addAction(load_im)
        file_menu.addAction(save_hist)
        file_menu.addMenu(load_menu)

        hist_menu =  menubar.addMenu('Histogram')

        bin_menu = QMenu('Binning', self) # drop down menu for binning options
        bin_options = QActionGroup(bin_menu)  # group together the options
        self.bin_actions = []
        for action_label in ['Automatic', 'Manual', 'No Update']:
            self.bin_actions.append(QAction(action_label, bin_menu, checkable=True, 
                            checked=action_label=='Automatic')) # default is auto
            bin_menu.addAction(self.bin_actions[-1])
            bin_options.addAction(self.bin_actions[-1])
        bin_options.setExclusive(True) # only one option checked at a time
        bin_options.triggered.connect(self.set_bins) # connect the signal
        hist_menu.addMenu(bin_menu)
            

        # button to initiate dir watcher
        self.dw_init_button = QPushButton('Initiate directory watcher', self)
        self.dw_init_button.clicked.connect(self.reset_DW) # function to start dir watcher
        self.dw_init_button.resize(self.dw_init_button.sizeHint())
        grid.addWidget(self.dw_init_button, 7,0, 1,2)

        # label to show status of dir watcher
        self.dw_status_label = QLabel('Stopped', self)  # should errors stop dir watcher???
        grid.addWidget(self.dw_status_label, 7,2, 1,1)

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        grid.addWidget(self.recent_label, 7,3, 1,4)

        # main subplot of histogram
        self.hist_canvas = pg.PlotWidget()
        self.hist_canvas.setTitle("Histogram of CCD counts")
        grid.addWidget(self.hist_canvas, 1,0, 6,8)  # plot spans 5 rows/columns
        
        # adjustable parameters: min/max counts, number of bins
        # min counts:
        min_counts_label = QLabel('Min. Counts: ', self)
        grid.addWidget(min_counts_label, 0,0, 1,1)
        self.min_counts_edit = QLineEdit(self)
        grid.addWidget(self.min_counts_edit, 0,1, 1,1)
        self.min_counts_edit.textChanged[str].connect(self.bins_text_edit)
        self.min_counts_edit.setValidator(double_validator)
        
        # max counts:
        max_counts_label = QLabel('Max. Counts: ', self)
        grid.addWidget(max_counts_label, 0,2, 1,1)
        self.max_counts_edit = QLineEdit(self)
        grid.addWidget(self.max_counts_edit, 0,3, 1,1)
        self.max_counts_edit.textChanged[str].connect(self.bins_text_edit)
        self.max_counts_edit.setValidator(double_validator)
        
        # number of bins
        num_bins_label = QLabel('# Bins: ', self)
        grid.addWidget(num_bins_label, 0,4, 1,1)
        self.num_bins_edit = QLineEdit(self)
        grid.addWidget(self.num_bins_edit, 0,5, 1,1)
        self.num_bins_edit.textChanged[str].connect(self.bins_text_edit)
        self.num_bins_edit.setValidator(double_validator)

        # user can set the threshold
        self.thresh_toggle = QPushButton('User Threshold: ', self)
        self.thresh_toggle.setCheckable(True)
        self.thresh_toggle.clicked[bool].connect(self.set_thresh)
        grid.addWidget(self.thresh_toggle, 0,6, 1,1)
        # user inputs threshold
        self.thresh_edit = QLineEdit(self)
        grid.addWidget(self.thresh_edit, 0,7, 1,1)
        self.thresh_edit.textChanged[str].connect(self.bins_text_edit)
        self.thresh_edit.setValidator(double_validator)
        
        # get user to set the image size in pixels
        size_label = QLabel('Image Size in Pixels: ', self)
        grid.addWidget(size_label, 0,8, 1,1)
        self.pic_size_edit = QLineEdit(self)
        grid.addWidget(self.pic_size_edit, 0,9, 1,1)
        self.pic_size_edit.setText(str(self.image_handler.pic_size)) # default
        self.pic_size_edit.textChanged[str].connect(self.pic_size_text_edit)
        self.pic_size_edit.setValidator(double_validator)
        
        # toggle to continuously plot images as they come in
        self.im_show_toggle = QPushButton('Auto-display last image', self)
        self.im_show_toggle.setCheckable(True)
        self.im_show_toggle.clicked[bool].connect(self.set_im_show)
        grid.addWidget(self.im_show_toggle, 0,10, 1,1)
        
        im_grid_pos = 8 # x grid position. leave enought space for the histogram
        # centre of ROI x position
        self.xc_label = QLabel('ROI x_c: ', self)
        grid.addWidget(self.xc_label, 7,im_grid_pos, 1,1)
        self.roi_x_edit = QLineEdit(self)
        grid.addWidget(self.roi_x_edit, 7,im_grid_pos+1, 1,1)
        self.roi_x_edit.textChanged[str].connect(self.roi_text_edit)
        self.roi_x_edit.setValidator(double_validator)
         
        # centre of ROI y position
        self.yc_label = QLabel('ROI y_c: ', self)
        grid.addWidget(self.yc_label, 7,im_grid_pos+2, 1,1)
        self.roi_y_edit = QLineEdit(self)
        grid.addWidget(self.roi_y_edit, 7,im_grid_pos+3, 1,1)
        self.roi_y_edit.textChanged[str].connect(self.roi_text_edit)
        self.roi_y_edit.setValidator(double_validator)
        
        # ROI size
        self.l_label = QLabel('ROI size: ', self)
        grid.addWidget(self.l_label, 7,im_grid_pos+4, 1,1)
        self.roi_l_edit = QLineEdit(self)
        grid.addWidget(self.roi_l_edit, 7,im_grid_pos+5, 1,1)
        self.roi_l_edit.textChanged[str].connect(self.roi_text_edit)
        self.roi_l_edit.setValidator(double_validator)
        
        # display last image if toggle is True
        self.im_canvas = pg.ImageView()
        grid.addWidget(self.im_canvas, 1,im_grid_pos, 6,8)
        self.roi = self.im_canvas.roi # get the ROI from the ROI plot
        self.roi.sigRegionChangeFinished.connect(self.user_roi) # signal emitted when user stops dragging ROI
        self.im_canvas.show()

        # choose main window position and dimensions: (xpos,ypos,width,height)
        self.setGeometry(100, 100, 1200, 700)
        self.setWindowTitle('Single Atom Image Analyser')
        self.setWindowIcon(QIcon('tempicon.png'))
        
    #### #### initiation functions #### #### 

    def init_DW(self):
        """Ask the user if they want to start the dir watcher or not"""
        dir_watcher_dirs = dir_watcher.get_dirs()
        text = "Loaded from config.dat:\n"
        text += dir_watcher.print_dirs(*dir_watcher_dirs)
        text += "\nStart the directory watcher with these settings?"
        reply = QMessageBox.question(self, 'Initiate the Directory Watcher',
            text, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
         
        if reply == QMessageBox.Yes:
            self.reset_DW()
            
    def reset_DW(self):
        """Initiate the dir watcher. If there is already one running, stop the 
        thread and delete the instance to ensure it doesn't run in the 
        background (which might overwrite files)."""
        if self.dir_watcher: # check if there is a current thread
            self.dir_watcher.observer.stop() # ensure that the old thread stops
            self.dir_watcher = None
            self.dw_status_label.setText("Stopped")
            self.dw_init_button.setText('Initiate directory watcher')

        else: 
            self.dir_watcher = dir_watcher()
            self.dir_watcher.event_handler.event_path.connect(self.update_plot)
            self.dw_status_label.setText("Running")
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Directory Watcher initiated with settings:\n"+
                "date\t\t\t--"+self.dir_watcher.date+"\n"+
                self.dir_watcher.print_dirs(*self.dir_watcher.dirs_list))
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            self.dw_init_button.setText('Stop directory watcher')
            self.setWindowTitle('Single Atom Image Analyser --- ' + self.dir_watcher.date)
                
    #### #### user input functions #### #### 
    
    def user_roi(self, pos):
        """The user drags an ROI and this updates the ROI centre and width"""
        x0, y0 = self.roi.pos()  # lower left corner of bounding rectangle
        xw, yw = self.roi.size() # widths
        l = int(0.5*(xw+yw))  # want a square ROI
        # note: setting the origing as bottom left but the image has origin top left
        xc, yc = int(x0 + l//2), self.image_handler.pic_size - int(y0 + l//2)  # centre
        self.image_handler.set_roi(dimensions=[xc, yc, l])
        self.xc_label.setText('ROI x_c = '+str(xc)) 
        self.yc_label.setText('ROI y_c = '+str(yc))
        self.l_label.setText('ROI size = '+str(l))
            
    def pic_size_text_edit(self, text):
        """Update the specified size of an image in pixels when the user 
        edits the text in the line edit widget"""
        self.image_handler.pic_size = int(text)
        
    
    def roi_text_edit(self, text):
        """Update the ROI position and size every time a text edit is made by
        the user to one of the line edit widgets"""
        [xc, yc, l] = [self.roi_x_edit.text(),
                            self.roi_y_edit.text(), self.roi_l_edit.text()]
        if any([v == '' for v in [xc, yc, l]]):
            xc, yc, l = 0, 0, -1 # default takes the whole image
        else:
            xc, yc, l = list(map(int, [xc, yc, l]))
        
        if (xc - l//2 < 0 or yc - l//2 < 0 
            or xc + l//2 > self.image_handler.pic_size 
            or yc + l//2 > self.image_handler.pic_size):
            l = 2*min([xc, yc])  # can't have the boundary go off the edge
        if int(l) == 0:
            l = -1 # can't have zero width
        
        self.image_handler.set_roi(dimensions=list(map(int, [xc, yc, l])))
        self.xc_label.setText('ROI x_c = '+str(xc)) 
        self.yc_label.setText('ROI y_c = '+str(yc))
        self.l_label.setText('ROI size = '+str(l))
        # update ROI on image canvas
        # note: setting the origing as bottom left but the image has origin top left
        self.roi.setPos(xc - l//2, self.image_handler.pic_size - yc - l//2)
        self.roi.setSize(l, l)
        
        
    def bins_text_edit(self, text):
        """Update the histogram bins every time a text edit is made by the user
        to one of the line edit widgets"""
        if self.bin_actions[1].isChecked(): # [auto, manual, no update]
            new_vals = [self.min_counts_edit.text(),
                            self.max_counts_edit.text(), self.num_bins_edit.text()]
                            
            # if the line edit widget is empty, take an estimate from histogram values
            if new_vals[0] == '' and self.image_handler.im_num > 0:
                new_vals[0] = min(self.image_handler.counts[:self.image_handler.im_num])
            if new_vals[1] == '' and self.image_handler.im_num > 0:
                new_vals[1] = max(self.image_handler.counts[:self.image_handler.im_num])
            elif not any([v == '' for v in new_vals[:2]]) and int(new_vals[1]) < int(new_vals[0]):
                # can't have max < min
                new_vals[1] = max(self.image_handler.counts[:self.image_handler.im_num])
            if new_vals[2] == '' and self.image_handler.im_num > 0:
                new_vals[2] = 20 + self.image_handler.im_num // 20
            if any([v == '' for v in new_vals]) and self.image_handler.im_num == 0:
                new_vals = [0, 1, 10]
            if int(new_vals[2]) < 1:
                # 0 bins causes value error
                new_vals[2] = 10
            min_bin, max_bin, num_bins = list(map(int, new_vals))
            
            # set the new values for the bins of the image handler
            self.image_handler.bin_array = np.linspace(min_bin, max_bin, num_bins)

            # set the new threshold if supplied
            if self.thresh_toggle.isChecked():
                try:
                    self.image_handler.thresh = float(self.thresh_edit.text())
                except ValueError: pass # user switched toggle before inputing text
                self.plot_current_hist(self.image_handler.histogram) # doesn't update thresh
            else:
                self.plot_current_hist(self.image_handler.hist_and_thresh) # updates thresh
            
    
    #### #### toggle functions #### #### 

    def set_thresh(self, toggle):
        """If the toggle is true, the user supplies the threshold value and it is
        kept constant using the image_handler.histogram() function. Otherwise,
        update the threshold with image_handler.hist_and_thresh()"""
        if toggle:
            try: # disconnect all slots because it might be connected several times
                self.dir_watcher.event_handler.event_path.disconnect()
            except Exception: pass # if already disconnected
            if self.dir_watcher:
                self.dir_watcher.event_handler.event_path.connect(self.update_plot_only)
        else:
            try: # disconnect all slots (including imshow...)
                self.dir_watcher.event_handler.event_path.disconnect()
            except Exception: pass # if already disconnected
            if self.dir_watcher:
                self.dir_watcher.event_handler.event_path.connect(self.update_plot)
        
    def set_im_show(self, toggle):
        """If the toggle is True, always update the widget with the last image.
        Note that disconnecting all slots means that this toggle might have to
        be reset when other buttons are pressed."""
        if toggle:
            self.dir_watcher.event_handler.event_path.connect(self.update_im)
        else:
            try: # note: it could have been connected several times... need while True: ... break
                self.dir_watcher.event_handler.event_path.disconnect(self.update_im)
            except Exception: pass # if it's already been disconnected 

    def swap_signals(self):
        """Disconnect the image_handler process signal from the dir_watcher event
        and (re)connect the update plot"""
        try: # disconnect all slots
            self.dir_watcher.event_handler.event_path.disconnect() 
        except Exception: pass
       
        if self.dir_watcher and self.thresh_toggle.isChecked():
            self.dir_watcher.event_handler.event_path.connect(self.update_plot_only)
        elif self.dir_watcher and not self.thresh_toggle.isChecked():
            self.dir_watcher.event_handler.event_path.connect(self.update_plot)
    
    def set_bins(self, action):
        """Check which of the bin action menu bar options is checked.
        If the toggle is Automatic, use automatic histogram binning.
        If the toggle is Manual, read in values from the line edit 
        widgets.
        If the toggle is No Update, disconnect the dir watcher new event signal
        from the plot update."""
        if self.bin_actions[1].isChecked(): # manual
            self.swap_signals()  # disconnect image handler, reconnect plot
            self.bins_text_edit('reset')            

        elif self.bin_actions[0].isChecked(): # automatic
            self.swap_signals()  # disconnect image handler, reconnect plot
            self.image_handler.bin_array = []
            self.plot_current_hist(self.image_handler.histogram)

        elif self.bin_actions[2].isChecked(): # no update
            try: # disconnect all slots
                self.dir_watcher.event_handler.event_path.disconnect()
            except Exception: pass # if it's already been disconnected 
            
            # just process the image and set the text of the most recent file
            if self.dir_watcher: # check that the dir watcher exists to prevent crash
                self.dir_watcher.event_handler.event_path.connect(self.image_handler.process)
                self.dir_watcher.event_handler.event_path.connect(self.recent_label.setText) # might need a better label
            
    #### #### canvas functions #### #### 
        
    def plot_current_hist(self, hist_function):
        """Reset the plot to show the current data stored in the image handler.
        hist_function is used to make the histogram and allows the toggling of
        different functions that may or may not update the threshold value."""
        # update the histogram and threshold estimate
        bins, occ, thresh = hist_function()
        
        self.hist_canvas.clear()
        self.hist_canvas.plot(bins, occ, stepMode=True,
                                fillLevel=0, brush = (250,250,250,250)) # histogram
        self.hist_canvas.plot([thresh]*2, [0, max(occ)], pen=1) # threshold line
        
    def update_im(self, event_path):
        """Receive the event path emitted from the system event handler signal
        display the image from the file in the image canvas"""
        im_vals = self.image_handler.load_full_im(event_path)
        self.im_canvas.setImage(im_vals)
        
        
    def update_plot(self, event_path):
        """Receive the event path emitted from the system event handler signal
        process the file in the event path with the image handler and update
        the figure"""
        # add the count
        t1 = time.time()
        self.image_handler.process(event_path)
        t2 = time.time()
        self.int_time = t2 - t1
        
        # display the name of the most recent file
        self.recent_label.setText('Just processed: '+os.path.basename(event_path))
        
        self.plot_current_hist(self.image_handler.hist_and_thresh) # update the displayed plot
        self.plot_time = time.time() - t2

    def update_plot_only(self, event_path):
        """Receive the event path emitted from the system event handler signal
        process the file in the event path with the image handler and update
        the figure but without changing the threshold value"""
        # add the count
        t1 = time.time()
        self.image_handler.process(event_path)
        t2 = time.time()
        self.int_time = t2 - t1
        
        # display the name of the most recent file
        self.recent_label.setText('Just processed: '+os.path.basename(event_path))
        
        self.plot_current_hist(self.image_handler.histogram) # update the displayed plot
        self.plot_time = time.time() - t2

#### #### save and load data functions #### ####

    def save_hist_data(self, trigger=None):
        """Prompt the user to give a directory to save the histogram data, then save"""
        try:
            save_file_name, _ = QFileDialog.getSaveFileName(self, 'Save File')
            self.image_handler.save_state(save_file_name)

            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("File saved to "+save_file_name)
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

        except OSError:
            pass # user cancelled - file not found

    def check_reset(self):
        """Ask the user if they would like to reset the current data stored"""
        reply = QMessageBox.question(self, 'Confirm Data Replacement',
            "Do you want to discard the current data before loading new data?", 
            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)
        
        if reply == QMessageBox.Cancel:
            return 0

        elif reply == QMessageBox.Yes:
            self.image_handler.reset_arrays() # gets rid of old data

        return 1

    def load_from_files(self, trigger=None):
        """Prompt the user to select image files to process, then sequentially process
        them and update the histogram"""
        if self.check_reset():
            try:
                self.recent_label.setText('Processing files...') # comes first otherwise not executed
                file_list, _ = QFileDialog.getOpenFileNames(self, 'Select Files')
                for file_name in file_list:
                    self.image_handler.process(file_name)
                    self.recent_label.setText('Just processed: '+os.path.basename(file_name))
            
                self.plot_current_hist(self.image_handler.histogram)

            except OSError:
                pass # user cancelled - file not found

    def load_from_csv(self, trigger=None):
        """Prompt the user to select a csv file to load histogram data from.
        It must have the specific layout that the image_handler saves in."""
        if self.check_reset():
            try:
                file_list, _ = QFileDialog.getOpenFileNames(self, 'Select A File')
                self.image_handler.load_from_csv(file_list[0])
                self.plot_current_hist(self.image_handler.histogram)

            except OSError:
                pass # user cancelled - file not found

    def load_image(self, trigger=None):
        """Prompt the user to select an image file to display"""
        try:
            file_list, _ = QFileDialog.getOpenFileNames(self, 'Select A File')
            if np.size(file_list) > 0: # avoid crash if the user cancelled
                self.update_im(file_list[0])

        except OSError:
            pass # user cancelled - file not found

    #### #### testing functions #### #### 
        
    def print_times(self, unit="s"):
        """Display the times measured for functions"""
        scale = 1
        if unit == "ms":
            scale *= 1e3
        elif unit == "us" or unit == "microseconds":
            scale *= 1e6
        else:
            unit = "s"
        if self.dir_watcher: # this is None if dir_watcher isn't initiated
            print("File copying event duration: %.4g "%(self.dir_watcher.event_handler.event_t*scale)+unit)
            print("Image processing duration: %.4g "%(self.int_time*scale)+unit)
            print("Image plotting duration: %.4g "%(self.plot_time*scale)+unit)
        else: 
            print("Initiate the directory watcher before testing timings")

    #### #### UI management functions #### #### 

    def closeEvent(self, event):
        """Prompt user to save data on closing"""
        reply = QMessageBox.question(self, 'Confirm Quit',
            "Save before you quit?", QMessageBox.Yes |
            QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)

        if reply == QMessageBox.Yes:
            # save current state:
            self.save_hist_data()
                            
            event.accept()
        
        elif reply == QMessageBox.No:
            event.accept()
            
        else:
            event.ignore()        

####    ####    ####    #### 
            
if __name__ == "__main__":
    # if running in Pylab/IPython then there may already be an app instance
    app = QApplication.instance()
    standalone = app is None # false if there is already an app instance
    if standalone: # if there isn't an instance, make one
        app = QApplication(sys.argv) 
        
    main_win = main_window()
    main_win.show()
    if standalone: # if an app instance was made, execute it
        sys.exit(app.exec_()) # when the window is closed, the python code also stops