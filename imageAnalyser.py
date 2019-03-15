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
from PyQt4.QtCore import QThread, pyqtSignal
from PyQt4.QtGui import (QApplication, QPushButton, QWidget, QLabel, QAction,
        QGridLayout, QMainWindow, QMessageBox, QLineEdit, 
        QDoubleValidator, QComboBox) 
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
            # get Dexter file number       ------------ but could we just sync at the start and then increment with every new file?
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
    
    @staticmethod 
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
        self.pic_size = 64              # number of pixels in an image
        self.thresh = 1                 # initial threshold for atom detection
        self.atom = np.zeros(self.n)    # deduce presence of an atom by comparison with threshold
        # file label list length < integrated counts so that no data is lost when extra spaces have to be appended
        self.files = np.array([None]*(self.n-10)) # labels of files. 
        self.im_num = 0                 # number of images processed
        self.im_vals = np.array([])     # the data from the last image is accessible to an image_handler instance
        self.bin_array = []             # if bins for the histogram are supplied, plotting can be faster
        
    def set_pic_size(self, im_name):
        """Set the pic size by looking at the number of columns in a file"""
        im_vals = np.genfromtxt(im_name, delimeter=self.delim)
        self.pic_size = int(np.size(im_vals[0]) - 1) # the first column of ASCII image is row #
        
    def load_full_im(self, im_name):
        """return an array with the values of the image"""
        #np.array(Image.open(im_name)) # for bmp images
        #np.genfromtxt(im_name, delimiter=self.delim)[:,1:] # first column gives column number
        return np.loadtxt(im_name, delimiter=self.delim,
                              usecols=range(1,self.pic_size-1))
        
    def process(self, im_name):
        """Get the data from an image """
        try:
            self.add_count(im_name)
            
        except IndexError: # this is a bad exception - the error might be from a bad ROI rather than reaching the end of the arrays
            self.counts = np.append(self.counts, np.zeros(self.n))
            self.mean_count = np.append(self.counts, np.zeros(self.n)) 
            self.std_count = np.append(self.counts, np.zeros(self.n))
            self.xc_list = np.append(self.xc_list, np.zeros(self.n))
            self.yc_list = np.append(self.yc_list, np.zeros(self.n))
            self.atom = np.append(self.atom, np.zeros(self.n))
            self.files = np.append(self.files, np.array([None]*self.n))
            
            self.add_count(im_name)
            
    def add_count(self, im_name):
        """Fill the next index of the counts and files arrays"""
        # could speed up by removing this if statement by having two different functions:
        # then the bool toggle to use the ROI changes which function is used.
        if self.roi_size > 0:
            self.im_vals = np.genfromtxt(im_name, delimiter=self.delim)[self.yc-self.roi_size//2:
            self.yc+self.roi_size//2, self.xc-self.roi_size//2:self.xc+self.roi_size//2]  # array of pixel values in ROI
        else:
            # note that ASCII file formats are variable... might have the last column is empty, comes out as NaN: [:,:-1]
            # might fail if trying to go really fast because the file hasn't been filled with data yet
            self.im_vals = self.load_full_im(im_name)
            if self.im_vals.size:
                pass # checking that the array isn't empty
            else:
                print("File was empty, waiting 0.01s and trying again")
                time.sleep(0.01)
                self.im_vals = self.load_full_im(im_name)
        
        # sum the counts in the image (which should already be an ROI)
        # self.counts[self.im_num] = np.sum(self.im_vals)
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
            
    def histogram(self):
        """Plot a histogram of the photon counts, determine a threshold for 
        single atom presence.
        """
        if np.size(self.bin_array) > 0: 
            occ, bins = np.histogram(self.counts[:self.im_num], self.bin_array) # fixed bins. 
        else:
            occ, bins = np.histogram(self.counts[:self.im_num]) # no bins provided, do automatic binning
            # could also have bins=20+self.im_num//10) # The number of bins increases with the number of image files processed.
        
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

        self.thresh = 0.5*(x1 + x2)  # between the two peaks
        
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
        # the first row is the header
        self.files = np.concatenate(data[1:,0], np.array([None]*self.n))
        self.counts = np.concatenate(data[1:,1], np.zeros(self.n))
        self.atom = np.concatenate(data[1:,2], np.zeros(self.n))
        self.xc_list = np.concatenate(data[1:,3], np.zeros(self.n))
        self.yc_list = np.concatenate(data[1:,4], np.zeros(self.n))
        self.mean_count = np.concatenate(data[1:,5], np.zeros(self.n)) 
        self.std_count = np.concatenate(data[1:,6], np.zeros(self.n))

        
    def save_state(self, save_file_name):
        """Save the processed data to csv"""
        self.histogram() # update the threshold estimate
         
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

        # button to initiate dir watcher
        dw_init_button = QPushButton('Initiate dir watcher', self)
        dw_init_button.clicked.connect(self.reset_DW) # function to start dir watcher
        dw_init_button.resize(dw_init_button.sizeHint())
        grid.addWidget(dw_init_button, 7,0, 1,2)

        # label to show status of dir watcher
        self.dw_status_label = QLabel('Stopped', self)  # should errors stop dir watcher???
        grid.addWidget(self.dw_status_label, 7,2, 1,1)

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        grid.addWidget(self.recent_label, 7,3, 1,4)

        # main subplot of histogram
        self.hist_canvas = pg.PlotWidget()
        self.hist_canvas.setTitle("Histogram of CCD counts")
        grid.addWidget(self.hist_canvas, 1,0, 5,7)  # plot spans 5 rows/columns
        
        # make sure user input is float:
        double_validator = QDoubleValidator()
        
        
#        # menubar allows you to load data
#        menubar = self.menuBar()
#        file_menu = menubar.addMenu('File')
#        
#        save_hist = QAction('Save histogram', self) # save current hist to csv
#        
#        load_menu = QMenu('Load Data', self)  # drop down menu for loading hist
#        load_dir = QAction('From Directory', self) # from dir of images
#        load_menu.addAction(load_dir)
#        load_csv = QAction('From csv', self) # from csv of hist data
#        load_menu.addAction(load_csv)
#        
#        file_menu.addAction(save_hist)
#        file_menu.addMenu(load_menu)
        
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
        
        # user chooses whether to use automatic or manual binning (default automatic)
        self.bins_toggle = QComboBox(self)
        self.bins_toggle.addItem('Auto Binning')
        self.bins_toggle.addItem('Manual Binning')
        # self.bins_toggle.addItem('No Update')
        self.bins_toggle.activated[str].connect(self.set_bins)
        grid.addWidget(self.bins_toggle, 0,6, 1,1)
        
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
        xc_label = QLabel('ROI x_c: ', self)
        grid.addWidget(xc_label, 7,im_grid_pos, 1,1)
        self.roi_x_edit = QLineEdit(self)
        grid.addWidget(self.roi_x_edit, 7,im_grid_pos+1, 1,1)
        self.roi_x_edit.textChanged[str].connect(self.roi_text_edit)
        self.roi_x_edit.setValidator(double_validator)
         
        # centre of ROI y position
        yc_label = QLabel('ROI y_c: ', self)
        grid.addWidget(yc_label, 7,im_grid_pos+2, 1,1)
        self.roi_y_edit = QLineEdit(self)
        grid.addWidget(self.roi_y_edit, 7,im_grid_pos+3, 1,1)
        self.roi_y_edit.textChanged[str].connect(self.roi_text_edit)
        self.roi_y_edit.setValidator(double_validator)
        
        # ROI size
        l_label = QLabel('ROI size: ', self)
        grid.addWidget(l_label, 7,im_grid_pos+4, 1,1)
        self.roi_l_edit = QLineEdit(self)
        grid.addWidget(self.roi_l_edit, 7,im_grid_pos+5, 1,1)
        self.roi_l_edit.textChanged[str].connect(self.roi_text_edit)
        self.roi_l_edit.setValidator(double_validator)
        
        # display last image if toggle is True
        self.im_canvas = pg.ImageView()
        
        grid.addWidget(self.im_canvas, 1,im_grid_pos, 5,8)
        self.im_canvas.show()

        
        self.roi_plot = self.im_canvas.getRoiPlot()
        # self.roi_plot.sigRegionChangeFinished.connect(self.update_roi) # signal emitted when user stops dragging ROI
        
        # choose main window position and dimensions: (xpos,ypos,width,height)
        self.setGeometry(100, 100, 1200, 700)
        self.setWindowTitle('Single Atom Image Analyser')
        
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
        """Initiate the dir watcher (restarts a new instance if there is already
        one running, since it crashes if there is an exception."""
        if self.dir_watcher: # check if there is a current thread
            self.dir_watcher.observer.stop()
            del self.dir_watcher    # ensure that the old thread stops
            self.dw_status_label.setText("Stopped")
            
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
                
    #### #### user input functions #### #### 
    
    def user_roi(self, pos):
        """The user drags an ROI and this updates the ROI centre and width"""
        x0, y0 = self.roi_plot.pos  # lower left corner of bounding rectangle
        xw, yw = self.roi_plot.size # widths
        l = int(0.5*(xw+yw))  # want a square ROI
        xc, yc = int(x0 + l//2), int(y0 + l//2)  # centre
        self.image_handler.set_roi([xc, yc, l])
        self.roi_x_edit.setText(str(xc))
        self.roi_y_edit.setText(str(yc))
        self.roi_l_edit.setText(str(l))
            
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
        
        if xc - l//2 < 0 or yc - l//2 < 0:
            l = 2*min([xc, yc])  # can't have the boundary go off the edge
        
        self.image_handler.set_roi(list(map(int, [xc, yc, l])))
        print("ROI: ", [xc, yc, l])
        
        # update ROI on image canvas
        self.roi_plot.setPos((xc - l//2, yc - l//2))
        self.roi_plot.setSize((l, l))
        
        
    def bins_text_edit(self, text):
        """Update the histogram bins every time a text edit is made by the user
        to one of the line edit widgets"""
        if str(self.bins_toggle.currentText()) == 'Manual Binning':
            new_vals = [self.min_counts_edit.text(),
                            self.max_counts_edit.text(), self.num_bins_edit.text()]
                            
            # if the line edit widget is empty, take an estimate from histogram values
            if new_vals[0] == '' and self.image_handler.im_num > 0:
                new_vals[0] = min(self.image_handler.counts[:self.image_handler.im_num])
            if new_vals[1] == '' and self.image_handler.im_num > 0:
                new_vals[1] = max(self.image_handler.counts[:self.image_handler.im_num])
            if new_vals[2] == '' and self.image_handler.im_num > 0:
                new_vals[2] = 20 + self.image_handler.im_num // 20
            if any([v == '' for v in new_vals]) and self.image_handler.im_num == 0:
                new_vals = [0, 1, 10]
            min_bin, max_bin, num_bins = list(map(int, new_vals))
            
            # set the new values for the bins of the image handler
            self.image_handler.bin_array = np.linspace(min_bin, max_bin, num_bins)
            # update the plot
            self.plot_current_hist()
    
    #### #### toggle functions #### #### 
        
    def set_im_show(self, toggle):
        """If the toggle is True, always update the widget with the last image"""
        if toggle:
            self.dir_watcher.event_handler.event_path.connect(self.update_im)
        else:
            try: 
                self.dir_watcher.event_handler.event_path.disconnect(self.update_im)
            except Exception: pass # if it's already been disconnected 
    
    def set_bins(self, toggle):
        """If the toggle is Auto Binning, use automatic histogram binning.
        If the toggle is Manual binning, read in values from the line edit 
        widgets.
        If the toggle is No Update, disconnect the dir watcher new event signal
        from the plot update."""
        if toggle == 'Manual Binning':
            self.bins_text_edit('reset')            
        elif toggle == 'Auto Binning':
            self.image_handler.bin_array = []
            self.plot_current_hist()
        elif toggle == 'No Update':
            pass
            # try: 
            #     self.dir_watcher.event_handler.event_path.disconnect(self.update_im)
            # except Exception: pass # if it's already been disconnected 
            
    #### #### canvas functions #### #### 
        
    def plot_current_hist(self):
        """Reset the plot to show the current data stored in the image handler"""
        # update the histogram and threshold estimate
        bins, occ, thresh = self.image_handler.histogram()
        
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
        the figure
        change bar heights for speed
        Also preferably imshow the file"""
        # add the count
        t1 = time.time()
        self.image_handler.process(event_path)
        t2 = time.time()
        self.int_time = t2 - t1
        
        # display the name of the most recent file
        self.recent_label.setText('Just processed: '+os.path.basename(event_path) )
        
        self.plot_current_hist()
        self.plot_time = time.time() - t2
        
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
            saved_file_list = os.listdir(self.dir_watcher.image_storage_path)
            measure = 0  # write a new file with an index one greater than the previous measure file
            for file_name in saved_file_list:
                if "measure" in file_name:
                    file_num = int(file_name.split('.')[0][-1])
                    if file_num > measure:
                        measure = file_num
            self.image_handler.save_state(os.path.join(self.dir_watcher.image_storage_path,
                            'measure'+str(measure)+'.csv'))
                            
            event.accept()
        
        elif reply == QMessageBox.No:
            event.accept()
            
        else:
            event.ignore()        

####    ####    ####    #### 
            
if __name__ == "__main__":
    # if running in IPython then creating an app instance isn't necessary...
    # app = QApplication(sys.argv)  
    main_win = main_window()
    main_win.show()
    # sys.exit(app.exec_())