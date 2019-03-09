"""Single Atom Image Analysis
Stefan Spence 26/02/19

 - watch the image_read_path directory for new images
 - save the new image with label into a dated subdirectory under image_storage_path
 - find the position of the atom in the image and integrate the signal
 - determine atom presence by comparison with a threshold count
 - plot a histogram of signal counts, which defines the threshold
 - save file references to easily re-analyse data

Assume that the first peak in the histogram is the background, 
then that the frequencies drop to zero and further peaks can be distinguished
after this point

Note that the images supplied should be an ROI so that we can integrate over them.
Assuming that image files are ASCII

Use Qt to send the signal from the watchdog to 
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
from PyQt4.QtGui import (QApplication, QPushButton, QWidget, QLabel,
        QGridLayout, QSizePolicy, QMainWindow, QMessageBox, QLineEdit, 
        QDoubleValidator) 
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
            
            # copy file with labeling: [species]_[date]_[Dexter file #]
            new_file_name = self.image_storage_path+r'\Species_'+self.date+'_'+self.dfn+'.'+event.src_path.split(".")[-1]
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
                self.image_storage_path = row.split('--')[-1] # where image files are saved
            elif "log file path" in row:
                self.log_file_path = row.split('--')[-1]      # where dat files of saved data and log files are saved
            elif "dexter sync file" in row:
                self.dexter_sync_file_name = row.split('--')[-1]   # where the txt from Dexter with the latest file # is saved
            elif "image read path" in row:
                self.image_read_path = row.split('--')[-1]    # where camera images are stored just after acquisition
        
        print("Image storage path: ", self.image_storage_path)
        print("Log file path: ", self.log_file_path)
        print("Dexter sync file: ", self.dexter_sync_file_name)
        print("Directory to watch for new files: ", self.image_read_path)
        
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
        
    def run(self):
        pass
        
    def save_config(self):
        with open('./config.dat', 'w+') as config_file:
            config_file.write('// list of required directories for SAIA\n')
            config_file.write('image storage path\t\t--'+self.image_storage_path+'\n')
            config_file.write('log file path\t\t--'+self.log_file_path+'\n')
            config_file.write('dexter sync file\t\t--'+self.dexter_sync_file_name+'\n')
            config_file.write('image read path\t\t--'+self.image_read_path+'\n')
            
            
####    ####    ####    ####
        
# convert an image into its pixel counts to put into a histogram
class image_handler:
    """load an ROI image centred on the atom, integrate the counts,
    then compare to the threshold
    for speed, make an array of counts with length n. If the number of images
    analysed exceeds (n-10) of this length then append another n"""
    def __init__(self):
        self.max_count = 2**16          # max count expected from the image
        self.n = 10000                  # length of array for storing counts
        self.counts = np.zeros(self.n)  # integrated counts from atom
        self.thresh = 1                 # initial threshold for atom detection
        self.atom = np.zeros(self.n)    # deduce presence of an atom by comparison with threshold
        # file label list length < integrated counts so that no data is lost when extra spaces have to be appended
        self.files  = [None]*(self.n-10)      # labels of files. 
        self.im_num = 0                 # number of images processed
        self.im_vals = np.array([])     # the data from the last image is accessible to an image_handler instance
        self.bin_array = []             # if bins for the histogram are supplied, plotting can be faster
        
    def process(self, im_name):
        try:
            self.add_count(im_name)
            
        except IndexError:
            self.counts = np.append(self.counts, np.zeros(self.n))
            self.files += [None]*self.n
            
            self.add_count(im_name)
            
    def add_count(self, im_name):
        """Fill the next index of the counts and files arrays"""
        # im_vals = np.array(Image.open(im_name)) # for bmp images
        self.im_vals = np.genfromtxt(im_name, delimiter=',')[:,:-1]  # ASCII image file: the last column is empty, comes out as NaN
        
        # sum the counts in the image (which should already be an ROI)
        self.counts[self.im_num] = np.sum(self.im_vals)
        
        # naming convention: [Species]_[date]_[Dexter file #]
        self.files[self.im_num] = im_name.split("_")[-1].split(".")[0]
        self.im_num += 1
        
        # xc, yc = np.where(im_vals == np.max(im_vals))
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
        
    def save_state(self, save_file_name):
        """Save the processed data to csv"""
        self.histogram() # update the threshold estimate
         
        # atom is present if the counts are above threshold
        self.atom[:self.im_num] = self.counts[:self.im_num] // self.thresh 
        
        out_arr = np.array((self.files[:self.im_num], 
                    self.counts[:self.im_num], self.atom[:self.im_num])).T
                    
        np.savetxt(save_file_name, out_arr, fmt='%s', delimiter=',',
                header='File, Counts, Atom Detected (threshold=%s)'%int(self.thresh))
            
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
        self.dw_status_label = QLabel('Stopped', self)
        grid.addWidget(self.dw_status_label, 7,2, 1,1)

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        grid.addWidget(self.recent_label, 7,3, 1,4)

        # main subplot of histogram
        self.hist_canvas = pg.PlotWidget()
        self.hist_canvas.setTitle("Histogram of CCD counts")
        grid.addWidget(self.hist_canvas, 1,0, 5,8)  # plot spans 5 rows/columns
        
        # adjustable parameters: min/max counts, number of bins
        bins_label = QLabel('Histogram Binning:', self)
        grid.addWidget(bins_label, 0,0, 1,1)
        
        # make sure user input is float:
        double_validator = QDoubleValidator()
        
        # min counts:
        min_counts_label = QLabel('Min. Counts: ', self)
        grid.addWidget(min_counts_label, 0,1, 1,1)
        self.min_counts_edit = QLineEdit(self)
        grid.addWidget(self.min_counts_edit, 0,2, 1,1)
        self.min_counts_edit.textChanged[str].connect(self.bins_text_edit)
        self.min_counts_edit.setValidator(double_validator)
        
        # max counts:
        max_counts_label = QLabel('Max. Counts: ', self)
        grid.addWidget(max_counts_label, 0,3, 1,1)
        self.max_counts_edit = QLineEdit(self)
        grid.addWidget(self.max_counts_edit, 0,4, 1,1)
        self.max_counts_edit.textChanged[str].connect(self.bins_text_edit)
        self.max_counts_edit.setValidator(double_validator)
        
        # number of bins
        num_bins_label = QLabel('# Bins: ', self)
        grid.addWidget(num_bins_label, 0,5, 1,1)
        self.num_bins_edit = QLineEdit(self)
        grid.addWidget(self.num_bins_edit, 0,6, 1,1)
        self.num_bins_edit.textChanged[str].connect(self.bins_text_edit)
        self.num_bins_edit.setValidator(double_validator)
        
        self.bins_toggle = QPushButton('Manual Binning', self)
        self.bins_toggle.setCheckable(True)
        self.bins_toggle.clicked[bool].connect(self.set_bins)
        grid.addWidget(self.bins_toggle, 0,7, 1,1)
        
        # display the last image taken
        # self.im_canvas = pg.ImageView()
        # self.im_canvas.show()

        self.setGeometry(150, 150, 700, 600)
        self.setWindowTitle('Single Atom Image Analyser')

    def init_DW(self):
        """Ask the user if they want to start the dir watcher or not"""
        reply = QMessageBox.question(self, 'Initiate the Directory Watcher',
            "Start the directory watcher?", QMessageBox.Yes | QMessageBox.No, 
                                    QMessageBox.No)

        if reply == QMessageBox.Yes:
            self.reset_DW()
            
    def reset_DW(self):
        """Initiate the dir watcher (restarts a new instance if there is already
        one running, since it crashes if there is an exception."""
        # if not self.dir_watcher: 
        self.dir_watcher = dir_watcher()
        self.dir_watcher.event_handler.event_path.connect(self.update_plot)
        self.dw_status_label.setText("Running")
        
    def bins_text_edit(self, text):
        """Update the histogram bins every time a text edit is made by the user
        to one of the line edit widgets"""
        if self.bins_toggle.isChecked():
            new_vals = [self.min_counts_edit.text(),
                            self.max_counts_edit.text(), self.num_bins_edit.text()]
                            
            # if the line edit widget is empty, take an estimate from histogram values
            if new_vals[0] == '' and self.image_handler.im_num > 0:
                new_vals[0] = min(self.image_handler.counts[:self.image_handler.im_num])
            if new_vals[1] == '' and self.image_handler.im_num > 0:
                new_vals[1] = max(self.image_handler.counts[:self.image_handler.im_num])
            if new_vals[2] == '' and self.image_handler.im_num > 0:
                new_vals[2] = 20 + self.image_handler.im_num // 20
            if self.image_handler.im_num == 0:
                new_vals = [0, 1, 10]
            min_bin, max_bin, num_bins = [float(x) for x in new_vals]
            
            # set the new values for the bins of the image handler
            self.image_handler.bin_array = np.linspace(min_bin, max_bin, num_bins)
            # update the plot
            self.plot_current_hist()
    
    def set_bins(self, toggle):
        """If the toggle is True, use automatic histogram binning.
        If the toggle is False, read in values from the line edit widgets."""
        if toggle:
            self.bins_text_edit('reset')            
        else:
            self.image_handler.bin_array = []
            self.plot_current_hist()
        
    def plot_current_hist(self):
        """Reset the plot to show the current data stored in the image handler"""
        # update the histogram and threshold estimate
        bins, occ, thresh = self.image_handler.histogram()
        
        self.hist_canvas.clear()
        self.hist_canvas.plot(bins, occ, stepMode=True,
                                fillLevel=0, brush = (250,250,250,250)) # histogram
        self.hist_canvas.plot([thresh]*2, [0, max(occ)], pen=1) # threshold line
        
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
            
            
if __name__ == "__main__":
    # if running in IPython then creating an app instance isn't necessary...
    # app = QApplication(sys.argv)  
    main_win = main_window()
    main_win.show()
    # sys.exit(app.exec_())