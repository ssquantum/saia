"""Single Atom Image Analysis
Stefan Spence 26/02/19

 - watch the image_read_path directory for new images
 - save the new image with label into a dated subdirectory under image_storage_path
 - delete the original file so that a new file with the same name can be created
 - set an ROI on the image and take an integrated count from the pixels
 - determine atom presence by comparison with a threshold count
 - plot a histogram of signal counts, which defines the threshold
 - save file references to easily re-analyse data

Assume that there are two peaks in the histogram 
Assume that image files are ASCII

Use Qt to send the signal from the watchdog to a real-time plot
"""
import os
import sys
import numpy as np
from astropy.stats import binom_conf_interval
import pyqtgraph as pg    # not as flexible as matplotlib but works a lot better with qt
# change directory to this file's location
os.chdir(os.path.dirname(os.path.realpath(__file__))) 
import imageHandler as ih # process images to build up a histogram
import histoHandler as hh # collect data from histograms together
import directoryWatcher as dw # use watchdog to get file creation events
import fitCurve as fc   # custom class to get best fit parameters using curve_fit
import time
# some python packages use PyQt4, some use PyQt5...
try:
    from PyQt4.QtCore import QThread, pyqtSignal, QEvent
    from PyQt4.QtGui import (QApplication, QPushButton, QWidget, QLabel, QAction,
            QGridLayout, QMainWindow, QMessageBox, QLineEdit, QIcon, QFileDialog,
            QDoubleValidator, QIntValidator, QComboBox, QMenu, QActionGroup, 
            QTabWidget, QVBoxLayout, QFont, QInputDialog) 
except ModuleNotFoundError:
    from PyQt5.QtCore import QThread, pyqtSignal, QEvent
    from PyQt5.QtGui import (QGridLayout, QMessageBox, QLineEdit, QIcon, 
            QFileDialog, QDoubleValidator, QIntValidator, QComboBox, QMenu, 
            QActionGroup, QVBoxLayout, QFont)
    from PyQt5.QtWidgets import (QApplication, QPushButton, QWidget, QTabWidget,
        QAction, QMainWindow, QLabel, QInputDialog)
          
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
        self.atomX = ['Cs ', 'Rb ']   # term labels for atoms
        self.c = [(255,127,14), (31,119,180)] # colours to plot in 
        self.image_handler = [ih.image_handler(i, self.atomX[i]) for i in range(len(self.atomX))] # class to process images
        self.histo_handler = [hh.histo_handler(i, self.atomX[i]) for i in range(len(self.atomX))] # class to process histograms
        pg.setConfigOption('background', 'w') # set graph background default white
        pg.setConfigOption('foreground', 'k') # set graph foreground default black
        self.date = time.strftime("%d %b %B %Y", time.localtime()).split(" ") # day short_month long_month year
        self.init_UI()  # make the widgets
        self.init_DW()  # ask the user if they want to start the dir watcher
        self.init_log() # write header to the log file that collects histograms
        self.t0 = time.time()  # time of initiation
        self.int_time = 0      # time taken to process an image
        self.plot_time = 0     # time taken to plot the graph

    def init_log(self):
        """Create a directory for today's date as a subdirectory in the log file path
        then write the header to the log file path defined in config.dat"""
        dir_watcher_dirs = dw.dir_watcher.get_dirs() # static method
        # make subdirectory if it doesn't already exist
        log_file_dir = dir_watcher_dirs[1]+'\%s\%s\%s'%(self.date[3],self.date[2],self.date[0]) 
        try:
            os.makedirs(log_file_dir, exist_ok=True)
        except PermissionError:  # couldn't access the path, start a log file here
            log_file_dir = '.\%s\%s\%s'%(self.date[3],self.date[2],self.date[0])
            os.makedirs(log_file_dir, exist_ok=True)

        # log is saved in a dated subdirectory and the file name also has the date
        # make a separate log file for each atomic species:
        self.log_file_names = [os.path.join(log_file_dir, 
                   X+self.date[0]+self.date[1]+self.date[3]+'.dat')  for X in self.atomX]
        # write the header to the log file
        for file_name in self.log_file_names:
            if not os.path.isfile(file_name): # don't overwrite if it already exists
                with open(file_name, 'w+') as f:
                    f.write('//Single Atom Image Analyser Log File: collects histogram data\n')
                    f.write('include --[]\n')
                    f.write('Histogram, '+', '.join(self.histo_handler[0].headers)+'\n')
       

    def init_UI(self):
        """Create all of the widget objects required"""
        # grid layout: central main plot, params above, dir watcher status at bottom
        self.centre_widget = QWidget()
        self.tabs = QTabWidget()                   # make tabs for each main display 
        self.centre_widget.layout = QVBoxLayout()
        self.centre_widget.layout.addWidget(self.tabs)
        self.centre_widget.setLayout(self.centre_widget.layout)
        self.setCentralWidget(self.centre_widget)
        
        # make sure user input is float:
        double_validator = QDoubleValidator()
        int_validator = QIntValidator()

        # change font size
        font = QFont()
        font.setPixelSize(12)

        #### menubar at top gives options ####
        menubar = self.menuBar()

        # file menubar allows you to save/load data
        file_menu = menubar.addMenu('File')
        load_im = QAction('Load Image', self) # display a loaded image
        load_im.triggered.connect(self.load_image)
        file_menu.addAction(load_im)
        
        # histogram menu saves/loads/resets histogram and gives binning options
        hist_menu =  menubar.addMenu('Histogram')

        save_hist = QAction('Save histogram', self) # save current hist to csv
        save_hist.triggered.connect(self.save_hist_data)
        hist_menu.addAction(save_hist)

        reset_hist = QAction('Reset histogram', self) # reset hist without loading new data
        reset_hist.triggered.connect(self.check_reset)
        hist_menu.addAction(reset_hist)
        
        load_menu = QMenu('Load histogram data', self)  # drop down menu for loading hist
        load_dir = QAction('From Files', self) # from image files
        load_dir.triggered.connect(self.load_from_files)
        load_menu.addAction(load_dir)
        load_csv = QAction('From csv', self) # from csv of hist data
        load_csv.triggered.connect(self.load_from_csv)
        load_menu.addAction(load_csv)
        
        hist_menu.addMenu(load_menu)

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

        # toggle whether to fix threshold at user specified value
        self.thresh_toggle = QAction('User Threshold', self, checkable=True)
        self.thresh_toggle.triggered.connect(self.set_thresh)
        hist_menu.addAction(self.thresh_toggle)

        # load plots from log files
        varplot_menu = menubar.addMenu('Plotting')
        # load histogram data into varplot
        load_varplot = QAction('Load from log file', self) 
        load_varplot.triggered.connect(self.load_from_log)
        varplot_menu.addAction(load_varplot)
        # choose which atoms to plot
        atom_menu = QMenu('Choose Species', self) # drop down menu for choosing atom
        self.atom_varplot_toggles = {}
        for X in self.atomX:
            self.atom_varplot_toggles[X] = QAction(X, atom_menu, checkable=True,
                                                checked=True)
            atom_menu.addAction(self.atom_varplot_toggles[X])
            self.atom_varplot_toggles[X].triggered.connect(self.update_varplot_axes)
        varplot_menu.addMenu(atom_menu)

        #### tab for settings  ####
        settings_tab = QWidget()
        settings_grid = QGridLayout()
        settings_tab.setLayout(settings_grid)
        self.tabs.addTab(settings_tab, "Settings")

        # get user to set the image size in pixels
        size_label = QLabel('Image size in pixels: ', self)
        settings_grid.addWidget(size_label, 0,0, 1,1)
        self.pic_size_edit = QLineEdit(self)
        settings_grid.addWidget(self.pic_size_edit, 0,1, 1,1)
        self.pic_size_edit.setText(str(self.image_handler[0].pic_size)) # default
        self.pic_size_edit.textChanged[str].connect(self.pic_size_text_edit)
        self.pic_size_edit.setValidator(double_validator)
        self.pic_size_edit.resize(self.pic_size_edit.sizeHint())

        # get image size from loading an image
        load_im_size = QPushButton('Load size from image', self)
        load_im_size.clicked.connect(self.load_im_size) # load image size from image
        load_im_size.resize(load_im_size.sizeHint())
        settings_grid.addWidget(load_im_size, 0,3, 1,1)

        # get ROI centre from loading an image
        load_roi = QPushButton('Get ROI from image', self)
        load_roi.clicked.connect(self.load_roi) # load roi centre from image
        load_roi.resize(load_roi.sizeHint())
        settings_grid.addWidget(load_roi, 0,5, 1,1)

        # get user to set ROI for Cs and Rb
        self.roi_edits = {} # make dictionary of QLineEdit 
        self.roi_label_text = ['xc: ', 'yc: ', 'size: ']
        for i, X in enumerate(self.atomX):
            for ii in range(len(self.roi_label_text)):
                text = X + self.roi_label_text[ii]  # text for the label identifies the QLineEdit
                new_label = QLabel(text, self)
                settings_grid.addWidget(new_label, i+1,2*ii, 1,1)
                self.roi_edits[text] = QLineEdit(self)
                settings_grid.addWidget(self.roi_edits[text], i+1,2*ii+1, 1,1)
                self.roi_edits[text].setText(str(ii//2)) # default
                self.roi_edits[text].textEdited[str].connect(self.roi_text_edit)
                self.roi_edits[text].setValidator(int_validator) # only numbers

        
        # change paths used for directory watcher
        path_label_text = ['Image Storage Path: ', 'Log File Path: ', 
                'Dexter Sync File: ', 'Image Read Path: ', 'Results Path: ']
        self.path_label_edits = []
        for i in range(len(path_label_text)):
            new_label = QLabel(path_label_text[i], self)
            settings_grid.addWidget(new_label, i+4,0, 1,1)
            self.path_label_edits.append(QLineEdit(self))
            self.path_label_edits[i].setObjectName(path_label_text[i])
            settings_grid.addWidget(self.path_label_edits[i], i+4,1, 1,3)
            self.path_label_edits[i].returnPressed.connect(self.path_text_edit)
        
        # button to initiate dir watcher
        self.dw_init_button = QPushButton('Initiate directory watcher', self)
        self.dw_init_button.clicked.connect(self.reset_DW) # function to start/stop dir watcher
        self.dw_init_button.resize(self.dw_init_button.sizeHint())
        settings_grid.addWidget(self.dw_init_button, i+5,0, 1,1)

        # label to show status of dir watcher
        self.dw_status_label = QLabel('Stopped', self)  # should errors stop dir watcher???
        settings_grid.addWidget(self.dw_status_label, i+5,1, 1,1)

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        settings_grid.addWidget(self.recent_label, i+6,0, 1,4)
        

        #### tab for histogram ####
        hist_tab = QWidget()
        hist_grid = QGridLayout()
        hist_tab.setLayout(hist_grid)
        self.tabs.addTab(hist_tab, "Histogram")

        # main subplot of histogram
        self.hist_canvas = [pg.PlotWidget() for i in range(len(self.atomX))]
        for i in range(len(self.hist_canvas)):
            self.hist_canvas[i].getAxis('bottom').tickFont = font
            self.hist_canvas[i].getAxis('left').tickFont = font # not doing anything...
        hist_grid.addWidget(self.hist_canvas[0], 1,0, 1,8)  # allocate space in the grid
        hist_grid.addWidget(self.hist_canvas[1], 3,0, 1,8)
        
        # adjustable parameters: min/max counts, number of bins
        self.hist_edits = {}
        self.hist_label_text = ['Min. Counts: ', 'Max. Counts: ', '# Bins: ', 'Threshold: ']
        for i, X in enumerate(self.atomX):
            for ii in range(len(self.hist_label_text)):
                text = X + self.hist_label_text[ii]
                new_label = QLabel(text, self)
                hist_grid.addWidget(new_label, i*2,2*ii, 1,1)
                self.hist_edits[text] = QLineEdit(self)
                hist_grid.addWidget(self.hist_edits[text], i*2,2*ii+1, 1,1)
                self.hist_edits[text].textChanged[str].connect(self.bins_text_edit)
                self.hist_edits[text].setValidator(double_validator)

        
        #### tab for current histogram statistics ####
        stat_tab = QWidget()
        stat_grid = QGridLayout()
        stat_tab.setLayout(stat_grid)
        self.tabs.addTab(stat_tab, 'Histogram Statistics')

        # user variable value
        user_var_label = QLabel('Variable value: ', self)
        stat_grid.addWidget(user_var_label, 0,0, 1,1)
        self.var_edit = QLineEdit(self)
        stat_grid.addWidget(self.var_edit, 0,1, 1,2)
        self.var_edit.setText('0')  # default
        self.var_edit.setValidator(double_validator) # only numbers

        self.stat_labels = {}  # dictionary of stat labels
        # get the list of labels from the histogram handler
        label_text = ['Counts above : below threshold'] + list(self.histo_handler[0].headers[1:])
        for ii in range(1, 1+len(label_text)):
            new_label = QLabel(label_text[ii-1], self) # description
            for i, X in enumerate(self.atomX):
                stat_grid.addWidget(new_label, ii,0, 1,1)
                self.stat_labels[X + label_text[ii-1]] = QLabel('', self) # value
                stat_grid.addWidget(self.stat_labels[X + label_text[ii-1]], ii,i+1, 1,1)
            
        # update statistics
        stat_update = QPushButton('Update statistics', self)
        stat_update.clicked[bool].connect(self.update_stats)
        stat_grid.addWidget(stat_update, ii+1,0, 1,1)

        # do Gaussian/Poissonian fit - peaks and widths
        fit_update = QPushButton('Get best fits', self)
        fit_update.clicked[bool].connect(self.update_fit)
        stat_grid.addWidget(fit_update, ii+1,1, 1,1)

        # quickly add the current histogram statistics to the plot
        add_to_plot = QPushButton('Add to plot', self)
        add_to_plot.clicked[bool].connect(self.add_stats_to_plot)
        stat_grid.addWidget(add_to_plot, ii+1,2, 1,1)

        #### tab for viewing images ####
        im_tab = QWidget()
        im_grid = QGridLayout()
        im_tab.setLayout(im_grid)
        self.tabs.addTab(im_tab, 'Image')
        # display the pic size widgets on this tab as well
        im_size_label = QLabel('Image Size in Pixels: ', self)
        im_grid.addWidget(im_size_label, 0,0, 1,1)
        self.pic_size_label = QLabel('', self)
        im_grid.addWidget(self.pic_size_label, 0,1, 1,1)
        self.pic_size_label.setText(str(self.image_handler[0].pic_size)) # default

        # toggle to continuously plot images as they come in
        self.im_show_toggle = QPushButton('Auto-display last image', self)
        self.im_show_toggle.setCheckable(True)
        self.im_show_toggle.clicked[bool].connect(self.set_im_show)
        im_grid.addWidget(self.im_show_toggle, 0,2, 1,1)
        
        # centre of ROI x position
        self.roi_labels = {}
        for i, X in enumerate(self.atomX):
            im_grid_pos = 0 # starting column. 
            for ii in range(len(self.roi_label_text)):
                text = X + self.roi_label_text[ii]
                self.roi_labels[text] = QLabel(text+str(ii//2), self) # default: 0,0,1
                im_grid.addWidget(self.roi_labels[text], 9+i,im_grid_pos, 1,1)
                im_grid_pos += 2

        # display last image if toggle is True
        im_widget = pg.GraphicsLayoutWidget() # containing widget
        viewbox = im_widget.addViewBox() # plot area to display image
        self.im_canvas = pg.ImageItem() # the image
        viewbox.addItem(self.im_canvas)
        im_grid.addWidget(im_widget, 1,0, 8,8)
        # make a ROIs that the user can drag. One for Cs, one for Rb
        self.rois = np.array((pg.ROI([0,0], [1,1], snapSize=1, scaleSnap=True, 
                                translateSnap=True, pen=pg.mkPen(color=self.c[0],width=4)),
                    pg.ROI([0,1], [1,1], snapSize=1, scaleSnap=True, 
                                translateSnap=True, pen=pg.mkPen(color=self.c[1],width=4))))
        for roi in self.rois:
            roi.addScaleHandle([1,1], [0.5,0.5]) # allow user to adjust ROI size
            viewbox.addItem(roi)
            roi.setZValue(10)   # make sure the ROI is drawn above the image
            roi.sigRegionChangeFinished.connect(self.user_roi) # signal emitted when user stops dragging ROI

        # make a histogram to control the intensity scaling
        self.im_hist = pg.HistogramLUTItem()
        self.im_hist.setImageItem(self.im_canvas)
        im_widget.addItem(self.im_hist)
        # self.im_canvas.show()


        #### tab for plotting variables ####
        plot_tab = QWidget()
        plot_grid = QGridLayout()
        plot_tab.setLayout(plot_grid)
        self.tabs.addTab(plot_tab, 'Plotting')

        # main plot
        self.varplot_canvas = pg.PlotWidget()
        self.varplot_canvas.getAxis('bottom').tickFont = font
        self.varplot_canvas.getAxis('left').tickFont = font
        plot_grid.addWidget(self.varplot_canvas, 0,1, 6,8)
        
        # x and y labels
        self.plot_labels = [QComboBox(self), QComboBox(self)]
        for i in range(len(self.plot_labels)):
            self.plot_labels[i].addItems(self.histo_handler[0].headers) # add options
            # connect buttons to update functions
            self.plot_labels[i].activated[str].connect(self.update_varplot_axes)
        # position labels in grid
        plot_grid.addWidget(self.plot_labels[0], 7,4, 1,1) # bottom middle
        plot_grid.addWidget(self.plot_labels[1], 2,0, 1,1) # middle left

        # button to clear plot data (it's still saved in the log file)
        clear_varplot = QPushButton('Clear plot', self)
        clear_varplot.clicked[bool].connect(self.clear_varplot)
        plot_grid.addWidget(clear_varplot, 7,0, 1,1)

        #### choose main window position and dimensions: (xpos,ypos,width,height)
        self.setGeometry(100, 100, 850, 700)
        self.setWindowTitle('Single Atom Image Analyser')
        self.setWindowIcon(QIcon('docs/tempicon.png'))
        
    #### #### initiation functions #### #### 

    def init_DW(self):
        """Ask the user if they want to start the dir watcher or not"""
        dir_watcher_dirs = dw.dir_watcher.get_dirs() # static method
        text = "Loaded from config.dat:\n"
        text += dw.dir_watcher.print_dirs(*dir_watcher_dirs) # static method
        text += "\nStart the directory watcher with these settings?"
        reply = QMessageBox.question(self, 'Initiate the Directory Watcher',
            text, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
         
        if reply == QMessageBox.Yes:
            self.reset_DW()

    def remove_im_files(self):
        """Ask the user if they want to remove image files from the read image
        path since the dir watcher only notices file created not modified"""
        text = 'The directory watcher only notices file creation events, not modifications.\n'
        text += 'Therefore the image read path must be emptied so new files can be created.\n'
        text += '\nDelete the following files from '+self.dir_watcher.image_read_path+"?\n"
        file_list = []
        for file_name in os.listdir(self.dir_watcher.image_read_path):
            if '.asc' in file_name:
                text += "\t - " + file_name + "\n"
                file_list.append(file_name)

        reply = QMessageBox.question(self, 'Remove Initial Image files?',
            text, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            for file_name in file_list:
                os.remove(os.path.join(self.dir_watcher.image_read_path, file_name))
            
    def reset_DW(self):
        """Initiate the dir watcher. If there is already one running, stop the 
        thread and delete the instance to ensure it doesn't run in the 
        background (which might overwrite files)."""
        if self.dir_watcher: # check if there is a current thread
            self.dir_watcher.observer.stop() # ensure that the old thread stops
            self.dir_watcher = None
            self.dw_status_label.setText("Stopped")
            self.dw_init_button.setText('Initiate directory watcher') # turns on
            self.recent_label.setText('')

        else: 
            # prompt user if they want to remove image files 
            self.dir_watcher = dw.dir_watcher() # instantiate dir watcher
            self.remove_im_files() # prompt to remove image files
            self.dir_watcher.event_handler.event_path.connect(self.update_plot) # default
            self.dw_status_label.setText("Running")

            # get current date
            self.date = self.dir_watcher.date
            date_str = ' '.join([self.date[0]]+self.date[2:])
            self.init_log() # make a new log file

            msg = QMessageBox() # pop up box to confirm it's started
            msg.setIcon(QMessageBox.Information)
            msg.setText("Directory Watcher initiated with settings:\n"+
                "date\t\t\t--" + date_str + "\n"+
                self.dir_watcher.print_dirs(*[d for d in self.dir_watcher.dirs_dict.values()]))
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            self.dw_init_button.setText('Stop directory watcher') # turns off
            # display current date on window title
            self.setWindowTitle('Single Atom Image Analyser --- ' + date_str)

            # set current file paths
            for i, label in enumerate(self.dir_watcher.dirs_dict.values()):
                self.path_label_edits[i].setText(label)

                
    #### #### user input functions #### #### 

    def path_text_edit(self, text=''):
        """The user finishes editing an edit text box by pressing return or clicking
        somewhere else, then the text is sent to this function to update the dir
        watcher's paths. It is the system event handler that processes file 
        creation events, so its paths must be updated. For consistency the dir 
        watcher's paths are also updated."""
        if self.dir_watcher:
            sender = self.sender()
            
            for i, label in enumerate(self.dir_watcher.dirs_dict):
                if i == 0 and label in sender.objectName(): # image storage path
                    new_image_storage_path = sender.text() + r'\%s\%s\%s'%(self.date[3],self.date[2],self.date[0])
                    # create image storage directory by date if it doesn't already exist
                    os.makedirs(new_image_storage_path, exist_ok=True) # requies Python version > 3.2
                    self.dir_watcher.event_handler.image_storage_path = new_image_storage_path # saves files
                    self.dir_watcher.image_storage_path = new_image_storage_path # for consistency
                    
                elif i == 1 and label in sender.objectName(): # log file path
                    self.dir_watcher.log_file_path = sender.text()

                elif i == 2 and label in sender.objectName(): # dexter sync file
                    self.dir_watcher.event_handler.dexter_sync_file_name = sender.text() # checks current file
                    self.dir_watcher.dexter_sync_file_name = sender.text() # for consistency

                elif i == 3 and label in sender.objectName(): # image read path
                    for path_edit in self.path_label_edits:
                        self.dir_watcher.dirs_dict[path_edit.objectName()] = path_edit.text()
                    self.dir_watcher.save_config()  # update the config file
                    self.reset_DW()  # stop the dir watcher
                    self.init_DW()   # prompt user to restart the dir watcher with new config
                    break # the loop can't continue if the dir watcher has been deleted

                elif i == 4 and label in sender.objectName(): # results path
                    self.dir_watcher.results_path = sender.text()
            
    
    def user_roi(self, pos):
        """The user drags an ROI and this updates the ROI centre and width"""
        roi_idx = np.where(self.rois == self.sender())[0][0]  # index of which roi we changed
        roi = self.rois[roi_idx]
        x0, y0 = roi.pos()  # lower left corner of bounding rectangle
        xw, yw = roi.size() # widths
        l = int(0.5*(xw+yw))  # want a square ROI
        # note: setting the origin as bottom left but the image has origin top left
        xc, yc = int(x0 + l//2), int(y0 + l//2)  # centre
        
        new_dim = [xc, yc, l]  # new dimensions for ROI
        self.image_handler[roi_idx].set_roi(dimensions=new_dim)
        for i in range(len(new_dim)):
            text = self.atomX[roi_idx] + self.roi_label_text[i]
            self.roi_labels[text].setText(text + str(new_dim[i]))
            self.roi_edits[text].setText(str(new_dim[i]))
            
    def pic_size_text_edit(self, text):
        """Update the specified size of an image in pixels when the user 
        edits the text in the line edit widget"""
        for i in range(len(self.image_handler)):
            self.image_handler[i].pic_size = int(text)
            self.pic_size_label.setText(str(self.image_handler[i].pic_size))

    def get_atom_idx(self, dict_items, sender):
        """Find the index of the atom term symbols list where the sender object
        matches an item in the dictionary"""
        for key, item in dict_items:
            if item == sender:
                return np.where([X == key[:3] for X in self.atomX])[0][0], key
        
    def roi_text_edit(self, text):
        """Update the ROI position and size every time a text edit is made by
        the user to one of the line edit widgets"""
        # find which ROI edit was changed
        roi_idx, _ = self.get_atom_idx(self.roi_edits.items(), self.sender())
        
        # [xc, yc, l]
        new_dim = [self.roi_edits[self.atomX[roi_idx]+dim].text() for dim in self.roi_label_text]
        
        if any([v == '' for v in new_dim]):
            new_dim = [0, 0, 1] # default takes the top left pixel
        else:
            new_dim = list(map(int, new_dim)) # crashes if the user inputs float
        
        if (new_dim[0] - new_dim[2]//2 < 0 or new_dim[1] - new_dim[2]//2 < 0 
            or new_dim[0] + new_dim[2]//2 > self.image_handler[roi_idx].pic_size 
            or new_dim[1] + new_dim[2]//2 > self.image_handler[roi_idx].pic_size):
            new_dim[2] = 2*min([new_dim[0], new_dim[1]])  # can't have the boundary go off the edge
        if int(new_dim[2]) == 0:
            new_dim[2] = 1 # can't have zero width
        
        self.image_handler[roi_idx].set_roi(dimensions=list(map(int, new_dim)))
        for i in range(len(new_dim)):
            text = self.atomX[roi_idx] + self.roi_label_text[i]
            self.roi_labels[text].setText(text + str(new_dim[i]))
        # update ROI on image canvas
        # note: setting the origin as top left because image is inverted
        self.rois[roi_idx].setPos(new_dim[0] - new_dim[2]//2, new_dim[1] - new_dim[2]//2) # xc-l//2, yc-l//2
        self.rois[roi_idx].setSize(new_dim[2], new_dim[2]) # l, l
        
        
    def bins_text_edit(self, text):
        """Update the histogram bins every time a text edit is made by the user
        to one of the line edit widgets"""
        if self.bin_actions[1].isChecked(): # [auto, manual, no update]
            # 0: Cs, 1: Rb
            idxkey = self.get_atom_idx(self.hist_edits.items(), self.sender())
            if idxkey:  # update just one atom
                idxs = [idxkey[0]]
                keys = [idxkey[1]]
            else:       # update both atoms
                idxs = range(len(self.atomX))
                keys = self.atomX
            for idx, key in list(map(list, zip(*(idxs, keys)))): # transpose iterables
                # min counts, max counts, # bins
                new_vals = [self.hist_edits[key[:3]+self.hist_label_text[0]].text(),
                            self.hist_edits[key[:3]+self.hist_label_text[1]].text(), 
                            self.hist_edits[key[:3]+self.hist_label_text[2]].text()]
                                
                # if the line edit widget is empty, take an estimate from histogram values
                if new_vals[0] == '' and self.image_handler[idx].im_num > 0:
                    new_vals[0] = min(self.image_handler[idx].counts[:self.image_handler[idx].im_num])
                if new_vals[1] == '' and self.image_handler[idx].im_num > 0:
                    new_vals[1] = max(self.image_handler[idx].counts[:self.image_handler[idx].im_num])
                elif not any([v == '' for v in new_vals[:2]]) and int(new_vals[1]) < int(new_vals[0]):
                    # can't have max < min
                    new_vals[1] = max(self.image_handler[idx].counts[:self.image_handler[idx].im_num])
                if new_vals[2] == '' and self.image_handler[idx].im_num > 0:
                    new_vals[2] = 20 + self.image_handler[idx].im_num // 20
                if any([v == '' for v in new_vals]) and self.image_handler[idx].im_num == 0:
                    new_vals = [0, 1, 10]
                if int(new_vals[2]) < 2:
                    # 0 bins causes value error
                    new_vals[2] = 10
                min_bin, max_bin, num_bins = list(map(int, new_vals))
                
                # set the new values for the bins of the image handler
                self.image_handler[idx].bin_array = np.linspace(min_bin, max_bin, num_bins)

                # set the new threshold if supplied
                if self.thresh_toggle.isChecked():
                    try:
                        # get threshold from QLineEdit
                        self.image_handler[idx].thresh = float(
                                self.hist_edits[key[:3]+self.hist_label_text[3]].text())
                        # show the new threshold in the label
                        self.stat_labels[key[:3]+'Threshold'].setText(
                                str(int(self.image_handler[idx].thresh)))
                    except ValueError: pass # user switched toggle before inputing text
                    self.plot_current_hist([x.histogram for x in self.image_handler]) # doesn't update thresh
                else:
                    self.plot_current_hist([x.hist_and_thresh for x in self.image_handler]) # updates thresh
            
    
    #### #### toggle functions #### #### 

    def update_stats(self, toggle=True):
        """Update the statistics from the current histograms
        image_handler uses a peak finding algorithm to get the peak positions and widths
        from these a threshold can be set. If the user thresh_toggle is checked then the
        threshold will not be updated."""
        new_stats = []
        for im_han in self.image_handler:
            if im_han.im_num > 0: # only update if a histogram exists
                if self.thresh_toggle.isChecked(): # using manual threshold
                    self.plot_current_hist([im_han.histogram]) # update hist and peak stats, keep thresh
                else:
                    self.plot_current_hist([im_han.hist_and_thresh]) # update hist and get peak stats


                atom_count = np.size(np.where(im_han.atom > 0)[0])  # images with counts above threshold
                empty_count = np.size(np.where(im_han.atom[:im_han.im_num] == 0)[0])
                loading_prob = np.around(atom_count/im_han.im_num, 4)
                # use the binomial distribution to get 1 sigma confidence intervals:
                conf = binom_conf_interval(atom_count, atom_count + empty_count, interval='jeffreys') 
                loading_err = np.around(conf[1] - conf[0], 4)

                peak_stats = [0]*7  # values if calculation fails
                if np.size(im_han.peak_counts) == 2:
                    peak_stats[0] = int(im_han.peak_counts[0]) # background position
                    peak_stats[1] = int(im_han.peak_widths[0]) # background width
                    peak_stats[2] = int(im_han.peak_counts[1]) # signal position
                    peak_stats[3] = int(im_han.peak_widths[1]) # signal width
                    peak_stats[4] = int(im_han.peak_counts[1] - 
                                            im_han.peak_counts[0]) # separation
                    peak_stats[5] = im_han.fidelity            # fidelity
                    peak_stats[6] = im_han.err_fidelity        # error in fidelity
                # set text on labels:
                # counts above:below threshold, images processed, loading probability, error in loading probability,
                # bg count, bg width, signal count, signal width, separation, fidelity, error in fidelity, threshold
                new_stats.append(([str(atom_count) + ' : ' + str(empty_count), str(im_han.im_num),
                    str(loading_prob), str(loading_err)] + list(map(str, peak_stats)) + [str(int(im_han.thresh))]))

        self.update_stat_labels(new_stats)
        return new_stats

    def fit_gaussians(self):
        """Update the histograms and fit two Gaussians, splitting the data at the threshold
        then use the fits to calculate histogram statistics, and set the threshold where the 
        fidelity is maximum"""
        results = []
        for idx, im_han in enumerate(self.image_handler):
            bins, occ, thresh = im_han.histogram()  # get histogram
            bin_mid = (bins[1] - bins[0]) * 0.5 # from edge of bin to middle
            diff = abs(bins - thresh)   # minimum is at the threshold
            thresh_i = np.argmin(diff)  # index of the threshold

            # split the histogram at the threshold value
            best_fits = [fc.fit(bins[:thresh_i]+bin_mid, occ[:thresh_i]),
                            fc.fit(bins[thresh_i:-1]+bin_mid, occ[thresh_i:])]

            for bf in best_fits:
                try:
                    bf.estGaussParam()         # get estimate of parameters
                    # parameters are: amplitude, centre, e^2 width
                    bf.getBestFit(bf.gauss)    # get best fit parameters
                except: return 0               # fit failed, do nothing

            # update image handler's values for peak parameters
            im_han.peak_heights = np.array((best_fits[0].ps[0], best_fits[1].ps[0]))
            im_han.peak_counts = np.array((best_fits[0].ps[1], best_fits[1].ps[1]))
            # e^2 width = 2 * sigma
            im_han.peak_widths = np.array((best_fits[0].ps[2], best_fits[1].ps[2]))/2.

            # update threshold to where fidelity is maximum
            if not self.thresh_toggle.isChecked(): # update thresh if not set by user
                im_han.search_fidelity(best_fits[0].ps[1], best_fits[1].ps[1], n=100)
            else:
                im_han.fidelity, im_han.err_fidelity = np.around(
                                im_han.get_fidelity(), 4) # round to 4 d.p.

            self.plot_current_hist([im_han.histogram]) # clear then update histogram plot
            for bf in best_fits:
                xs = np.linspace(min(bf.x), max(bf.x), 100) # interpolate
                self.hist_canvas[idx].plot(xs, bf.gauss(xs, *bf.ps), pen='b') # plot best fit

            # update atom statistics
            im_han.atom[:im_han.im_num] = im_han.counts[:im_han.im_num] // im_han.thresh   # update atom presence
            atom_count = np.size(np.where(im_han.atom > 0)[0]) # images with counts above threshold
            empty_count = np.size(np.where(im_han.atom[:im_han.im_num] == 0)[0]) # images with counts below threshold
            loading_prob = np.around(atom_count/im_han.im_num, 4) # loading probability
            # use the binomial distribution to get 1 sigma confidence intervals:
            conf = binom_conf_interval(atom_count, atom_count + empty_count, interval='jeffreys') 
            loading_err = np.around(conf[1] - conf[0], 4)

            # counts above : below threshold, images processed, loading probability, error in loading probability, 
            # bg count, bg width, signal count, signal width, separation, fidelity, error in fidelity, threshold
            # note: Gaussian beam waist is 2x standard deviation
            results.append([str(atom_count) + ' : ' + str(empty_count), str(im_han.im_num),
                str(loading_prob), str(loading_err), "%.0f"%best_fits[0].ps[1], "%.0f"%(best_fits[0].ps[2]/2.), 
                "%.0f"%best_fits[1].ps[1], "%.0f"%(best_fits[1].ps[2]/2.), "%.0f"%(best_fits[1].ps[1] - 
                best_fits[0].ps[1]), str(im_han.fidelity), str(im_han.err_fidelity),
                str(int(im_han.thresh))])
        return results

    def update_fit(self, toggle=True):
        """Fit Gaussians to the peaks and use it to get a better estimate of the 
        peak centres and widths. The peaks are not quite Poissonian because of the
        bias from dark counts. Use the fits to get histogram statistics, then set 
        the threshold to maximise fidelity. Iterate until the threshold converges."""
        # only update if a histogram exists
        if self.image_handler[0].im_num > 0 and self.image_handler[1].im_num > 0: 
            # store the previous values
            oldthresh = np.array([x.thresh for x in self.image_handler])
            diff = np.ones(len(oldthresh))        # convergence criterion
            for i in range(20):          # shouldn't need many iterations
                if any(diff < 0.001):
                    break
                new_stats = self.fit_gaussians()
                diff = abs(oldthresh - np.array([x.thresh for x in self.image_handler])) / oldthresh

            if new_stats: # fit_gaussians returns 0 if the fit fails
                self.update_stat_labels(new_stats)
            
    
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

            self.bins_text_edit('reset') # update histogram
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
        if self.dir_watcher:
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
            for im_han in self.image_handler:
                im_han.bin_array = []
            if self.thresh_toggle.isChecked():
                self.plot_current_hist([x.histogram for x in self.image_handler])
            else:
                self.plot_current_hist([x.hist_and_thresh for x in self.image_handler])
            self.bins_text_edit('reset') # to update threshold

        elif self.bin_actions[2].isChecked(): # no update
            try: # disconnect all slots
                self.dir_watcher.event_handler.event_path.disconnect()
            except Exception: pass # if it's already been disconnected 
            
            # just process the image and set the text of the most recent file
            if self.dir_watcher: # check that the dir watcher exists to prevent crash
                for im_han in self.image_handler:
                    self.dir_watcher.event_handler.event_path.connect(im_han.process)
                self.dir_watcher.event_handler.event_path.connect(self.recent_label.setText) # might need a better label
            
    #### #### canvas functions #### #### 
        
    def plot_current_hist(self, hist_functions):
        """Reset the plots to show the current data stored in the image handler.
        hist_functions (must be a list) is used to make the histogram and allows 
        the toggling of different functions that may or may not update the 
        threshold value."""
        # update the histogram and threshold estimate
        for hf in hist_functions:
            bins, occ, thresh = hf()
            idx = np.where([x.histogram == hf or x.hist_and_thresh == hf
                    for x in self.image_handler])[0][0]
            
            self.hist_canvas[idx].clear()
            self.hist_canvas[idx].plot(bins, occ, stepMode=True, pen='k',
                                    fillLevel=0, brush = (220,220,220,220)) # histogram
            self.hist_canvas[idx].plot([thresh]*2, [0, max(occ)], pen='r') # threshold line

    
    def update_im(self, event_path):
        """Receive the event path emitted from the system event handler signal
        display the image from the file in the image canvas"""
        im_vals = self.image_handler[0].load_full_im(event_path)
        self.im_canvas.setImage(im_vals)
        self.im_hist.setLevels(np.min(im_vals), np.max(im_vals))
        
        
    def update_plot(self, event_path):
        """Receive the event path emitted from the system event handler signal
        process the file in the event path with the image handler and update
        the figure"""
        # add the count
        t1 = time.time()
        for im_han in self.image_handler:
            im_han.process(event_path)
        t2 = time.time()
        self.int_time = t2 - t1
        
        # display the name of the most recent file
        self.recent_label.setText('Just processed: '+os.path.basename(event_path))
        
        self.plot_current_hist([x.hist_and_thresh for x in self.image_handler]) # update the displayed plot

        self.plot_time = time.time() - t2

    def update_plot_only(self, event_path):
        """Receive the event path emitted from the system event handler signal
        process the file in the event path with the image handler and update
        the figure but without changing the threshold value"""
        # add the count
        t1 = time.time()
        for im_han in self.image_handler:
            im_han.process(event_path)
        t2 = time.time()
        self.int_time = t2 - t1
        
        # display the name of the most recent file
        self.recent_label.setText('Just processed: '+os.path.basename(event_path))
        
        self.plot_current_hist([x.histogram for x in self.image_handler]) # update the displayed plot
        self.plot_time = time.time() - t2

    def add_stats_to_plot(self, toggle=True):
        """Take the current histogram statistics from the Histogram Statistics labels
        and add the values to the variable plot, saving the parameters to the log
        file at the same time. If any of the labels are empty, do nothing and return -1"""
        stats = self.get_stats() # get statistics from histogram statistics tab labels (list of strings)

        for i in range(len(stats)): # should loop over the atoms
            if not any([s == '' for s in stats[i]]): # only add stats if the fit is successful
                # append current statistics to the histogram handler's list
                if np.size(self.histo_handler[i].vals):
                    self.histo_handler[i].vals = np.append(self.histo_handler[i].vals,
                            [list(map(float, stats[i]))], axis=0)
                else:  # if the array is empty we can't append to it
                    self.histo_handler[i].vals = np.array([list(map(float, stats[i]))]) # needs the right shape
                
                self.update_varplot_axes()  # update the plot with the new values
        
                hist_num = np.size(self.histo_handler[i].vals) // len(self.histo_handler[i].headers) - 1 # index for histograms
                # append histogram stats to log file:
                with open(self.log_file_names[i], 'a') as f:
                    f.write(','.join([str(hist_num)] + stats[i]) + '\n')
                
        return hist_num

    def add_to_varplot(self, hist_han):
        """The user selects which variable they want to display on the plot
        The variables are read from the x and y axis QComboBoxes
        Then the plot is updated with statistics from the histo_handler"""
        if np.size(hist_han.vals) > 0:
            xi = np.where(hist_han.headers == str(self.plot_labels[0].currentText()))[0][0]
            hist_han.xvals = hist_han.vals[:, xi] # set x values
            
            y_label = str(self.plot_labels[1].currentText())
            yi = np.where(hist_han.headers == y_label)[0][0]
            hist_han.yvals = hist_han.vals[:, yi] # set y values

            try:
                self.varplot_canvas.plot(hist_han.xvals, hist_han.yvals, 
                                            pen=None, symbol='o', symbolBrush=self.c[hist_han.i])
                # add error bars if available:
                if 'Loading probability'in y_label or 'Fidelity' in y_label:
                    # estimate sensible beam width at the end of the errorbar
                    if np.size(hist_han.xvals)//2:
                        beam_width = 0.1*(hist_han.xvals[1]-hist_han.xvals[0])
                    else:
                        beam_width = 0.2
                    # add widget for errorbars
                    err_bars = pg.ErrorBarItem(x=hist_han.xvals, 
                                          y=hist_han.yvals, 
                                          height=hist_han.vals[:,yi+1],
                                          beam=beam_width)
                    self.varplot_canvas.addItem(err_bars)
            except Exception: pass # probably wrong length of arrays

    def update_varplot_axes(self, label=''):
        """If the user has set the toggle for the given atom, then plot its 
        histogram statistics on the varplot"""
        self.varplot_canvas.clear()  # remove previous data
        for i in range(len(self.histo_handler)):
            if self.atom_varplot_toggles[self.atomX[i]].isChecked():
                self.add_to_varplot(self.histo_handler[i])
        

    def clear_varplot(self):
        """Clear the plot of histogram statistics by resetting the histo_handler.
        The data is not lost since it has been appended to the log file."""
        for hist_han in self.histo_handler:
            hist_han.__init__ () # empty the stored arrays
        self.varplot_canvas.clear()    # clear the displayed plot


    #### #### save and load data functions #### ####

    def update_stat_labels(self, args):
        """Set the text of the histogram statistics labels to the given arguments.
        args is the results for all the atoms: [[Cs], [Rb]]
        The labels are: 'Counts above : below threshold', 'Number of images processed', 
        'Loading probability', 'Error in loading probability', Background peak count', 'Background peak width', 
        'Signal peak count', 'Signal peak width', 'Separation', 'Fidelity', 'Error in fidelity', 'Threshold'"""
        for i in range(len(self.histo_handler)):
            for ii, label in enumerate(['Counts above : below threshold']+list(self.histo_handler[i].headers[1:])):
                self.stat_labels[self.atomX[i]+label].setText(args[i][ii])


    def get_stats(self):
        """Return the histogram statistics currently stored in the statistics labels
        The data type is a list of strings for all atoms: [[Cs], [Rb]]
        user variable, images processed, loading probability, error in loading probability, 
        bg count, bg width, signal count, signal width, separation, threshold"""
        stats = [[self.var_edit.text()] for i in range(len(self.atomX))]
        for i in range(len(self.histo_handler)): # loop over atoms
            # the first stat label is the ratio bg : signal 
            for label in list(self.histo_handler[i].headers[1:]):
                stats[i].append(self.stat_labels[self.atomX[i]+label].text())

        return stats


    def get_default_path(self, default_path='', option='hist'):
        """If the directory watcher is active, set its results path attribute as the
        default path when a file browser is opened.
        default_path: set the default path if the directory watcher isn't running
        option: 'hist' takes the results path where histograms are stored
                'im' takes the image storage path for loading images
                'log' take the path where log files are saved"""
        date = time.strftime("%Y %B %d", time.localtime()).split(" ") # day, long month, year
        
        if self.dir_watcher and option=='hist':    # make results path the default
            default_path = self.dir_watcher.results_path + r'\%s\%s\%s'%(date[0], date[1], date[2]) 
        elif self.dir_watcher and option=='im':
            default_path = self.dir_watcher.image_storage_path
        elif option=='log':
            default_path = os.path.dirname(self.log_file_names[0])
        
        return default_path


    def load_im_size(self):
        """Get the user to select an image file and then use this to get the image size"""
        default_path = self.get_default_path(option='im')
        try:
            if 'PyQt4' in sys.modules:
                file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')
            elif 'PyQt5' in sys.modules:
                file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')

            for im_han in self.image_handler:
                im_han.set_pic_size(file_name) # sets image handler's pic size
                self.pic_size_edit.setText(str(im_han.pic_size)) # update loaded value
                self.pic_size_label.setText(str(im_han.pic_size)) # update loaded value

        except OSError:
            pass # user cancelled - file not found


    def load_roi(self):
        """Get the user to select an image file and then use this to get the ROI centre"""
        default_path = self.get_default_path(option='im')
        try:
            if 'PyQt4' in sys.modules:
                file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')
            elif 'PyQt5' in sys.modules:
                file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')

            # get pic size from this image in case the user forgot to set it
            for i in range(len(self.atomX)):  # loop over atomic species
                self.image_handler[i].set_pic_size(file_name) # sets image handler's pic size
                self.pic_size_edit.setText(str(self.image_handler[i].pic_size)) # update loaded value
                self.pic_size_label.setText(str(self.image_handler[i].pic_size)) # update loaded value
                # get the position of the max count
                self.image_handler[i].set_roi(im_name=file_name) # sets xc and yc
                self.roi_edits[self.atomX[i]+self.roi_label_text[0]].setText(str(self.image_handler[i].xc)) # update loaded value
                self.roi_edits[self.atomX[i]+self.roi_label_text[1]].setText(str(self.image_handler[i].yc)) 
                self.roi_edits[self.atomX[i]+self.roi_label_text[2]].setText(str(self.image_handler[i].roi_size))
                self.roi_labels[self.atomX[i]+self.roi_label_text[0]].setText(str(self.image_handler[i].xc))
                self.roi_labels[self.atomX[i]+self.roi_label_text[1]].setText(str(self.image_handler[i].yc))
                self.roi_labels[self.atomX[i]+self.roi_label_text[2]].setText(str(self.image_handler[i].roi_size))
                self.rois[i].setPos(self.image_handler[i].xc - self.image_handler[i].roi_size//2, 
                            self.image_handler[i].yc - self.image_handler[i].roi_size//2) # set ROI in image display
                self.rois[i].setSize(self.image_handler[i].roi_size, self.image_handler[i].roi_size)

        except OSError:
            pass # user cancelled - file not found


    def save_hist_data(self, trigger=None, atoms=range(2)):
        """Prompt the user to give a directory to save the histogram data, then save
        atoms specifies which histograms to save, referring to the indices of self.atomX"""
        default_path = self.get_default_path()
        try:
            if 'PyQt4' in sys.modules:
                save_file_name = QFileDialog.getSaveFileName(self, 'Save File', default_path, 'csv(*.csv);;all (*)')
            elif 'PyQt5' in sys.modules:
                save_file_name, _ = QFileDialog.getSaveFileName(self, 'Save File', default_path, 'csv(*.csv);;all (*)')
            
            # don't update the threshold  - trust the user to have already set it
            for i in atoms:
                # save separate histograms for each atom, with the atom name at the start of the file
                self.image_handler[i].save_state(os.path.join(os.path.dirname(save_file_name), 
                        self.atomX[i].replace(' ','')+os.path.basename(save_file_name))) # save histogram
            
            hist_num = self.add_stats_to_plot()
            
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("The following files were saved to directory "+os.path.dirname(save_file_name)+" \n"+
                    "\n".join([self.atomX[i].replace(' ','')+os.path.basename(save_file_name) for i in atoms])+
                    "\n\nand histogram %s was appended to the log files."%hist_num)
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

        except OSError:
            pass # user cancelled - file not found

    def get_choice_idx(self, choice):
        """return the indices where choices match self.atomX
        If the choice is 'None' then idxs remains empty"""
        idxs = [] # indices of the atomic species to include
        if 'All' in choice:
            idxs = range(len(self.atomX))
        else: # find the indexes corresponding to the chosen atoms
            for i in range(len(self.atomX)):
                if self.atomX[i] in choice:
                    idxs.append(i)
        return idxs

    def check_reset(self):
        """Ask the user if they would like to save/reset the current data stored"""
        options = ['None'] + [save+X for save in ['Save first, reset ', 'reset '] 
                                                    for X in ['All ']+self.atomX]

        choice, ok = QInputDialog.getItem(self, 'Confirm Data Replacement',
            "Choose the atomic species that you want to reset:", options, 
            current=0, editable=False)

        idxs = self.get_choice_idx(choice)
        if ok:
            if 'Save' in choice: # prompt user for file name then save
                self.save_hist_data(atoms=idxs)
            for i in idxs:
                self.image_handler[i].reset_arrays() # get rid of old data
                self.hist_canvas[i].clear() # remove old histogram from display

        return choice, ok, idxs


    def load_from_files(self, trigger=None):
        """Prompt the user to select image files to process, then sequentially process
        them and update the histogram"""
        default_path = self.get_default_path(option='im')
        _, ok, _ = self.check_reset() # ask the user to select which atom
            
        if ok: # False if the user cancelled
            try:
                self.recent_label.setText('Processing files...') # comes first otherwise not executed
                if 'PyQt4' in sys.modules:
                    file_list = QFileDialog.getOpenFileNames(self, 'Select Files', default_path, 'Images(*.asc);;all (*)')
                elif 'PyQt5' in sys.modules:
                    file_list, _ = QFileDialog.getOpenFileNames(self, 'Select Files', default_path, 'Images(*.asc);;all (*)')

                for file_name in file_list:
                    for im_han in self.image_handler:
                        im_han.process(file_name)
                    self.recent_label.setText('Just processed: '+os.path.basename(file_name)) # only updates at end of loop
            
                self.update_stats()
                if self.recent_label.text == 'Processing files...':
                    self.recent_label.setText('Finished Processing')

            except OSError:
                pass # user cancelled - file not found

    def load_from_csv(self, trigger=None):
        """Prompt the user to select a csv file to load histogram data from.
        It must have the specific layout that the image_handler saves in."""
        default_path = self.get_default_path()
        _, ok, _ = self.check_reset() # ask the user to select which atom
            
        if ok:
            try:
                # the implementation of QFileDialog changed...
                if 'PyQt4' in sys.modules: 
                    file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'csv(*.csv);;all (*)')
                elif 'PyQt5' in sys.modules:
                    file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'csv(*.csv);;all (*)')
                
                for im_han in self.image_handler:
                    im_han.load_from_csv(file_name)
                self.update_stats()

            except OSError:
                pass # user cancelled - file not found

    def load_image(self, trigger=None):
        """Prompt the user to select an image file to display"""
        default_path = self.get_default_path(option='im')
        
        try:
            if 'PyQt4' in sys.modules:
                file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')
            elif 'PyQt5' in sys.modules:
                file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')
            if file_name:  # avoid crash if the user cancelled
                self.update_im(file_name)

        except OSError:
            pass # user cancelled - file not found

    def load_from_log(self, trigger=None):
        """Prompt the user to select the log file then pass it to the histohandler"""
        default_path = self.get_default_path(option='log')
        
        try:
            # the implementation of QFileDialog changed...
            if 'PyQt4' in sys.modules: 
                file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'dat(*.dat);;all (*)')
            elif 'PyQt5' in sys.modules:
                file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'dat(*.dat);;all (*)')

            for i in range(len(self.atomX)):
                if self.atomX[i] in file_name:
                    self.histo_handler[i].load_from_log(file_name)
            self.update_varplot_axes()

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
        reply = QMessageBox.question(self, 'Confirm Action',
            "Save before closing?", QMessageBox.Save |
            QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Cancel)

        if reply == QMessageBox.Save:
            self.save_hist_data()         # save current state
            
            if self.dir_watcher:          # make sure that the directory watcher stops
                self.dir_watcher.observer.stop()   
            event.accept()
        
        elif reply == QMessageBox.Discard:
            if self.dir_watcher: # make sure that the directory watcher stops
                self.dir_watcher.observer.stop()
            event.accept()
            
        else:
            event.ignore()        

####    ####    ####    #### 

def run():
    """Initiate an app to run the program
    if running in Pylab/IPython then there may already be an app instance"""
    app = QApplication.instance()
    standalone = app is None # false if there is already an app instance
    if standalone: # if there isn't an instance, make one
        app = QApplication(sys.argv) 
        
    main_win = main_window()
    main_win.showMaximized() # display app over the full screen
    if standalone: # if an app instance was made, execute it
        sys.exit(app.exec_()) # when the window is closed, the python code also stops

            
if __name__ == "__main__":
    run()