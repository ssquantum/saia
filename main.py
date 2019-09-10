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
__version__ = '2.2'
import os
import sys
import time
import numpy as np
from astropy.stats import binom_conf_interval
import pyqtgraph as pg    # not as flexible as matplotlib but works a lot better with qt
# some python packages use PyQt4, some use PyQt5...
try:
    from PyQt4.QtCore import QThread, pyqtSignal, QEvent, QRegExp
    from PyQt4.QtGui import (QApplication, QPushButton, QWidget, QLabel, QAction,
            QGridLayout, QMainWindow, QMessageBox, QLineEdit, QIcon, QFileDialog,
            QDoubleValidator, QIntValidator, QComboBox, QMenu, QActionGroup, 
            QTabWidget, QVBoxLayout, QFont, QInputDialog, QRegExpValidator) 
except ModuleNotFoundError:
    from PyQt5.QtCore import QThread, pyqtSignal, QEvent, QRegExp
    from PyQt5.QtGui import (QGridLayout, QMessageBox, QLineEdit, QIcon, 
            QFileDialog, QDoubleValidator, QIntValidator, QComboBox, QMenu, 
            QActionGroup, QVBoxLayout, QFont, QRegExpValidator)
    from PyQt5.QtWidgets import (QApplication, QPushButton, QWidget, QTabWidget,
        QAction, QMainWindow, QLabel, QInputDialog)
# change directory to this file's location
os.chdir(os.path.dirname(os.path.realpath(__file__))) 
import imageHandler as ih # process images to build up a histogram
import histoHandler as hh # collect data from histograms together
import directoryWatcher as dw # use watchdog to get file creation events
import fitCurve as fc   # custom class to get best fit parameters using curve_fit
####    ####    ####    ####

# main GUI window contains all the widgets                
class main_window(QMainWindow):
    """Main GUI window managing an instance of SAIA.

    Use Qt to produce the window where the histogram plot is shown.
    A simple interface allows the user to control the limits of the plot 
    and number of bins in the histogram. Separate tabs are made for 
    settings, multirun options, the histogram, histogram statistics,
    displaying an image, and plotting histogram statistics.
    Separate imageHandler instances are used to analyse two different ROIs.
     - The imageHandler module manages variables associated with individual 
        files and creates the histogram
     - The histoHandler module manages variables associated with the 
        collection of files in several histograms
     - The directoryWatcher module manages file moving, saving, and naming.
     - The fitCurve module stores common functions for curve fitting.
    This GUI was produced with help from http://zetcode.com/gui/pyqt5/.
    Keyword arguments:
    config_file  -- if absolute path to a config file that contains 
        the directories for the directoryWatcher is supplied, then use 
        this instead of the default './config/config.dat'.
    pop_up       -- control whether a pop-up window asks the user to 
        initiate the directoryWatcher. 
        0: don't initiate the directoryWatcher.
        1: tacitly initiate the directoryWatcher.
        2: pop-up window asks the user if they want to initiate.
    """
    def __init__(self, config_file='./config/config.dat', pop_up=2):
        super().__init__()
        self.bias = 697   # bias off set from EMCCD
        self.Nr   = 8.8   # read-out noise from EMCCD
        self.dir_watcher = None  # a button will initiate the dir watcher
        self.atomX = ['Cs ', 'Rb ']   # term labels for atoms
        self.c = [(255,127,14), (31,119,180)] # colours to plot in 
        self.image_handler = [ih.image_handler(i, self.atomX[i]) for i in range(len(self.atomX))] # class to process images
        self.histo_handler = [hh.histo_handler(i, self.atomX[i]) for i in range(len(self.atomX))] # class to process histograms
        self.hist_num = 0 # ID number for the next histogram 
        pg.setConfigOption('background', 'w') # set graph background default white
        pg.setConfigOption('foreground', 'k') # set graph foreground default black
        self.date = time.strftime("%d %b %B %Y", time.localtime()).split(" ") # day short_month long_month year
        self.init_UI(config_file)  # make the widgets
        self.init_DW(pop_up)  # ask the user if they want to start the dir watcher
        self.init_log() # write header to the log file that collects histograms
        self.t0 = time.time()  # time of initiation
        self.int_time = 0      # time taken to process an image
        self.plot_time = 0     # time taken to plot the graph

    def init_log(self):
        """Create a directory for today's date as a subdirectory in the log file path
        then write the header to the log file path defined in config.dat"""
        dir_watcher_dict = dw.dir_watcher.get_dirs(self.config_edit.text()) # static method
        # make subdirectory if it doesn't already exist
        log_file_dir = dir_watcher_dict['Log File Path: ']+'\%s\%s\%s'%(self.date[3],self.date[2],self.date[0]) 
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
                    f.write('#'+', '.join(self.histo_handler[0].stats_dict.keys())+'\n')
       

    def init_UI(self, config_file='./config/config.dat'):
        """Create all of the widget objects required"""
        # grid layout: central main plot, params above, dir watcher status at bottom
        self.centre_widget = QWidget()
        self.tabs = QTabWidget()                   # make tabs for each main display 
        self.centre_widget.layout = QVBoxLayout()
        self.centre_widget.layout.addWidget(self.tabs)
        self.centre_widget.setLayout(self.centre_widget.layout)
        self.setCentralWidget(self.centre_widget)
        
        # validators for user input
        reg_exp = QRegExp(r'([0-9]+(\.[0-9]+)?,?)+')
        comma_validator = QRegExpValidator(reg_exp) # floats and commas
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

        save_hist = QAction('Save histograms', self) # save current hist to csv
        save_hist.triggered.connect(self.save_hist_data)
        hist_menu.addAction(save_hist)

        reset_hist = QAction('Reset histogram', self) # reset hist without loading new data
        reset_hist.triggered.connect(self.check_reset)
        hist_menu.addAction(reset_hist)
        
        load_menu = QMenu('Load histogram data', self)  # drop down menu for loading hist
        load_dir = QAction('From Files', self) # from image files
        load_dir.triggered.connect(self.load_from_files)
        load_menu.addAction(load_dir)
        load_fnums = QAction('From File Numbers', self) # from image file numbers
        load_fnums.triggered.connect(self.load_from_file_nums)
        load_menu.addAction(load_fnums)
        load_csv = QAction('From csv', self) # from csv of hist data
        load_csv.triggered.connect(self.load_from_csv)
        load_menu.addAction(load_csv)
        
        hist_menu.addMenu(load_menu)

        bin_menu = QMenu('Binning', self) # drop down menu for binning options
        bin_options = QActionGroup(bin_menu)  # group together the options
        self.bin_actions = []
        for action_label in ['Automatic', 'Manual', 'No Display', 'No Update']:
            self.bin_actions.append(QAction(action_label, bin_menu, checkable=True, 
                            checked=action_label=='Automatic')) # default is auto
            bin_menu.addAction(self.bin_actions[-1])
            bin_options.addAction(self.bin_actions[-1])
        bin_options.setExclusive(True) # only one option checked at a time
        bin_options.triggered.connect(self.set_bins) # connect the signal
        hist_menu.addMenu(bin_menu)

        
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

        # EMCCD bias offset
        bias_offset_label = QLabel('EMCCD bias offset: ', self)
        settings_grid.addWidget(bias_offset_label, 4,0, 1,1)
        self.bias_offset_edit = QLineEdit(self)
        settings_grid.addWidget(self.bias_offset_edit, 4,1, 1,1)
        self.bias_offset_edit.setText(str(self.bias)) # default
        self.bias_offset_edit.editingFinished.connect(self.CCD_stat_edit)
        self.bias_offset_edit.setValidator(double_validator) # only floats

        # EMCCD readout noise
        read_noise_label = QLabel('EMCCD read-out noise: ', self)
        settings_grid.addWidget(read_noise_label, 5,0, 1,1)
        self.read_noise_edit = QLineEdit(self)
        settings_grid.addWidget(self.read_noise_edit, 5,1, 1,1)
        self.read_noise_edit.setText(str(self.Nr)) # default
        self.read_noise_edit.editingFinished.connect(self.CCD_stat_edit)
        self.read_noise_edit.setValidator(double_validator) # only floats

        
        # change paths used for directory watcher by reloading config file:
        config_label = QLabel('Config File: ')
        settings_grid.addWidget(config_label, 6,0, 1,1)
        self.config_edit = QLineEdit(self)
        settings_grid.addWidget(self.config_edit, 6,1, 1,1)
        self.config_edit.setText(config_file)
        self.config_edit.editingFinished.connect(self.path_text_edit)


        # show paths from the current config file
        self.path_label_text = ['Image Storage Path: ', 'Log File Path: ', 
                'Dexter Sync File: ', 'Image Read Path: ', 'Results Path: ']
        self.path_label = {}
        for i in range(len(self.path_label_text)):
            new_label = QLabel(self.path_label_text[i], self)
            settings_grid.addWidget(new_label, i+7,0, 1,1)
            self.path_label[self.path_label_text[i]] = QLabel('', self)
            settings_grid.addWidget(self.path_label[self.path_label_text[i]], i+7,1, 1,1)
        
        # button to initiate dir watcher
        self.dw_init_button = QPushButton('Initiate directory watcher', self)
        self.dw_init_button.clicked.connect(self.reset_DW) # function to start/stop dir watcher
        self.dw_init_button.resize(self.dw_init_button.sizeHint())
        settings_grid.addWidget(self.dw_init_button, i+8,0, 1,1)

        # toggle to choose whether the dir watcher is active or passive
        self.dw_mode = QPushButton('Active', self, checkable=True, checked=True)
        self.dw_mode.clicked[bool].connect(self.dw_mode_switch)
        settings_grid.addWidget(self.dw_mode, i+8,1, 1,1)

        # label to show status of dir watcher
        self.dw_status_label = QLabel('Stopped', self)  # should errors stop dir watcher???
        settings_grid.addWidget(self.dw_status_label, i+8,2, 1,1)

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        settings_grid.addWidget(self.recent_label, i+9,0, 1,4)
        
        #### tab for multi-run settings ####
        multirun_tab = QWidget()
        multirun_grid = QGridLayout()
        multirun_tab.setLayout(multirun_grid)
        self.tabs.addTab(multirun_tab, "Multirun")

        # dictionary for multirun settings
        self.mr = {'# omit':0, '# hist':100, 'var list':[], 
                'prefix':'0', 'o':0, 'h':0, 'v':0, 
                'measure':0}

        # user chooses an ID as a prefix for the histogram files
        measure_label = QLabel('Measure prefix: ', self)
        multirun_grid.addWidget(measure_label, 0,0, 1,1)
        self.measure_edit = QLineEdit(self)
        multirun_grid.addWidget(self.measure_edit, 0,1, 1,1)
        self.measure_edit.setText(str(self.mr['prefix']))
        
        # user chooses a variable to include in the multi-run
        entry_label = QLabel('User variable: ', self)
        multirun_grid.addWidget(entry_label, 1,0, 1,1)
        self.entry_edit = QLineEdit(self)
        multirun_grid.addWidget(self.entry_edit, 1,1, 1,1)
        self.entry_edit.returnPressed.connect(self.add_var_to_multirun)
        self.entry_edit.setValidator(comma_validator)
        # add the current variable to list
        add_var_button = QPushButton('Add to list', self)
        add_var_button.clicked.connect(self.add_var_to_multirun)
        add_var_button.resize(add_var_button.sizeHint())
        multirun_grid.addWidget(add_var_button, 1,2, 1,1)
        # display current list of user variables
        var_list_label = QLabel('Current list: ', self)
        multirun_grid.addWidget(var_list_label, 2,0, 1,1)
        self.multirun_vars = QLabel('', self)
        multirun_grid.addWidget(self.multirun_vars, 2,1, 1,1)
        # clear the current list of user variables
        clear_vars_button = QPushButton('Clear list', self)
        clear_vars_button.clicked.connect(self.clear_multirun_vars)
        clear_vars_button.resize(clear_vars_button.sizeHint())
        multirun_grid.addWidget(clear_vars_button, 2,2, 1,1)
        
        # choose how many files to omit before starting the next histogram
        omit_label = QLabel('Omit the first N files: ', self)
        multirun_grid.addWidget(omit_label, 3,0, 1,1)
        self.omit_edit = QLineEdit(self)
        multirun_grid.addWidget(self.omit_edit, 3,1, 1,1)
        self.omit_edit.setText(str(self.mr['# omit'])) # default
        self.omit_edit.setValidator(int_validator)

        # choose how many files to have in one histogram
        hist_size_label = QLabel('# files in the histogram: ', self)
        multirun_grid.addWidget(hist_size_label, 4,0, 1,1)
        self.multirun_hist_size = QLineEdit(self)
        multirun_grid.addWidget(self.multirun_hist_size, 4,1, 1,1)
        self.multirun_hist_size.setText(str(self.mr['# hist'])) # default
        self.multirun_hist_size.setValidator(int_validator)

        # choose the directory to save histograms and measure files to
        multirun_dir_button = QPushButton('Choose directory to save to: ', self)
        multirun_grid.addWidget(multirun_dir_button, 5,0, 1,1)
        multirun_dir_button.clicked.connect(self.choose_multirun_dir)
        multirun_dir_button.resize(multirun_dir_button.sizeHint())
        # default directory is the results folder
        self.multirun_save_dir = QLabel(self.get_default_path(option='hist'), self)
        multirun_grid.addWidget(self.multirun_save_dir, 5,1, 1,1)

        # start/abort the multirun
        self.multirun_switch = QPushButton('Start', self, checkable=True)
        self.multirun_switch.clicked[bool].connect(self.multirun_go)
        multirun_grid.addWidget(self.multirun_switch, 6,1, 1,1)
        # pause/restart the multirun
        self.multirun_pause = QPushButton('Resume', self)
        self.multirun_pause.clicked.connect(self.multirun_resume)
        multirun_grid.addWidget(self.multirun_pause, 6,2, 1,1)

        # display current progress
        self.multirun_progress = QLabel(
            'User variable: , omit 0 of 0 files, 0 of 100 histogram files, 0% complete')
        multirun_grid.addWidget(self.multirun_progress, 712,0, 1,3)

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
        
        # toggle whether to fix threshold at user specified value
        self.thresh_toggle = QAction('User Threshold', self, checkable=True)
        self.thresh_toggle.triggered.connect(self.set_thresh)
        hist_menu.addAction(self.thresh_toggle)

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
        user_var_label = QLabel('User Variable: ', self)
        stat_grid.addWidget(user_var_label, 0,0, 1,1)
        self.var_edit = QLineEdit(self)
        self.var_edit.editingFinished.connect(self.set_user_var)
        stat_grid.addWidget(self.var_edit, 0,1, 1,1)
        self.var_edit.setText('0')  # default
        self.var_edit.setValidator(double_validator) # only numbers

        self.stat_labels = {}  # dictionary of stat labels
        # get the list of labels from the histogram handler
        for i, label_text in enumerate(self.histo_handler[0].stats_dict.keys()):
            new_label = QLabel(label_text, self) # description
            stat_grid.addWidget(new_label, i+1,0, 1,1)
            for ii, X in enumerate(self.atomX): 
                self.stat_labels[X+label_text] = QLabel('', self) # value
                stat_grid.addWidget(self.stat_labels[X+label_text], i+1,ii+1, 1,1)
            
        # update statistics
        stat_update = QPushButton('Update statistics', self)
        stat_update.clicked[bool].connect(self.update_stats)
        stat_grid.addWidget(stat_update, i+2,0, 1,1)

        # do Gaussian/Poissonian fit - peaks and widths
        fit_update = QPushButton('Get best fits', self)
        fit_update.clicked[bool].connect(self.update_fit)
        stat_grid.addWidget(fit_update, i+2,1, 1,1)

        # do a Gaussian fit just to the background peak
        self.fit_bg_button = QPushButton('Fit background', self)
        self.fit_bg_button.clicked[bool].connect(self.fit_bg_gaussian)
        stat_grid.addWidget(self.fit_bg_button, i+2,2, 1,1)

        # quickly add the current histogram statistics to the plot
        add_to_plot = QPushButton('Add to plot', self)
        add_to_plot.clicked[bool].connect(self.add_stats_to_plot)
        stat_grid.addWidget(add_to_plot, i+3,1, 1,1)

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
            self.plot_labels[i].addItems(list(self.histo_handler[0].stats_dict.keys())) 
            # connect buttons to update functions
            self.plot_labels[i].activated[str].connect(self.update_varplot_axes)
        # empty text box for the user to write their xlabel
        self.plot_labels.append(QLineEdit(self))
        # position labels in grid
        plot_grid.addWidget(self.plot_labels[0], 7,3, 1,1) # bottom middle
        plot_grid.addWidget(self.plot_labels[1], 2,0, 1,1) # middle left
        plot_grid.addWidget(self.plot_labels[2], 7,4, 1,1) # bottom middle
        # button to clear plot data (it's still saved in the log file)
        clear_varplot = QPushButton('Clear plot', self)
        clear_varplot.clicked[bool].connect(self.clear_varplot)
        plot_grid.addWidget(clear_varplot, 7,0, 1,1)
        # button to save plot data to separate file (it's also in the log file)
        save_varplot = QPushButton('Save plot data', self)
        save_varplot.clicked[bool].connect(self.save_varplot)
        plot_grid.addWidget(save_varplot, 5,0, 1,1)

        #### choose main window position and dimensions: (xpos,ypos,width,height)
        self.setGeometry(100, 100, 850, 700)
        self.setWindowTitle('Single Atom Image Analyser')
        self.setWindowIcon(QIcon('docs/tempicon.png'))
        
    #### #### initiation functions #### #### 

    def init_DW(self, pop_up=2):
        """Ask the user if they want to start the dir watcher or not
        Keyword arguments:
        pop_up       -- control whether a pop-up window asks the user to 
            initiate the directoryWatcher. 
            0: don't initiate the directoryWatcher.
            1: tacitly initiate the directoryWatcher.
            2: pop-up window asks the user if they want to initiate."""
        dir_watcher_dict = dw.dir_watcher.get_dirs(self.config_edit.text()) # static method
        if pop_up == 2: # make pop_up window ask whether you want to initiate
            pad = 0 # make the message box wider by padding out the first line
            for fp in dir_watcher_dict.values():
                if len(fp) > pad:
                    pad = len(fp)
            text = "Loaded from config file."+''.join(['  ']*pad)+".\n"
            text += dw.dir_watcher.print_dirs(dir_watcher_dict.items()) # static method
            text += "\nStart the directory watcher with these settings?"
            reply = QMessageBox.question(self, 'Initiate the Directory Watcher',
                text, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                self.reset_DW() # takes the config file from config_edit
        elif pop_up == 1:
            self.reset_DW()
        elif pop_up == 0:
            pass

    def remove_im_files(self):
        """Ask the user if they want to remove image files from the read image
        path since the dir watcher only notices file created not modified"""
        text = 'The directory watcher only notices file creation events, not modifications.\n'
        text += 'Therefore the image read path must be emptied so new files can be created.\n'
        text += '\nDelete the following files from '+self.dir_watcher.image_read_path+"?\n"
        file_list = []
        for file_name in os.listdir(self.dir_watcher.image_read_path):
            if '.asc' in file_name:
                file_list.append(file_name)
                if len(file_list) < 10:
                    text += "\t - " + file_name + "\n"
        text += '(Total %s files found.)\n'%len(file_list)

        reply = QMessageBox.question(self, 'Remove Initial Image files?',
            text, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            for file_name in file_list:
                os.remove(os.path.join(self.dir_watcher.image_read_path, file_name))

    def dw_mode_switch(self):
        """Change the dw_mode switch so that when in active mode it reads active,
        when in passive mode it reads passive"""
        if self.dw_mode.isChecked():
            self.dw_mode.setText('Active')
        else:
            self.dw_mode.setText('Passive')
            
    def reset_DW(self):
        """Initiate the dir watcher. If there is already one running, stop the 
        thread and delete the instance to ensure it doesn't run in the 
        background (which might overwrite files)."""
        if self.dir_watcher: # check if there is a current thread
            self.print_times("ms")  # prints performance of dir_watcher
            self.dir_watcher.observer.stop() # ensure that the old thread stops
            self.dir_watcher = None
            self.dw_status_label.setText("Stopped")
            self.dw_init_button.setText('Initiate directory watcher') # turns on
            self.recent_label.setText('')

        else: 
            self.dir_watcher = dw.dir_watcher(
                    config_file=self.config_edit.text(),
                    active=self.dw_mode.isChecked()) # instantiate dir watcher
            self.remove_im_files() # prompt to remove image files
            self.dir_watcher.event_handler.event_path.connect(self.update_plot) # default
            self.dir_watcher.event_handler.sync_dexter() # get the current Dexter file number
            self.dw_status_label.setText("Running")
            # get current date
            self.date = self.dir_watcher.date
            date_str = ' '.join([self.date[0]]+self.date[2:])
            self.init_log() # make a new log file
            pad = 0 # make the message box wider by padding out the first line
            for fp in self.dir_watcher.dirs_dict.values():
                if len(fp) > pad:
                    pad = len(fp)
            msg = QMessageBox() # pop up box to confirm it's started
            msg.setIcon(QMessageBox.Information)
            msg.setText(
                "Directory Watcher initiated in " + self.dw_mode.text()
                + " mode with settings:" + ''.join([' ']*pad) + ".\n\n" + 
                "date\t\t\t--" + date_str + "\n\n" +
                self.dir_watcher.print_dirs(self.dir_watcher.dirs_dict.items()))
            msg.setStandardButtons(QMessageBox.Ok)
            msg.setFixedSize(msg.sizeHint())
            msg.exec_()
            self.dw_init_button.setText('Stop directory watcher') # turns off
            # display current date on window title
            self.setWindowTitle('Single Atom Image Analyser --- ' + date_str)

            # set current file paths
            for key, value in self.dir_watcher.dirs_dict.items():
                self.path_label[key].setText(value)

    #### #### user input functions #### #### 

    def set_user_var(self, text=''):
        """When the user finishes editing the var_edit line edit, update the displayed 
        user variable and assign it in the temp_vals of the histo_handler"""
        for i in range(len(self.histo_handler)):
            self.histo_handler[i].temp_vals['User variable'] = self.var_edit.text()
            self.stat_labels[self.atomX[i]+'User variable'].setText(self.var_edit.text())

    def path_text_edit(self, text=''):
        """The user finishes editing an edit text box by pressing return or clicking
        somewhere else, then the text is sent to this function to reload the config file. 
        The dir watcher is not updated unless the 'initiate dir watcher' button is used."""
        dw_dict = dw.dir_watcher.get_dirs(self.config_edit.text())
        for key, value in dw_dict.items():
            self.path_label[key].setText(value)
            
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

    def CCD_stat_edit(self):
        """Update the values used for the EMCCD bias offset and readout noise"""
        if self.bias_offset_edit.text(): # check the label isn't empty
            self.bias = float(self.bias_offset_edit.text())
        if self.read_noise_edit.text():
            self.Nr = float(self.read_noise_edit.text())

    def add_var_to_multirun(self):
        """When the user hits enter or the 'Add to list' button, add the 
        text from the entry edit to the list of user variables that will 
        be used for the multi-run. For speed, you can enter a range in 
        the form start,stop,step,repeat. If the multi-run has already
        started, do nothing."""
        if not self.multirun_switch.isChecked():
            new_var = list(map(float, [v for v in self.entry_edit.text().split(',') if v]))
            if np.size(new_var) == 1: # just entered a single variable
                self.mr['var list'].append(new_var[0])
                # empty the text edit so that it's quicker to enter a new variable
                self.entry_edit.setText('') 

            elif np.size(new_var) == 3: # range, with no repeats
                self.mr['var list'] += list(np.arange(new_var[0], new_var[1], new_var[2]))
            elif np.size(new_var) == 4: # range, with repeats
                self.mr['var list'] += list(np.arange(new_var[0], new_var[1],
                                            new_var[2]))*int(new_var[3])
            # display the whole list
            self.multirun_vars.setText(','.join(list(map(str, self.mr['var list']))))

    def clear_multirun_vars(self):
        """Reset the list of user variables to be used in the multi-run.
        If the multi-run is already running, don't do anything"""
        if not self.multirun_switch.isChecked():
            self.mr['var list'] = []
            self.multirun_vars.setText('')

    def choose_multirun_dir(self):
        """Allow the user to choose the directory where the histogram .csv
        files and the measure .dat file will be saved as part of the multi-run"""
        default_path = self.get_default_path(option='hist')
        try:
            dir_path = QFileDialog.getExistingDirectory(self, "Select Directory", default_path)
            self.multirun_save_dir.setText(dir_path)
        except OSError:
            pass # user cancelled - file not found
        
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

    def get_correlation(self, atom1, atom2, out_type='str'):
        """Given atom arrays 1 and 2 containing the comparison of the counts to the threshold value,
        output the cases where the image contains: no atoms, only atom1, only atom2, both atoms
        out_type: str - give the number in each case as a string
                index - give the indexes satisfying each case.
        NB: this doesn't crash when atom1, atom2 have different length, but it is assuming that the
        indexing refers to the same files"""
        no_atom = np.where(np.isin(np.where(atom1 == 0)[0], np.where(atom2 == 0)[0]))[0] # neither atom present
        only_1  = np.where(np.isin(np.where(atom1 >  0)[0], np.where(atom2 == 0)[0]))[0] # only atom1 present
        only_2  = np.where(np.isin(np.where(atom1 == 0)[0], np.where(atom2 >  0)[0]))[0] # only atom2 present
        both    = np.where(np.isin(np.where(atom1 >  0)[0], np.where(atom2 >  0)[0]))[0] # both atoms present
        if out_type == 'str':
            return list(map(str, map(np.size, [no_atom, only_1, only_2, both])))
        elif out_type == 'index':
            return [no_atom, only_1, only_2, both]
        
    def dappend(self, key, value):
        """Shorthand for appending a new value to the array in the histo_handler
        dictionary with the given key. Also update the temp values so that they
        are always the most recent value. Then update the label in the statistics
        tab to display the new value."""
        for idx, hh in enumerate(self.histo_handler):
            hh.stats_dict[key] = np.append(
                    hh.stats_dict[key], 
                    np.array(value).astype(hh.stats_dict[key].dtype))
            hh.temp_vals[key] = value
            self.stat_labels[self.atomX[idx]+key].setText(str(value))

    def update_stats(self, toggle=True):
        """Update the statistics from the current histogram in order to save them
        image_handler uses a peak finding algorithm to get the peak positions and widths
        from these a threshold can be set. If the user thresh_toggle is checked then the
        threshold will not be updated.
        The histo_handler stores temporary values that we might not yet want to add to
        the plot."""
        atom_list = [] 
        for i in range(len(self.image_handler)):
            if self.image_handler[i].im_num > 0: # only update if a histogram exists
                if self.thresh_toggle.isChecked(): # using manual threshold
                    self.plot_current_hist([self.image_handler[i].histogram]) # update hist and peak stats, keep thresh
                else:
                    self.plot_current_hist([self.image_handler[i].hist_and_thresh]) # update hist and get peak stats

                atom_array = self.image_handler[i].atom[:self.image_handler[i].im_num]    # images with counts above threshold
                atom_list.append(atom_array)
                above_idxs = np.where(atom_array > 0)[0] # index of images with counts above threshold
                atom_count = np.size(above_idxs)  # number of images with counts above threshold
                above = self.image_handler[i].counts[above_idxs] # counts above threshold
                below_idxs = np.where(atom_array == 0)[0] # index of images with counts below threshold
                empty_count = np.size(below_idxs) # number of images with counts below threshold
                below = self.image_handler[i].counts[below_idxs] # counts below threshold
                # use the binomial distribution to get 1 sigma confidence intervals:
                conf = binom_conf_interval(atom_count, atom_count + empty_count, interval='jeffreys')
                loading_prob = atom_count/self.image_handler[i].im_num # fraction of images above threshold
                uplperr = conf[1] - loading_prob # 1 sigma confidence above mean
                lolperr = loading_prob - conf[0] # 1 sigma confidence below mean
                # store the calculated histogram statistics as temp, don't add to plot
                self.histo_handler[i].temp_vals['Hist ID'] = int(self.hist_num)
                file_list = [x for x in self.image_handler[i].files if x]
                self.histo_handler[i].temp_vals['Start file #'] = min(map(int, file_list))
                self.histo_handler[i].temp_vals['End file #'] = max(map(int, file_list))
                self.histo_handler[i].temp_vals['ROI xc ; yc ; size'] = ' ; '.join([self.roi_edits[self.atomX[i]+label].text()
                        for label in self.roi_label_text])
                self.histo_handler[i].temp_vals['User variable'] = float(self.var_edit.text())
                self.histo_handler[i].temp_vals['Number of images processed'] = self.image_handler[i].im_num
                self.histo_handler[i].temp_vals['Counts above : below threshold'] = str(atom_count) + ' : ' + str(empty_count)
                self.histo_handler[i].temp_vals['Loading probability'] = np.around(loading_prob, 4)
                self.histo_handler[i].temp_vals['Error in Loading probability'] = np.around((uplperr+lolperr)*0.5, 4)
                self.histo_handler[i].temp_vals['Lower Error in Loading probability'] = np.around(lolperr, 4)
                self.histo_handler[i].temp_vals['Upper Error in Loading probability'] = np.around(uplperr, 4)
                if np.size(self.image_handler[i].peak_counts) == 2:
                    self.histo_handler[i].temp_vals['Background peak count'] = int(self.image_handler[i].peak_counts[0])
                    # assume bias offset is self.bias, readout noise standard deviation Nr
                    if self.Nr**2+self.image_handler[i].peak_counts[0]-self.bias > 0:
                        self.histo_handler[i].temp_vals['sqrt(Nr^2 + Nbg)'] = int((self.Nr**2+self.image_handler[i].peak_counts[0]-self.bias)**0.5)
                    else: # don't take the sqrt of a -ve number
                        self.histo_handler[i].temp_vals['sqrt(Nr^2 + Nbg)'] = 0
                    bgw = self.image_handler[i].peak_widths[0] # fitted background peak width
                    self.histo_handler[i].temp_vals['Background peak width'] = int(bgw)
                    self.histo_handler[i].temp_vals['Error in Background peak count'] = np.around(self.image_handler[i].peak_widths[0] / empty_count**0.5, 2)
                    self.histo_handler[i].temp_vals['Background mean'] = np.around(np.mean(below), 1)
                    self.histo_handler[i].temp_vals['Background standard deviation'] = np.around(np.std(below, ddof=1), 1)
                    self.histo_handler[i].temp_vals['Signal peak count'] = int(self.image_handler[i].peak_counts[1])
                    # assume bias offset is self.bias, readout noise standard deviation Nr
                    if self.Nr**2+self.image_handler[i].peak_counts[1]-self.bias > 0:
                        self.histo_handler[i].temp_vals['sqrt(Nr^2 + Ns)'] = int((self.Nr**2+self.image_handler[i].peak_counts[1]-self.bias)**0.5)
                    else: # don't take the sqrt of a -ve number
                        self.histo_handler[i].temp_vals['sqrt(Nr^2 + Ns)'] = 0
                    siw = self.image_handler[i].peak_widths[1] # fitted signal peak width
                    self.histo_handler[i].temp_vals['Signal peak width'] = int(siw)
                    self.histo_handler[i].temp_vals['Error in Signal peak count'] = np.around(self.image_handler[i].peak_widths[1] / atom_count**0.5, 2)
                    self.histo_handler[i].temp_vals['Signal mean'] = np.around(np.mean(above), 1)
                    self.histo_handler[i].temp_vals['Signal standard deviation'] = np.around(np.std(above, ddof=1), 1)
                    sep = self.image_handler[i].peak_counts[1] - self.image_handler[i].peak_counts[0] # separation of fitted peaks
                    self.histo_handler[i].temp_vals['Separation'] = int(sep)
                    seperr = np.sqrt(self.image_handler[i].peak_widths[0]**2 / empty_count
                                    + self.image_handler[i].peak_widths[1]**2 / atom_count) # propagated error in separation
                    self.histo_handler[i].temp_vals['Error in Separation'] = np.around(seperr, 2)
                    self.histo_handler[i].temp_vals['Fidelity'] = self.image_handler[i].fidelity
                    self.histo_handler[i].temp_vals['Error in Fidelity'] = self.image_handler[i].err_fidelity
                    self.histo_handler[i].temp_vals['S/N'] = np.around(sep / np.sqrt(bgw**2 + siw**2), 2)
                    # fractional error in the error is 1/sqrt(2N - 2)
                    self.histo_handler[i].temp_vals['Error in S/N'] = np.around(
                        self.histo_handler[i].temp_vals['S/N'] * np.sqrt((seperr/sep)**2 + 
                        (bgw**2/(2*empty_count - 2) + siw**2/(2*atom_count - 2))/(bgw**2 + siw**2)), 2)
                else:
                    for key in ['Background peak count', 'sqrt(Nr^2 + Nbg)', 'Background peak width', 
                    'Error in Background peak count', 'Signal peak count', 'sqrt(Nr^2 + Ns)', 
                    'Signal peak width', 'Error in Signal peak count', 'Separation', 'Error in Separation', 
                    'Fidelity', 'Error in Fidelity', 'S/N', 'Error in S/N']:
                        self.histo_handler[i].temp_vals[key] = 0
                self.histo_handler[i].temp_vals['Threshold'] = int(self.image_handler[i].thresh)
            
        # calculate correlations:  - assuming only two atoms!
        if self.histo_handler[i].temp_vals['Background peak count']:
            no_atom, only_1, only_2, both = self.get_correlation(atom_list[0], atom_list[1], out_type='str')
            self.histo_handler[0].temp_vals['No atom'] = no_atom # neither atom present
            self.histo_handler[1].temp_vals['No atom'] = no_atom # neither atom present
            self.histo_handler[0].temp_vals['Single atom'] = only_1 # just one atom present
            self.histo_handler[1].temp_vals['Single atom'] = only_2 # just one atom present
            self.histo_handler[0].temp_vals['Both atoms'] = both # both atoms present
            self.histo_handler[1].temp_vals['Both atoms'] = both # both atoms present
            
        # display the new statistics in the labels
        for idx, hh in enumerate(self.histo_handler):
            for key, val in hh.temp_vals.items():
                self.stat_labels[self.atomX[idx]+key].setText(str(val))


    def fit_gaussians(self, store_stats=False):
        """Update the histogram and fit two Gaussians, splitting the data at the threshold
        then use the fits to calculate histogram statistics, and set the threshold where the 
        fidelity is maximum. If the store_stats Boolean is True, append the calculated values
        to the histo_handler's statistics dictionary."""
        atom_list = []
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
                    # parameters are: amplitude, centre, standard deviation
                    bf.getBestFit(bf.gauss)    # get best fit parameters
                except: return 0               # fit failed, do nothing
            # update image handler's values for peak parameters
            im_han.peak_heights = np.array((best_fits[0].ps[0], best_fits[1].ps[0]))
            im_han.peak_counts = np.array((best_fits[0].ps[1], best_fits[1].ps[1]))
            im_han.peak_widths = np.array((best_fits[0].ps[2], best_fits[1].ps[2]))

            # update threshold to where fidelity is maximum
            if not self.thresh_toggle.isChecked(): # update thresh if not set by user
                im_han.search_fidelity(best_fits[0].ps[1], best_fits[0].ps[2], 
                                                            best_fits[1].ps[1], n=100)
            else:
                im_han.fidelity, im_han.err_fidelity = np.around(
                                im_han.get_fidelity(), 4) # round to 4 d.p.

            self.plot_current_hist([im_han.histogram]) # clear then update histogram plot
            for bf in best_fits:
                xs = np.linspace(min(bf.x), max(bf.x), 100) # interpolate
                self.hist_canvas[idx].plot(xs, bf.gauss(xs, *bf.ps), pen='b') # plot best fit

            # update atom statistics
            im_han.atom[:im_han.im_num] = im_han.counts[:im_han.im_num] // im_han.thresh   # update atom presence
            atom_array = im_han.atom[:im_han.im_num]
            atom_list.append(atom_array)
            above_idxs = np.where(atom_array > 0)[0] # index of images with counts above threshold
            atom_count = np.size(above_idxs)  # number of images with counts above threshold
            above = im_han.counts[above_idxs] # counts above threshold
            below_idxs = np.where(atom_array == 0)[0] # index of images with counts below threshold
            empty_count = np.size(below_idxs) # number of images with counts below threshold
            below = im_han.counts[below_idxs] # counts below threshold
            loading_prob = atom_count/im_han.im_num # loading probability
            # use the binomial distribution to get 1 sigma confidence intervals:
            conf = binom_conf_interval(atom_count, atom_count + empty_count, interval='jeffreys') 
            uplperr = conf[1] - loading_prob # 1 sigma confidence above mean
            lolperr = loading_prob - conf[0] # 1 sigma confidence below mean
            # store the calculated histogram statistics as temp, don't add to plot
        if store_stats:
            for idx, im_han in enumerate(self.image_handler):
                self.histo_handler[idx].temp_vals['Hist ID'] = int(self.hist_num)
                self.histo_handler[idx].temp_vals['User variable'] = float(self.var_edit.text())
                file_list = [x for x in self.image_handler[idx].files if x]
                self.histo_handler[idx].temp_vals['Start file #'] = min(map(int, file_list))
                self.histo_handler[idx].temp_vals['End file #'] = max(map(int, file_list))
                self.histo_handler[idx].temp_vals['ROI xc ; yc ; size'] = ' ; '.join([self.roi_edits[self.atomX[i]+label].text()
                                        for label in self.roi_label_text])
                self.histo_handler[idx].temp_vals['Number of images processed'] = self.image_handler[idx].im_num
                self.histo_handler[idx].temp_vals['Counts above : below threshold'] = str(atom_count) + ' : ' + str(empty_count)
                self.histo_handler[idx].temp_vals['Loading probability'] = np.around(loading_prob, 4)
                self.histo_handler[idx].temp_vals['Error in Loading probability'] = np.around((uplperr + lolperr)*0.5, 4)
                self.histo_handler[idx].temp_vals['Lower Error in Loading probability'] = np.around(lolperr, 4)
                self.histo_handler[idx].temp_vals['Upper Error in Loading probability'] = np.around(uplperr, 4)
                self.histo_handler[idx].temp_vals['Background peak count'] = int(best_fits[0].ps[1])
                # assume bias offset is self.bias, readout noise standard deviation Nr
                if self.Nr**2+best_fits[0].ps[1]-self.bias > 0:
                    self.histo_handler[idx].temp_vals['sqrt(Nr^2 + Nbg)'] = int((self.Nr**2+best_fits[0].ps[1]-self.bias)**0.5)
                else: # don't take the sqrt of a -ve number
                    self.histo_handler[idx].temp_vals['sqrt(Nr^2 + Nbg)'] = 0
                bgw = best_fits[0].ps[2] # fitted background peak width
                self.histo_handler[idx].temp_vals['Background peak width'] = int(bgw)
                self.histo_handler[idx].temp_vals['Error in Background peak count'] = np.around(best_fits[0].ps[2] / empty_count**0.5, 2)
                self.histo_handler[idx].temp_vals['Background mean'] = np.around(np.mean(below), 1)
                self.histo_handler[idx].temp_vals['Background standard deviation'] = np.around(np.std(below, ddof=1), 1)
                self.histo_handler[idx].temp_vals['Signal peak count'] = int(best_fits[1].ps[1])
                # assume bias offset is self.bias, readout noise standard deviation Nr
                if self.Nr**2+best_fits[1].ps[1]-self.bias > 0:
                    self.histo_handler[idx].temp_vals['sqrt(Nr^2 + Ns)'] = int((self.Nr**2+best_fits[1].ps[1]-self.bias)**0.5)
                else:
                    self.histo_handler[idx].temp_vals['sqrt(Nr^2 + Ns)'] = 0
                siw = best_fits[1].ps[2] # fitted signal peak width
                self.histo_handler[idx].temp_vals['Signal peak width'] = int(siw)
                self.histo_handler[idx].temp_vals['Error in Signal peak count'] = np.around(best_fits[1].ps[2] / atom_count**0.5, 2)
                self.histo_handler[idx].temp_vals['Signal mean'] = np.around(np.mean(above), 1)
                self.histo_handler[idx].temp_vals['Signal standard deviation'] = np.around(np.std(above, ddof=1), 1)
                sep = best_fits[1].ps[1] - best_fits[0].ps[1] # separation of fitted peak centres
                self.histo_handler[idx].temp_vals['Separation'] = int(sep)
                seperr = np.sqrt(best_fits[0].ps[2]**2 / empty_count + best_fits[1].ps[2]**2 / atom_count) # error in separation
                self.histo_handler[idx].temp_vals['Error in Separation'] = np.around(seperr, 2)
                self.histo_handler[idx].temp_vals['Fidelity'] = self.image_handler[idx].fidelity
                self.histo_handler[idx].temp_vals['Error in Fidelity'] = self.image_handler[idx].err_fidelity
                self.histo_handler[idx].temp_vals['S/N'] = np.around(sep / np.sqrt(bgw**2 + siw**2), 2)
                # fractional error in the error is 1/sqrt(2N - 2)
                self.histo_handler[idx].temp_vals['Error in S/N'] = np.around(
                            self.histo_handler[idx].temp_vals['S/N'] * np.sqrt((seperr/sep)**2 + 
                            (bgw**2/(2*empty_count - 2) + siw**2/(2*atom_count - 2))/(bgw**2 + siw**2)), 2)
                self.histo_handler[idx].temp_vals['Threshold'] = int(self.image_handler[idx].thresh)
                if self.histo_handler[idx].temp_vals['Background peak count']:
                    no_atom, only_1, only_2, both = self.get_correlation(atom_list[0], atom_list[1], out_type='str')
                    self.histo_handler[0].temp_vals['No atom'] = no_atom # neither atom present
                    self.histo_handler[1].temp_vals['No atom'] = no_atom # neither atom present
                    self.histo_handler[0].temp_vals['Single atom'] = only_1 # just one atom present
                    self.histo_handler[1].temp_vals['Single atom'] = only_2 # just one atom present
                    self.histo_handler[0].temp_vals['Both atoms'] = both # both atoms present
                    self.histo_handler[1].temp_vals['Both atoms'] = both # both atoms present
                # display the new statistics in the labels
                for key, val in self.histo_handler[idx].temp_vals.items():
                    self.stat_labels[self.atomX[idx]+key].setText(str(val))
        return 1 # fit successful
        
    def fit_bg_gaussian(self, store_stats=False):
        """Assume that there is only one peak in the histogram as there is no single
        atom signal. Fit a Gaussian to this peak."""
        for i in range(len(self.image_handler)):
            n = self.image_handler[i].im_num # number of images processed
            c = self.image_handler[i].counts[:n] # integrated counts
            bins, occ, _ = self.image_handler[i].histogram()  # get histogram
            bin_mid = (bins[1] - bins[0]) * 0.5 # from edge of bin to middle
            # make a Gaussian fit to the peak
            best_fit = fc.fit(bins[:-1]+bin_mid, occ)
            try:
                best_fit.estGaussParam()
                # parameters are: amplitude, centre, standard deviation
                best_fit.getBestFit(best_fit.gauss)    # get best fit parameters
            except: return 0               # fit failed, do nothing
            # calculate mean and std dev from the data
            mu, sig = np.mean(c), np.std(c, ddof=1)
            # best_fit.ps = [best_fit.ps[0], mu, sig] # use the peak from the fit
            lperr = np.around(binom_conf_interval(0, n, interval='jeffreys')[1], 4) # upper 1 sigma confidence
            # update image handler's values for peak parameters
            self.image_handler[i].peak_heights = np.array((best_fit.ps[0], 0))
            self.image_handler[i].peak_counts = np.array((best_fit.ps[1], 0))
            self.image_handler[i].peak_widths = np.array((best_fit.ps[2], 0))
            self.plot_current_hist([self.image_handler[i].histogram]) # clear then update histogram plot
            xs = np.linspace(min(best_fit.x), max(best_fit.x), 100) # interpolate
            self.hist_canvas.plot(xs, best_fit.gauss(xs, *best_fit.ps), pen='b') # plot best fit
            # store the calculated histogram statistics as temp, don't add to plot
            self.histo_handler[i].temp_vals['Hist ID'] = int(self.hist_num)
            file_list = [x for x in self.image_handler[i].files if x]
            self.histo_handler[i].temp_vals['Start file #'] = min(map(int, file_list))
            self.histo_handler[i].temp_vals['End file #'] = max(map(int, file_list))
            self.histo_handler[i].temp_vals['ROI xc ; yc ; size'] = ' ; '.join([self.roi_edits[self.atomX[i]+label].text()
                        for label in self.roi_label_text])
            self.histo_handler[i].temp_vals['User variable'] = float(self.var_edit.text())
            self.histo_handler[i].temp_vals['Number of images processed'] = n
            self.histo_handler[i].temp_vals['Counts above : below threshold'] = '0 : ' + str(n)
            self.histo_handler[i].temp_vals['Loading probability'] = 0
            self.histo_handler[i].temp_vals['Error in Loading probability'] = lperr
            self.histo_handler[i].temp_vals['Lower Error in Loading probability'] = 0
            self.histo_handler[i].temp_vals['Upper Error in Loading probability'] = lperr
            self.histo_handler[i].temp_vals['Background peak count'] = int(best_fit.ps[1])
            # assume bias offset is self.bias, readout noise standard deviation Nr
            if self.Nr**2+mu-self.bias:
                self.histo_handler[i].temp_vals['sqrt(Nr^2 + Nbg)'] = int((self.Nr**2+mu-self.bias)**0.5)
            else: # don't take the sqrt of a -ve number
                self.histo_handler[i].temp_vals['sqrt(Nr^2 + Nbg)'] = 0
            self.histo_handler[i].temp_vals['Background peak width'] = int(best_fit.ps[2])
            self.histo_handler[i].temp_vals['Error in Background peak count'] = np.around(best_fit.ps[2] / n**0.5, 4)
            self.histo_handler[i].temp_vals['Background mean'] = np.around(mu, 1)
            self.histo_handler[i].temp_vals['Background standard deviation'] = np.around(sig, 1)
            self.histo_handler[i].temp_vals['Signal peak count'] = 0
            self.histo_handler[i].temp_vals['sqrt(Nr^2 + Ns)'] = 0
            self.histo_handler[i].temp_vals['Signal peak width'] = 0
            self.histo_handler[i].temp_vals['Error in Signal peak count'] = 0
            self.histo_handler[i].temp_vals['Signal mean'] = 0
            self.histo_handler[i].temp_vals['Signal standard deviation'] = 0
            self.histo_handler[i].temp_vals['Separation'] = 0
            self.histo_handler[i].temp_vals['Error in Separation'] = 0
            self.histo_handler[i].temp_vals['Fidelity'] = 0
            self.histo_handler[i].temp_vals['Error in Fidelity'] = 0
            self.histo_handler[i].temp_vals['S/N'] = 0
            self.histo_handler[i].temp_vals['Error in S/N'] = 0
            self.histo_handler[i].temp_vals['Threshold'] = int(self.image_handler[i].thresh)
            # display the new statistics in the labels
            for key, val in self.histo_handler[i].temp_vals.items():
                self.stat_labels[self.atomX[i]+key].setText(str(val))

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
                success = self.fit_gaussians()
                diff = abs(oldthresh - np.array([x.thresh for x in self.image_handler])) / oldthresh

            if success: # fit_gaussians returns 0 if the fit fails
                self.fit_gaussians(store_stats=True) # add new stats to histo_handler
        return success
            
    
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

    def multirun_go(self, toggle):
        """Initiate the multi-run: omit N files, save a histogram of M files, and
        repeat for the user variables in the list. If the button is pressed during
        the multi-run, save the current histogram, save the measure file, then
        return to normal operation of the dir_watcher"""
        if toggle and np.size(self.mr['var list']) > 0:
            self.check_reset()
            self.plot_current_hist([x.histogram for x in self.image_handler])
            try: # disconnect all slots
                self.dir_watcher.event_handler.event_path.disconnect() 
            except Exception: pass # already disconnected
            if self.dir_watcher:
                if self.multirun_save_dir.text() == '':
                    self.choose_multirun_dir()
                self.dir_watcher.event_handler.event_path.connect(self.multirun_step)
                self.mr['# omit'] = int(self.omit_edit.text()) # number of files to omit
                self.mr['# hist'] = int(self.multirun_hist_size.text()) # number of files in histogram                
                self.mr['o'], self.mr['h'], self.mr['v'] = 0, 0, 0 # counters for different stages of multirun
                self.mr['prefix'] = self.measure_edit.text() # prefix for histogram files 
                self.multirun_switch.setText('Abort')
                self.clear_varplot() # varplot cleared so it only has multirun data
                self.multirun_progress.setText(       # update progress label
                    'User variable: %s, omit %s of %s files, %s of %s histogram files, 0%% complete'%(
                        self.mr['var list'][self.mr['v']], self.mr['o'], self.mr['# omit'],
                        self.mr['h'], self.mr['# hist']))
            else: # If dir_watcher isn't running, can't start multirun.
                self.multirun_switch.setChecked(False)
        else: # cancel the multi-run
            self.set_bins() # reconnect the dir_watcher
            self.multirun_switch.setText('Start') # reset button text
            self.multirun_progress.setText(       # update progress label
                'Stopped at - User variable: %s, omit %s of %s files, %s of %s histogram files, %.3g%% complete'%(
                    self.mr['var list'][self.mr['v']], self.mr['o'], self.mr['# omit'],
                    self.mr['h'], self.mr['# hist'], 100 * ((self.mr['# omit'] + self.mr['# hist']) * 
                    self.mr['v'] + self.mr['o'] + self.mr['h']) / (self.mr['# omit'] + self.mr['# hist']) / 
                    np.size(self.mr['var list'])))

    def multirun_resume(self):
        """If the button is clicked, resume the multi-run where it was left off.
        If the multirun is already running, do nothing."""
        if not self.multirun_switch.isChecked(): 
            self.multirun_switch.setChecked(True)
            self.multirun_switch.setText('Abort')
            try: # disconnect all slots
                self.dir_watcher.event_handler.event_path.disconnect() 
            except Exception: pass # already disconnected
            if self.dir_watcher:
                self.dir_watcher.event_handler.event_path.connect(self.multirun_step)
    
    def set_bins(self, action=None):
        """Check which of the bin action menu bar options is checked.
        If the toggle is Automatic, use automatic histogram binning.
        If the toggle is Manual, read in values from the line edit 
        widgets.
        If the toggle is No Display, disconnect the dir watcher new event signal
        from the plot update. Still processes files but doesn't show on histogram
        If the toggle is No Update, disconnect the dir watcher new event signal
        from the image handler entirely. Files are copied but not processed for
        the histogram."""
        if not self.multirun_switch.isChecked(): # don't interrupt multirun
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
            elif self.bin_actions[2].isChecked() or self.bin_actions[3].isChecked(): # No Display or No Update
                try: # disconnect all slots
                    self.dir_watcher.event_handler.event_path.disconnect()
                except Exception: pass # if it's already been disconnected 
                # just process the image and set the text of the most recent file
                if self.dir_watcher: # check that the dir watcher exists to prevent crash
                    # set the text of the most recent file
                    self.dir_watcher.event_handler.event_path.connect(self.recent_label.setText) # might need a better label
                    # just process the image
                    if self.bin_actions[2].isChecked():
                        for im_han in self.image_handler:
                            self.dir_watcher.event_handler.event_path.connect(im_han.process)                    
            
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

    def multirun_step(self, event_path):
        """Receive event paths emitted from the system event handler signal
        for the first '# omit' events, only save the files
        then for '# hist' events, add files to a histogram,
        save the histogram 
        repeat this for the user variables in the multi-run list,
        then return to normal operation as set by the histogram binning"""
        if self.mr['v'] < np.size(self.mr['var list']):
            if self.mr['o'] < self.mr['# omit']: # don't process, just copy
                self.recent_label.setText('Just omitted: '+os.path.basename(event_path))
                self.mr['o'] += 1 # increment counter
            elif self.mr['h'] < self.mr['# hist']: # add to histogram
                # add the count to the histogram
                t1 = time.time()
                for im_han in self.image_handler:
                    im_han.process(event_path)
                t2 = time.time()
                self.int_time = t2 - t1
                # display the name of the most recent file
                self.recent_label.setText('Just processed: '+os.path.basename(event_path))
                self.plot_current_hist([x.hist_and_thresh for x in self.image_handler]) # update the displayed plot
                self.plot_time = time.time() - t2
                self.mr['h'] += 1 # increment counter

            if self.mr['o'] == self.mr['# omit'] and self.mr['h'] == self.mr['# hist']:
                self.mr['o'], self.mr['h'] = 0, 0 # reset counters
                self.var_edit.setText(str(self.mr['var list'][self.mr['v']])) # set user variable
                self.bins_text_edit(text='reset') # set histogram bins 
                success = self.update_fit()       # get best fit
                if not success:                   # if fit fails, use peak search
                    self.update_stats()
                    print(
                        '\nWarning: multi-run fit failed at ' +
                        self.mr['prefix'] + '_' + str(self.mr['v']) + '.csv')
                self.save_hist_data(
                    save_file_name=os.path.join(
                        self.multirun_save_dir.text(), self.mr['prefix']) 
                            + '_' + str(self.mr['v']) + '.csv', 
                    confirm=False)# save histogram
                for im_han in self.image_handler:
                    im_han.reset_arrays() # clear histogram
                self.mr['v'] += 1 # increment counter
            
        if self.mr['v'] == np.size(self.mr['var list']):
            self.save_varplot(
                save_file_name=os.path.join(
                    self.multirun_save_dir.text(), self.mr['prefix']) 
                        + '.dat', 
                confirm=False)# save measure file
            # reconnect previous signals to dir_watcher
            self.multirun_switch.setChecked(False) # reset multi-run button
            self.multirun_switch.setText('Start')  # reset multi-run button text
            self.set_bins() # reconnects dir_watcher with given histogram binning settings
            self.mr['o'], self.mr['h'], self.mr['v'] = 0, 0, 0 # reset counters
            self.mr['measure'] += 1 # completed a measure successfully
            self.mr['prefix'] = str(self.mr['measure']) # suggest new measure as file prefix
            self.measure_edit.setText(self.mr['prefix'])

        self.multirun_progress.setText( # update progress label
            'User variable: %s, omit %s of %s files, %s of %s histogram files, %.3g%% complete'%(
                self.mr['var list'][self.mr['v']], self.mr['o'], self.mr['# omit'],
                self.mr['h'], self.mr['# hist'], 100 * ((self.mr['# omit'] + self.mr['# hist']) * 
                self.mr['v'] + self.mr['o'] + self.mr['h']) / (self.mr['# omit'] + self.mr['# hist']) / 
                np.size(self.mr['var list'])))

    def add_stats_to_plot(self, toggle=True):
        """Take the current histogram statistics from the Histogram Statistics labels
        and add the values to the variable plot, saving the parameters to the log
        file at the same time. If any of the labels are empty, replace them with 0."""
        # append current statistics to the histogram handler's list
        for idx in range(len(self.histo_handler)):
            for key in self.histo_handler[idx].temp_vals.keys():
                self.dappend(key, self.stat_labels[self.atomX+key].text() 
                                if self.stat_labels[self.atomX+key].text() else 0)
            # append histogram stats to log file:
            with open(self.log_file_names[idx], 'a') as f:
                f.write(','.join(list(map(str, 
                        self.histo_handler[idx].temp_vals.values()))) + '\n')
        self.update_varplot_axes()  # update the plot with the new values
        self.hist_num = np.size(self.histo_handler[0].stats_dict['Hist ID'])
        
    def add_to_varplot(self, hist_han):
        """The user selects which variable they want to display on the plot
        The variables are read from the x and y axis QComboBoxes
        Then the plot is updated with statistics from the histo_handler"""
        if np.size(hist_han.vals) > 0:
            hist_han.xvals = hist_han.stats_dict[str(
                    self.plot_labels[0].currentText())] # set x values
            
            y_label = str(self.plot_labels[1].currentText())
            hist_han.yvals = hist_han.stats_dict[y_label] # set y values

            try:
                self.varplot_canvas.plot(hist_han.xvals, hist_han.yvals, 
                                    pen=None, symbol='o', symbolBrush=self.c[hist_han.i])
                # add error bars if available:
                if ('Loading probability'in y_label or 'Fidelity' in y_label
                    or 'Background peak count' in y_label or 'Signal peak count' in y_label):
                    # estimate sensible beam width at the end of the errorbar
                    if np.size(hist_han.xvals)//2:
                        beam_width = 0.1*(hist_han.xvals[1]-hist_han.xvals[0])
                    else:
                        beam_width = 0.2
                    # add widget for errorbars
                    err_bars = pg.ErrorBarItem(x=hist_han.xvals, 
                                    y=hist_han.yvals, 
                                    height=hist_han.stats_dict['Error in '+y_label],
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
        self.hist_num = 0


    #### #### save and load data functions #### ####

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


    def save_hist_data(self, trigger=None, atoms=range(2), save_file_name='', confirm=True):
        """Prompt the user to give a directory to save the histogram data, then save
        atoms specifies which histograms to save, referring to the indices of self.atomX"""
        default_path = self.get_default_path()
        try:
            if not save_file_name and 'PyQt4' in sys.modules:
                save_file_name = QFileDialog.getSaveFileName(self, 'Save File', default_path, 'csv(*.csv);;all (*)')
            elif not save_file_name and 'PyQt5' in sys.modules:
                save_file_name, _ = QFileDialog.getSaveFileName(self, 'Save File', default_path, 'csv(*.csv);;all (*)')
            self.add_stats_to_plot()
            if save_file_name:
                # don't update the threshold  - trust the user to have already set it
                for i in atoms:
                    # save separate histograms for each atom, with the atom name at the start of the file
                    self.image_handler[i].save_state(
                        os.path.join(os.path.dirname(save_file_name), 
                            self.atomX[i].replace(' ','')+os.path.basename(save_file_name)),
                        hist_header=list(self.histo_handler[i].temp_vals.keys()),
                        hist_stats=list(self.histo_handler[i].temp_vals.values())) 
                try: 
                    hist_num = self.histo_handler.stats_dict['Hist ID'][-1]
                except IndexError: # if there are no values in the stats_dict yet
                    hist_num = -1
                if confirm:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Information)
                    msg.setText("The following files were saved to directory "+os.path.dirname(save_file_name)+" \n"+
                            "\n - ".join([self.atomX[i].replace(' ','')+os.path.basename(save_file_name) for i in atoms])+
                            "\n\nand histogram %s was appended to the log files."%hist_num)
                    msg.setStandardButtons(QMessageBox.Ok)
                    msg.exec_()
                return 1
        except OSError:
            return 0 # user cancelled - file not found

    def save_varplot(self, save_file_name='', confirm=True):
        """Save the data in the current plot, which is held in the histoHandler's
        dictionary and saved in the log file, to a new file."""
        default_path = self.get_default_path('log')
        try:
            if not save_file_name and 'PyQt4' in sys.modules:
                save_file_name = QFileDialog.getSaveFileName(
                    self, 'Save File', default_path, 'dat(*.dat);;all (*)')
            elif not save_file_name and 'PyQt5' in sys.modules:
                save_file_name, _ = QFileDialog.getSaveFileName(
                    self, 'Save File', default_path, 'dat(*.dat);;all (*)')

            for idx in range(len(self.histo_handler)):
                atom_file_name = os.path.join(os.path.dirname(save_file_name), 
                            self.atomX[idx].replace(' ','')+os.path.basename(save_file_name))
                with open(atom_file_name, 'w+') as f:
                    f.write('#Single Atom Image Analyser Log File: collects histogram data\n')
                    f.write('#include --[]\n')
                    f.write('#'+', '.join(self.histo_handler[idx].stats_dict.keys())+'\n')
                    for i in range(len(self.histo_handler[idx].stats_dict['Hist ID'])):
                        f.write(','.join(list(map(str, [v[i] for v in 
                            self.histo_handler[idx].stats_dict.values()])))+'\n')
            if confirm:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText("Plot data saved to files "+", ".join([
                    os.path.join(os.path.dirname(save_file_name), 
                        X.replace(' ','')+os.path.basename(save_file_name)) 
                        for X in self.atomX]))
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
                if self.save_hist_data(atoms=idxs):
                    # only reset the histograms if the save was successful
                    for i in idxs:
                        self.image_handler[i].reset_arrays() # get rid of old data
                        self.hist_canvas[i].clear() # remove old histogram from display
            else:
                for i in idxs:
                    self.image_handler[i].reset_arrays() # get rid of old data
                    self.hist_canvas[i].clear() # remove old histogram from display
        return choice, ok, idxs

    def load_empty_hist(self):
        """Prompt the user with options to save the data and then reset the 
        histogram"""
        reply = QMessageBox.question(self, 'Confirm reset', 
            'Save the current histogram before resetting?',
            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
            QMessageBox.Cancel)
        if reply == QMessageBox.Cancel:
            return 0
        elif reply == QMessageBox.Yes:
            self.save_hist_data()  # prompt user for file name then save
            for idx in range(len(self.image_handler)):
                self.image_handler[idx].reset_arrays() # get rid of old data
                self.hist_canvas[idx].clear() # remove old histogram from display
        elif reply == QMessageBox.No:
            for idx in range(len(self.image_handler)):
                self.image_handler[idx].reset_arrays() # get rid of old data
                self.hist_canvas[idx].clear() # remove old histogram from display


    def load_from_files(self, trigger=None):
        """Prompt the user to select image files to process, then sequentially process
        them and update the histogram"""
        default_path = self.get_default_path(option='im')
        _, ok, _ = self.check_reset() # ask the user to select which atom
        if ok: # False if the user cancelled
            try:
                self.recent_label.setText('Processing files...') # comes first otherwise not executed
                if 'PyQt4' in sys.modules:
                    file_list = QFileDialog.getOpenFileNames(self, 
                        'Select Files', default_path, 'Images(*.asc);;all (*)')
                elif 'PyQt5' in sys.modules:
                    file_list, _ = QFileDialog.getOpenFileNames(self, 
                        'Select Files', default_path, 'Images(*.asc);;all (*)')
                for file_name in file_list:
                    for im_han in self.image_handler:
                        try:
                            im_han.process(file_name)
                            self.recent_label.setText(
                                'Just processed: '+os.path.basename(file_name)) # only updates at end of loop
                        except:
                            print("\n WARNING: failed to load "+file_name)
                self.update_stats()
                if self.recent_label.text == 'Processing files...':
                    self.recent_label.setText('Finished Processing')
            except OSError:
                pass # user cancelled - file not found

    def load_from_file_nums(self, trigger=None):
        """Prompt the user to enter a range of image file numbers.
        Use these to select the image files from the current image storage path.
        Sequentially process the images then update the histogram"""
        default_range = ''
        image_storage_path = self.path_label['Image Storage Path: '].text() + '\%s\%s\%s'%(self.date[3],self.date[2],self.date[0])  
        date = self.date[0]+self.date[1]+self.date[3]
        if self.image_handler[0].im_num > 0: # defualt load all files in folder
            default_range = '0 - ' + str(self.image_handler[0].im_num)
        text, ok = QInputDialog.getText( # user inputs the range
            self, 'Choose file numbers to load from','Range of file numbers: ',
            text=default_range)
        if ok and text and image_storage_path: # if user cancels or empty text, do nothing
            for file_range in text.split(','):
                minmax = file_range.split('-')
                if np.size(minmax) == 1: # only entered one file number
                    file_list = [
                        os.path.join(image_storage_path, 
                            '_' + date + '_' + minmax[0].replace(' ','') + '.asc')]
                if np.size(minmax) == 2:
                    file_list = [
                        os.path.join(image_storage_path, 
                            '_' + date + '_' + dfn + '.asc') for dfn in list(map(str, 
                            range(int(minmax[0]), int(minmax[1]))))] 
            for file_name in file_list:
                try:
                    for im_han in self.image_handler:
                        im_han.process(file_name)
                    self.recent_label.setText(
                        'Just processed: '+os.path.basename(file_name)) # only updates at end of loop
                except:
                    print("\n WARNING: failed to load "+file_name) # probably file size was wrong
            self.update_stats()
            if self.recent_label.text == 'Processing files...':
                self.recent_label.setText('Finished Processing')

    def load_from_csv(self, trigger=None):
        """Prompt the user to select a csv file to load histogram data from.
        It must have the specific layout that the image_handler saves in."""
        default_path = self.get_default_path()
        _, ok, _ = self.check_reset() # ask the user to select which atom
        if ok:
            try:
                for im_han in self.image_handler: # load separate csv files for each atom
                    if 'PyQt4' in sys.modules: 
                        file_name = QFileDialog.getOpenFileName(self, 'Select File for '+im_han.X, 
                                                            default_path, 'csv(*.csv);;all (*)')
                    elif 'PyQt5' in sys.modules:
                        file_name, _ = QFileDialog.getOpenFileName(self, 'Select File for '+im_han.X, 
                                                            default_path, 'csv(*.csv);;all (*)')
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
            for hist_han in self.histo_handler:
                if 'PyQt4' in sys.modules: 
                    file_name = QFileDialog.getOpenFileName(self, 'Select File for '+hist_han.X, 
                                                            default_path, 'dat(*.dat);;all (*)')
                elif 'PyQt5' in sys.modules:
                    file_name, _ = QFileDialog.getOpenFileName(self, 'Select File for '+hist_han.X, 
                                                            default_path, 'dat(*.dat);;all (*)')

                    hist_han.load_from_log(file_name)
            self.update_varplot_axes()
        except OSError:
            pass # user cancelled - file not found


    #### #### testing functions #### #### 
        
    def print_times(self, unit="s"):
        """Display the times measured for functions"""
        scale = 1
        if unit == "ms" or unit == "milliseconds":
            scale *= 1e3
        elif unit == "us" or unit == "microseconds":
            scale *= 1e6
        else:
            unit = "s"
        if self.dir_watcher: # this is None if dir_watcher isn't initiated
            print("\nFile processing event duration: %.4g "%(
                self.dir_watcher.event_handler.event_t*scale)+unit)
            print("Most recent idle time between events: %.4g "%(
                self.dir_watcher.event_handler.idle_t*scale)+unit)
            print("Image processing duration: %.4g "%(
                self.int_time*scale)+unit)
            print("Image plotting duration: %.4g "%(
                self.plot_time*scale)+unit)
            print("File writing duration: %.4g "%(
                self.dir_watcher.event_handler.write_t*scale)+unit)
            print("File copying duration: %.4g "%(
                self.dir_watcher.event_handler.copy_t*scale)+unit)
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