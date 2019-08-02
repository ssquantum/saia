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
            QTabWidget, QVBoxLayout, QFont) 
except ImportError:
    from PyQt5.QtCore import QThread, pyqtSignal, QEvent
    from PyQt5.QtGui import (QGridLayout, QMessageBox, QLineEdit, QIcon, 
            QFileDialog, QDoubleValidator, QIntValidator, QComboBox, QMenu, 
            QActionGroup, QVBoxLayout, QFont)
    from PyQt5.QtWidgets import (QApplication, QPushButton, QWidget, QTabWidget,
        QAction, QMainWindow, QLabel)
          
####    ####    ####    ####

# main GUI window contains all the widgets                
class main_window(QMainWindow):
    """Use Qt to produce the window where the histogram plot is shown
    have a simple interface for the user to control the limits of the plot 
    and number of bins in the histogram
    
    This GUI was produced with help from http://zetcode.com/gui/pyqt5/"""
    def __init__(self):
        super().__init__()
        self.bias = 698.5 # bias off set from EMCCD
        self.Nr   = 8.8   # read-out noise from EMCCD
        self.dir_watcher = None  # a button will initiate the dir watcher
        self.image_handler = ih.image_handler() # class to process images
        self.histo_handler = hh.histo_handler() # class to process histograms
        self.hist_num = 0 # ID number for the next histogram 
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
        dir_watcher_dict = dw.dir_watcher.get_dirs(self.config_edit.text()) # static method
        # make subdirectory if it doesn't already exist
        log_file_dir = dir_watcher_dict['Log File Path: ']+'\%s\%s\%s'%(self.date[3],self.date[2],self.date[0]) 
        try:
            os.makedirs(log_file_dir, exist_ok=True)
        except PermissionError:  # couldn't access the path, start a log file here
            log_file_dir = '.\%s\%s\%s'%(self.date[3],self.date[2],self.date[0])
            os.makedirs(log_file_dir, exist_ok=True)

        # log is saved in a dated subdirectory and the file name also has the date
        self.log_file_name = os.path.join(log_file_dir, 
                   'log'+self.date[0]+self.date[1]+self.date[3]+'.dat')  
        # write the header to the log file
        if not os.path.isfile(self.log_file_name): # don't overwrite if it already exists
            with open(self.log_file_name, 'w+') as f:
                f.write('#Single Atom Image Analyser Log File: collects histogram data\n')
                f.write('#include --[]\n')
                f.write('#'+', '.join(self.histo_handler.stats_dict.keys())+'\n')
       

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
        reset_hist.triggered.connect(self.load_empty_hist)
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

        # load plots from log files
        varplot_menu = menubar.addMenu('Plotting')

        load_varplot = QAction('Load from log file', self)
        load_varplot.triggered.connect(self.load_from_log)
        varplot_menu.addAction(load_varplot)

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
        self.pic_size_edit.setText(str(self.image_handler.pic_size)) # default
        self.pic_size_edit.textChanged[str].connect(self.pic_size_text_edit)
        self.pic_size_edit.setValidator(double_validator)

        # get image size from loading an image
        load_im_size = QPushButton('Load size from image', self)
        load_im_size.clicked.connect(self.load_im_size) # load image size from image
        load_im_size.resize(load_im_size.sizeHint())
        settings_grid.addWidget(load_im_size, 0,2, 1,1)

        # get ROI centre from loading an image
        load_roi = QPushButton('Get ROI from image', self)
        load_roi.clicked.connect(self.load_roi) # load roi centre from image
        load_roi.resize(load_im_size.sizeHint())
        settings_grid.addWidget(load_roi, 1,2, 1,1)

        # get user to set ROI:
        # centre of ROI x position
        roi_xc_label = QLabel('ROI x_c: ', self)
        settings_grid.addWidget(roi_xc_label, 1,0, 1,1)
        self.roi_x_edit = QLineEdit(self)
        settings_grid.addWidget(self.roi_x_edit, 1,1, 1,1)
        self.roi_x_edit.setText('0')  # default
        self.roi_x_edit.textEdited[str].connect(self.roi_text_edit)
        self.roi_x_edit.setValidator(int_validator) # only numbers
        
        # centre of ROI y position
        roi_yc_label = QLabel('ROI y_c: ', self)
        settings_grid.addWidget(roi_yc_label, 2,0, 1,1)
        self.roi_y_edit = QLineEdit(self)
        settings_grid.addWidget(self.roi_y_edit, 2,1, 1,1)
        self.roi_y_edit.setText('0')  # default
        self.roi_y_edit.textEdited[str].connect(self.roi_text_edit)
        self.roi_y_edit.setValidator(int_validator) # only numbers
        
        # ROI size
        roi_l_label = QLabel('ROI size: ', self)
        settings_grid.addWidget(roi_l_label, 3,0, 1,1)
        self.roi_l_edit = QLineEdit(self)
        settings_grid.addWidget(self.roi_l_edit, 3,1, 1,1)
        self.roi_l_edit.setText('1')  # default
        self.roi_l_edit.textEdited[str].connect(self.roi_text_edit)
        self.roi_l_edit.setValidator(int_validator) # only numbers

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
        self.config_edit.setText('./config.dat')
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
        self.dw_mode = QPushButton('Active', self, checkable=True)
        self.dw_mode.setChecked(True)
        self.dw_mode.clicked[bool].connect(self.dw_mode_switch)
        settings_grid.addWidget(self.dw_mode, i+8,1, 1,1)

        # label to show status of dir watcher
        self.dw_status_label = QLabel('Stopped', self)  # should errors stop dir watcher???
        settings_grid.addWidget(self.dw_status_label, i+8,2, 1,1)

        # label to show last file analysed
        self.recent_label = QLabel('', self)
        settings_grid.addWidget(self.recent_label, i+9,0, 1,4)
        

        #### tab for histogram ####
        hist_tab = QWidget()
        hist_grid = QGridLayout()
        hist_tab.setLayout(hist_grid)
        self.tabs.addTab(hist_tab, "Histogram")

        # main subplot of histogram
        self.hist_canvas = pg.PlotWidget()
        self.hist_canvas.getAxis('bottom').tickFont = font
        self.hist_canvas.getAxis('left').tickFont = font
        self.hist_canvas.setTitle("Histogram of CCD counts")
        hist_grid.addWidget(self.hist_canvas, 1,0, 6,8)  # allocate space in the grid
        
        # adjustable parameters: min/max counts, number of bins
        # min counts:
        min_counts_label = QLabel('Min. Counts: ', self)
        hist_grid.addWidget(min_counts_label, 0,0, 1,1)
        self.min_counts_edit = QLineEdit(self)
        hist_grid.addWidget(self.min_counts_edit, 0,1, 1,1)
        self.min_counts_edit.textChanged[str].connect(self.bins_text_edit)
        self.min_counts_edit.setValidator(double_validator)
        
        # max counts:
        max_counts_label = QLabel('Max. Counts: ', self)
        hist_grid.addWidget(max_counts_label, 0,2, 1,1)
        self.max_counts_edit = QLineEdit(self)
        hist_grid.addWidget(self.max_counts_edit, 0,3, 1,1)
        self.max_counts_edit.textChanged[str].connect(self.bins_text_edit)
        self.max_counts_edit.setValidator(double_validator)
        
        # number of bins
        num_bins_label = QLabel('# Bins: ', self)
        hist_grid.addWidget(num_bins_label, 0,4, 1,1)
        self.num_bins_edit = QLineEdit(self)
        hist_grid.addWidget(self.num_bins_edit, 0,5, 1,1)
        self.num_bins_edit.textChanged[str].connect(self.bins_text_edit)
        self.num_bins_edit.setValidator(double_validator)

        # user can set the threshold
        self.thresh_toggle = QPushButton('User Threshold: ', self)
        self.thresh_toggle.setCheckable(True)
        self.thresh_toggle.clicked[bool].connect(self.set_thresh)
        hist_grid.addWidget(self.thresh_toggle, 0,6, 1,1)
        # user inputs threshold
        self.thresh_edit = QLineEdit(self)
        hist_grid.addWidget(self.thresh_edit, 0,7, 1,1)
        self.thresh_edit.textChanged[str].connect(self.bins_text_edit)
        self.thresh_edit.setValidator(double_validator)
        
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
        for i, label_text in enumerate(self.histo_handler.stats_dict.keys()):
            new_label = QLabel(label_text, self) # description
            stat_grid.addWidget(new_label, i+1,0, 1,1)
            self.stat_labels[label_text] = QLabel('', self) # value
            stat_grid.addWidget(self.stat_labels[label_text], i+1,1, 1,1)
            
        # update statistics
        stat_update = QPushButton('Update statistics', self)
        stat_update.clicked[bool].connect(self.update_stats)
        stat_grid.addWidget(stat_update, i+2,0, 1,1)

        # do Gaussian/Poissonian fit - peaks and widths
        fit_update = QPushButton('Get best fit', self)
        fit_update.clicked[bool].connect(self.update_fit)
        stat_grid.addWidget(fit_update, i+2,1, 1,1)

        # do a Gaussian fit just to the background peak
        fit_bg = QPushButton('Fit background', self)
        fit_bg.clicked[bool].connect(self.fit_bg_gaussian)
        stat_grid.addWidget(fit_bg, i+2,2, 1,1)

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
        self.pic_size_label.setText(str(self.image_handler.pic_size)) # default

        # toggle to continuously plot images as they come in
        self.im_show_toggle = QPushButton('Auto-display last image', self)
        self.im_show_toggle.setCheckable(True)
        self.im_show_toggle.clicked[bool].connect(self.set_im_show)
        im_grid.addWidget(self.im_show_toggle, 0,2, 1,1)
        
        im_grid_pos = 0 # starting column. 
        # centre of ROI x position
        self.xc_label = QLabel('ROI x_c: 0', self)
        im_grid.addWidget(self.xc_label, 7,im_grid_pos, 1,1)
        
        # centre of ROI y position
        self.yc_label = QLabel('ROI y_c: 0', self)
        im_grid.addWidget(self.yc_label, 7,im_grid_pos+2, 1,1)
        
        # ROI size
        self.l_label = QLabel('ROI size: 1', self)
        im_grid.addWidget(self.l_label, 7,im_grid_pos+4, 1,1)
        
        # display last image if toggle is True
        im_widget = pg.GraphicsLayoutWidget() # containing widget
        viewbox = im_widget.addViewBox() # plot area to display image
        self.im_canvas = pg.ImageItem() # the image
        viewbox.addItem(self.im_canvas)
        im_grid.addWidget(im_widget, 1,im_grid_pos, 6,8)
        # make an ROI that the user can drag
        self.roi = pg.ROI([0,0], [1,1]) 
        self.roi.addScaleHandle([1,1], [0.5,0.5]) # allow user to adjust ROI size
        viewbox.addItem(self.roi)
        self.roi.setZValue(10)   # make sure the ROI is drawn above the image
        self.roi.sigRegionChangeFinished.connect(self.user_roi) # signal emitted when user stops dragging ROI
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
            self.plot_labels[i].addItems(list(self.histo_handler.stats_dict.keys())) # add options
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

        #### choose main window position and dimensions: (xpos,ypos,width,height)
        self.setGeometry(100, 100, 850, 700)
        self.setWindowTitle('Single Atom Image Analyser')
        self.setWindowIcon(QIcon('docs/tempicon.png'))
        
    #### #### initiation functions #### #### 

    def init_DW(self):
        """Ask the user if they want to start the dir watcher or not"""
        dir_watcher_dict = dw.dir_watcher.get_dirs(self.config_edit.text()) # static method
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
            # prompt user if they want to remove image files
            self.dir_watcher = dw.dir_watcher(config_file=self.config_edit.text(),
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
            msg.setText("Directory Watcher initiated in " + self.dw_mode.text()
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
        self.histo_handler.temp_vals['User variable'] = self.var_edit.text()
        self.stat_labels['User variable'].setText(self.var_edit.text())

    def path_text_edit(self, text=''):
        """The user finishes editing an edit text box by pressing return or clicking
        somewhere else, then the text is sent to this function to reload the config file. 
        The dir watcher is not updated unless the 'initiate dir watcher' button is used."""
        dw_dict = dw.dir_watcher.get_dirs(self.config_edit.text())
        for key, value in dw_dict.items():
            self.path_label[key].setText(value)
    
    def user_roi(self, pos):
        """The user drags an ROI and this updates the ROI centre and width"""
        x0, y0 = self.roi.pos()  # lower left corner of bounding rectangle
        xw, yw = self.roi.size() # widths
        l = int(0.5*(xw+yw))  # want a square ROI
        # note: setting the origin as bottom left but the image has origin top left
        xc, yc = int(x0 + l//2), int(y0 + l//2)  # centre
        self.image_handler.set_roi(dimensions=[xc, yc, l])
        self.xc_label.setText('ROI x_c = '+str(xc)) 
        self.yc_label.setText('ROI y_c = '+str(yc))
        self.l_label.setText('ROI size = '+str(l))
        self.roi_x_edit.setText(str(xc))
        self.roi_y_edit.setText(str(yc))
        self.roi_l_edit.setText(str(l))
            
    def pic_size_text_edit(self, text):
        """Update the specified size of an image in pixels when the user 
        edits the text in the line edit widget"""
        self.image_handler.pic_size = int(text)
        self.pic_size_label.setText(str(self.image_handler.pic_size))

    def CCD_stat_edit(self):
        """Update the values used for the EMCCD bias offset and readout noise"""
        if self.bias_offset_edit.text(): # check the label isn't empty
            self.bias = float(self.bias_offset_edit.text())
        if self.read_noise_edit.text():
            self.Nr = float(self.read_noise_edit.text())
        
    def roi_text_edit(self, text):
        """Update the ROI position and size every time a text edit is made by
        the user to one of the line edit widgets"""
        xc, yc, l = [self.roi_x_edit.text(),
                            self.roi_y_edit.text(), self.roi_l_edit.text()]
        if any([v == '' for v in [xc, yc, l]]):
            xc, yc, l = 0, 0, 1 # default takes the top left pixel
        else:
            xc, yc, l = list(map(int, [xc, yc, l])) # crashes if the user inputs float
        
        if (xc - l//2 < 0 or yc - l//2 < 0 
            or xc + l//2 > self.image_handler.pic_size 
            or yc + l//2 > self.image_handler.pic_size):
            l = 2*min([xc, yc])  # can't have the boundary go off the edge
        if int(l) == 0:
            l = 1 # can't have zero width
        
        self.image_handler.set_roi(dimensions=list(map(int, [xc, yc, l])))
        self.xc_label.setText('ROI x_c = '+str(xc)) 
        self.yc_label.setText('ROI y_c = '+str(yc))
        self.l_label.setText('ROI size = '+str(l))
        # update ROI on image canvas
        # note: setting the origin as top left because image is inverted
        self.roi.setPos(xc - l//2, yc - l//2)
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
            if int(new_vals[2]) < 2:
                # 0 bins causes value error
                new_vals[2] = 10
            min_bin, max_bin, num_bins = list(map(int, new_vals))
            
            # set the new values for the bins of the image handler
            self.image_handler.bin_array = np.linspace(min_bin, max_bin, num_bins)

            # set the new threshold if supplied
            if self.thresh_toggle.isChecked():
                try:
                    self.image_handler.thresh = float(self.thresh_edit.text())
                    self.stat_labels['Threshold'].setText(str(int(self.image_handler.thresh)))
                except ValueError: pass # user switched toggle before inputing text
                self.plot_current_hist(self.image_handler.histogram) # doesn't update thresh
            else:
                self.plot_current_hist(self.image_handler.hist_and_thresh) # updates thresh
            
    
    #### #### toggle functions #### #### 

    def dappend(self, key, value):
        """Shorthand for appending a new value to the array in the histo_handler
        dictionary with the given key. Also update the temp values so that they
        are always the most recent value. Then update the label in the statistics
        tab to display the new value."""
        self.histo_handler.stats_dict[key] = np.append(
                self.histo_handler.stats_dict[key], 
                np.array(value).astype(self.histo_handler.stats_dict[key].dtype))
        self.histo_handler.temp_vals[key] = value
        self.stat_labels[key].setText(str(value))

    def update_stats(self, toggle=True):
        """Update the statistics from the current histogram in order to save them
        image_handler uses a peak finding algorithm to get the peak positions and widths
        from these a threshold can be set. If the user thresh_toggle is checked then the
        threshold will not be updated.
        The histo_handler stores temporary values that we might not yet want to add to
        the plot."""
        if self.image_handler.im_num > 0: # only update if a histogram exists
            if self.thresh_toggle.isChecked(): # using manual threshold
                self.plot_current_hist(self.image_handler.histogram) # update hist and peak stats, keep thresh
            else:
                self.plot_current_hist(self.image_handler.hist_and_thresh) # update hist and get peak stats

            above_idxs = np.where(self.image_handler.atom > 0)[0] # index of images with counts above threshold
            atom_count = np.size(above_idxs)  # number of images with counts above threshold
            above = self.image_handler.counts[above_idxs] # counts above threshold
            below_idxs = np.where(self.image_handler.atom[:self.image_handler.im_num] == 0)[0] # index of images with counts below threshold
            empty_count = np.size(below_idxs) # number of images with counts below threshold
            below = self.image_handler.counts[below_idxs] # counts below threshold
            # use the binomial distribution to get 1 sigma confidence intervals:
            conf = binom_conf_interval(atom_count, atom_count + empty_count, interval='jeffreys') 

            # store the calculated histogram statistics as temp, don't add to plot
            self.histo_handler.temp_vals['Hist ID'] = int(self.hist_num)
            self.histo_handler.temp_vals['User variable'] = float(self.var_edit.text())
            self.histo_handler.temp_vals['Number of images processed'] = self.image_handler.im_num
            self.histo_handler.temp_vals['Counts above : below threshold'] = str(atom_count) + ' : ' + str(empty_count)
            self.histo_handler.temp_vals['Loading probability'] = np.around(atom_count/self.image_handler.im_num, 4)
            self.histo_handler.temp_vals['Error in Loading probability'] = np.around(conf[1] - conf[0], 4)
            if np.size(self.image_handler.peak_counts) == 2:
                self.histo_handler.temp_vals['Background peak count'] = int(self.image_handler.peak_counts[0])
                # assume bias offset is self.bias, readout noise standard deviation Nr
                if self.Nr**2+self.image_handler.peak_counts[0]-self.bias > 0:
                    self.histo_handler.temp_vals['sqrt(Nr^2 + Nbg)'] = int((self.Nr**2+self.image_handler.peak_counts[0]-self.bias)**0.5)
                else: # don't take the sqrt of a -ve number
                    self.histo_handler.temp_vals['sqrt(Nr^2 + Nbg)'] = 0
                self.histo_handler.temp_vals['Background peak width'] = int(self.image_handler.peak_widths[0])
                self.histo_handler.temp_vals['Error in Background peak count'] = np.around(self.image_handler.peak_widths[0] / empty_count**0.5, 2)
                self.histo_handler.temp_vals['Background mean'] = np.around(np.mean(below), 1)
                self.histo_handler.temp_vals['Background standard deviation'] = np.around(np.std(below, ddof=1), 1)
                self.histo_handler.temp_vals['Signal peak count'] = int(self.image_handler.peak_counts[1])
                # assume bias offset is self.bias, readout noise standard deviation Nr
                if self.Nr**2+self.image_handler.peak_counts[1]-self.bias > 0:
                    self.histo_handler.temp_vals['sqrt(Nr^2 + Ns)'] = int((self.Nr**2+self.image_handler.peak_counts[1]-self.bias)**0.5)
                else: # don't take the sqrt of a -ve number
                    self.histo_handler.temp_vals['sqrt(Nr^2 + Ns)'] = 0
                self.histo_handler.temp_vals['Signal peak width'] = int(self.image_handler.peak_widths[1])
                self.histo_handler.temp_vals['Error in Signal peak count'] = np.around(self.image_handler.peak_widths[1] / atom_count**0.5, 2)
                self.histo_handler.temp_vals['Signal mean'] = np.around(np.mean(above), 1)
                self.histo_handler.temp_vals['Signal standard deviation'] = np.around(np.std(above, ddof=1), 1)
                self.histo_handler.temp_vals['Separation'] = int(self.image_handler.peak_counts[1] - 
                                                                        self.image_handler.peak_counts[0])
                self.histo_handler.temp_vals['Error in Separation'] = np.around(np.sqrt(self.image_handler.peak_widths[0]**2 / empty_count
                                + self.image_handler.peak_widths[1]**2 / atom_count), 2)
                self.histo_handler.temp_vals['Fidelity'] = self.image_handler.fidelity
                self.histo_handler.temp_vals['Error in Fidelity'] = self.image_handler.err_fidelity
            else:
                for key in ['Background peak count', 'sqrt(Nr^2 + Nbg)', 'Background peak width', 
                'Error in Background peak count', 'Signal peak count', 'sqrt(Nr^2 + Ns)', 
                'Signal peak width', 'Error in Signal peak count', 'Separation', 'Fidelity', 'Error in Fidelity']:
                    self.histo_handler.temp_vals[key] = 0

            self.histo_handler.temp_vals['Threshold'] = int(self.image_handler.thresh)

            # display the new statistics in the labels
            for key, val in self.histo_handler.temp_vals.items():
                self.stat_labels[key].setText(str(val))

    def fit_gaussians(self, store_stats=False):
        """Update the histogram and fit two Gaussians, splitting the data at the threshold
        then use the fits to calculate histogram statistics, and set the threshold where the 
        fidelity is maximum. If the store_stats Boolean is True, append the calculated values
        to the histo_handler's statistics dictionary."""
        bins, occ, thresh = self.image_handler.histogram()  # get histogram
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
        self.image_handler.peak_heights = np.array((best_fits[0].ps[0], best_fits[1].ps[0]))
        self.image_handler.peak_counts = np.array((best_fits[0].ps[1], best_fits[1].ps[1]))
        self.image_handler.peak_widths = np.array((best_fits[0].ps[2], best_fits[1].ps[2]))

        # update threshold to where fidelity is maximum
        if not self.thresh_toggle.isChecked(): # update thresh if not set by user
            self.image_handler.search_fidelity(best_fits[0].ps[1], best_fits[1].ps[1], n=100)
        else:
            self.image_handler.fidelity, self.image_handler.err_fidelity = np.around(
                            self.image_handler.get_fidelity(), 4) # round to 4 d.p.

        self.plot_current_hist(self.image_handler.histogram) # clear then update histogram plot
        for bf in best_fits:
            xs = np.linspace(min(bf.x), max(bf.x), 100) # interpolate
            self.hist_canvas.plot(xs, bf.gauss(xs, *bf.ps), pen='b') # plot best fit

        # update atom statistics
        self.image_handler.atom[:self.image_handler.im_num] = self.image_handler.counts[
            :self.image_handler.im_num] // self.image_handler.thresh   # update atom presence
        above_idxs = np.where(self.image_handler.atom > 0)[0] # index of images with counts above threshold
        atom_count = np.size(above_idxs)  # number of images with counts above threshold
        above = self.image_handler.counts[above_idxs] # counts above threshold
        below_idxs = np.where(self.image_handler.atom[:self.image_handler.im_num] == 0)[0] # index of images with counts below threshold
        empty_count = np.size(below_idxs) # number of images with counts below threshold
        below = self.image_handler.counts[below_idxs] # counts below threshold
        loading_prob = np.around(atom_count/self.image_handler.im_num, 4) # loading probability
        # use the binomial distribution to get 1 sigma confidence intervals:
        conf = binom_conf_interval(atom_count, atom_count + empty_count, interval='jeffreys') 
        loading_err = np.around(conf[1] - conf[0], 4)

        # store the calculated histogram statistics as temp, don't add to plot
        if store_stats:
            self.histo_handler.temp_vals['Hist ID'] = int(self.hist_num)
            self.histo_handler.temp_vals['User variable'] = float(self.var_edit.text())
            self.histo_handler.temp_vals['Number of images processed'] = self.image_handler.im_num
            self.histo_handler.temp_vals['Counts above : below threshold'] = str(atom_count) + ' : ' + str(empty_count)
            self.histo_handler.temp_vals['Loading probability'] = loading_prob
            self.histo_handler.temp_vals['Error in Loading probability'] = loading_err
            self.histo_handler.temp_vals['Background peak count'] = int(best_fits[0].ps[1])
            # assume bias offset is self.bias, readout noise standard deviation Nr
            if self.Nr**2+best_fits[0].ps[1]-self.bias > 0:
                self.histo_handler.temp_vals['sqrt(Nr^2 + Nbg)'] = int((self.Nr**2+best_fits[0].ps[1]-self.bias)**0.5)
            else: # don't take the sqrt of a -ve number
                self.histo_handler.temp_vals['sqrt(Nr^2 + Nbg)'] = 0
            self.histo_handler.temp_vals['Background peak width'] = int(best_fits[0].ps[2])
            self.histo_handler.temp_vals['Error in Background peak count'] = np.around(best_fits[0].ps[2] / empty_count**0.5, 2)
            self.histo_handler.temp_vals['Background mean'] = np.around(np.mean(below), 1)
            self.histo_handler.temp_vals['Background standard deviation'] = np.around(np.std(below, ddof=1), 1)
            self.histo_handler.temp_vals['Signal peak count'] = int(best_fits[1].ps[1])
            # assume bias offset is self.bias, readout noise standard deviation Nr
            if self.Nr**2+best_fits[1].ps[1]-self.bias > 0:
                self.histo_handler.temp_vals['sqrt(Nr^2 + Ns)'] = int((self.Nr**2+best_fits[1].ps[1]-self.bias)**0.5)
            else:
                self.histo_handler.temp_vals['sqrt(Nr^2 + Ns)'] = 0
            self.histo_handler.temp_vals['Signal peak width'] = int(best_fits[1].ps[2])
            self.histo_handler.temp_vals['Error in Signal peak count'] = np.around(best_fits[1].ps[2] / atom_count**0.5, 2)
            self.histo_handler.temp_vals['Signal mean'] = np.around(np.mean(above), 1)
            self.histo_handler.temp_vals['Signal standard deviation'] = np.around(np.std(above, ddof=1), 1)
            self.histo_handler.temp_vals['Separation'] = int(best_fits[1].ps[1] - best_fits[0].ps[1])
            self.histo_handler.temp_vals['Error in Separation'] = np.around(np.sqrt(best_fits[0].ps[2]**2 / empty_count
                        + best_fits[1].ps[2]**2 / atom_count), 2)
            self.histo_handler.temp_vals['Fidelity'] = self.image_handler.fidelity
            self.histo_handler.temp_vals['Error in Fidelity'] = self.image_handler.err_fidelity
            self.histo_handler.temp_vals['Threshold'] = int(self.image_handler.thresh)

            # display the new statistics in the labels
            for key, val in self.histo_handler.temp_vals.items():
                self.stat_labels[key].setText(str(val))
        
        return 1 # fit successful

    def update_fit(self, toggle=True):
        """Fit Gaussians to the peaks and use it to get a better estimate of the 
        peak centres and widths. The peaks are not quite Poissonian because of the
        bias from dark counts. Use the fits to get histogram statistics, then set 
        the threshold to maximise fidelity. Iterate until the threshold converges."""
        if self.image_handler.im_num > 0: # only update if a histogram exists
            oldthresh = self.image_handler.thresh # store the last value
            diff = 1                              # convergence criterion
            for i in range(20):          # shouldn't need many iterations
                if diff < 0.0015:
                    break
                success = self.fit_gaussians()
                diff = abs(oldthresh - self.image_handler.thresh) / float(oldthresh)

            if success: # fit_gaussians returns 0 if it fails
                self.fit_gaussians(store_stats=True) # add new stats to histo_handler
                
            
    def fit_bg_gaussian(self, store_stats=False):
        """Assume that there is only one peak in the histogram as there is no single
        atom signal. Fit a Gaussian to this peak."""
        n = self.image_handler.im_num # number of images processed
        c = self.image_handler.counts[:n] # integrated counts

        bins, occ, _ = self.image_handler.histogram()  # get histogram
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

        # update image handler's values for peak parameters
        self.image_handler.peak_heights = np.array((best_fit.ps[0], 0))
        self.image_handler.peak_counts = np.array((best_fit.ps[1], 0))
        self.image_handler.peak_widths = np.array((best_fit.ps[2], 0))

        self.plot_current_hist(self.image_handler.histogram) # clear then update histogram plot
        
        xs = np.linspace(min(best_fit.x), max(best_fit.x), 100) # interpolate
        self.hist_canvas.plot(xs, best_fit.gauss(xs, *best_fit.ps), pen='b') # plot best fit

        # store the calculated histogram statistics as temp, don't add to plot
        self.histo_handler.temp_vals['Hist ID'] = int(self.hist_num)
        self.histo_handler.temp_vals['User variable'] = float(self.var_edit.text())
        self.histo_handler.temp_vals['Number of images processed'] = n
        self.histo_handler.temp_vals['Counts above : below threshold'] = '0 : ' + str(n)
        self.histo_handler.temp_vals['Loading probability'] = 0
        self.histo_handler.temp_vals['Error in Loading probability'] = 0
        self.histo_handler.temp_vals['Background peak count'] = int(best_fit.ps[1])
        # assume bias offset is self.bias, readout noise standard deviation Nr
        if self.Nr**2+mu-self.bias:
            self.histo_handler.temp_vals['sqrt(Nr^2 + Nbg)'] = int((self.Nr**2+mu-self.bias)**0.5)
        else: # don't take the sqrt of a -ve number
            self.histo_handler.temp_vals['sqrt(Nr^2 + Nbg)'] = 0
        self.histo_handler.temp_vals['Background peak width'] = int(best_fit.ps[2])
        self.histo_handler.temp_vals['Error in Background peak count'] = np.around(best_fit.ps[2] / n**0.5, 4)
        self.histo_handler.temp_vals['Background mean'] = np.around(mu, 1)
        self.histo_handler.temp_vals['Background standard deviation'] = np.around(sig, 1)
        self.histo_handler.temp_vals['Signal peak count'] = 0
        self.histo_handler.temp_vals['sqrt(Nr^2 + Ns)'] = 0
        self.histo_handler.temp_vals['Signal peak width'] = 0
        self.histo_handler.temp_vals['Error in Signal peak count'] = 0
        self.histo_handler.temp_vals['Signal mean'] = 0
        self.histo_handler.temp_vals['Signal standard deviation'] = 0
        self.histo_handler.temp_vals['Separation'] = 0
        self.histo_handler.temp_vals['Error in Separation'] = 0
        self.histo_handler.temp_vals['Fidelity'] = 0
        self.histo_handler.temp_vals['Error in Fidelity'] = 0
        self.histo_handler.temp_vals['Threshold'] = int(self.image_handler.thresh)

        # display the new statistics in the labels
        for key, val in self.histo_handler.temp_vals.items():
            self.stat_labels[key].setText(str(val))
            
    
    def update_varplot_axes(self, label=''):
        """The user selects which variable they want to display on the plot
        The variables are read from the x and y axis QComboBoxes
        Then the plot is updated"""
        if np.size(self.histo_handler.stats_dict['Hist ID']) > 0:
            self.histo_handler.xvals = self.histo_handler.stats_dict[
                                str(self.plot_labels[0].currentText())] # set x values
            
            y_label = str(self.plot_labels[1].currentText())
            self.histo_handler.yvals = self.histo_handler.stats_dict[y_label] # set y values
            
            self.varplot_canvas.clear()  # remove previous data
            try:
                self.varplot_canvas.plot(self.histo_handler.xvals, 
                            self.histo_handler.yvals, pen=None, symbol='o')
                # add error bars if available:
                if ('Loading probability' in y_label or 'Fidelity' in y_label 
        or 'Background peak count' in y_label or 'Signal peak count' in y_label):
                    # add widget for errorbars
                    # estimate sensible beam width at the end of the errorbar
                    if np.size(self.histo_handler.xvals)//2:
                        beam_width = 0.1*(self.histo_handler.xvals[1]
                                                - self.histo_handler.xvals[0])
                    else:
                        beam_width = 0.2
                    err_bars = pg.ErrorBarItem(x=self.histo_handler.xvals, 
                        y=self.histo_handler.yvals, 
                        height=self.histo_handler.stats_dict['Error in '+y_label],
                        beam=beam_width)
                    self.varplot_canvas.addItem(err_bars)
            except Exception: pass # probably wrong length of arrays

    def clear_varplot(self):
        """Clear the plot of histogram statistics by resetting the histo_handler.
        The data is not lost since it has been appended to the log file."""
        self.histo_handler.__init__ () # empty the stored arrays
        self.varplot_canvas.clear()    # clear the displayed plot
        self.hist_num = 0


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
            self.image_handler.bin_array = []
            if self.thresh_toggle.isChecked():
                self.plot_current_hist(self.image_handler.histogram)
            else:
                self.plot_current_hist(self.image_handler.hist_and_thresh)

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
        self.hist_canvas.plot(bins, occ, stepMode=True, pen='k',
                                fillLevel=0, brush = (220,220,220,220)) # histogram
        self.hist_canvas.plot([thresh]*2, [0, max(occ)], pen='r') # threshold line

    
    def update_im(self, event_path):
        """Receive the event path emitted from the system event handler signal
        display the image from the file in the image canvas"""
        im_vals = self.image_handler.load_full_im(event_path)
        self.im_canvas.setImage(im_vals)
        self.im_hist.setLevels(np.min(im_vals), np.max(im_vals))
        
        
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

    def add_stats_to_plot(self, toggle=True):
        """Take the current histogram statistics from the Histogram Statistics labels
        and add the values to the variable plot, saving the parameters to the log
        file at the same time. If any of the labels are empty, do nothing and return -1"""
        # append current statistics to the histogram handler's list
        for key in self.stat_labels.keys():
            self.dappend(key, self.stat_labels[key].text())
        
        self.update_varplot_axes()  # update the plot with the new values
        self.hist_num = np.size(self.histo_handler.stats_dict['Hist ID']) # index for histograms
        # append histogram stats to log file:
        with open(self.log_file_name, 'a') as f:
            f.write(','.join(self.histo_handler.temp_vals.values()) + '\n')


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
            default_path = os.path.dirname(self.log_file_name)
        
        return default_path


    def load_im_size(self):
        """Get the user to select an image file and then use this to get the image size"""
        default_path = self.get_default_path(option='im')
        try:
            if 'PyQt4' in sys.modules:
                file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')
            elif 'PyQt5' in sys.modules:
                file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'Images (*.asc);;all (*)')

            self.image_handler.set_pic_size(file_name) # sets image handler's pic size
            self.pic_size_edit.setText(str(self.image_handler.pic_size)) # update loaded value
            self.pic_size_label.setText(str(self.image_handler.pic_size)) # update loaded value

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
            self.image_handler.set_pic_size(file_name) # sets image handler's pic size
            self.pic_size_edit.setText(str(self.image_handler.pic_size)) # update loaded value
            self.pic_size_label.setText(str(self.image_handler.pic_size)) # update loaded value
            # get the position of the max count
            self.image_handler.set_roi(im_name=file_name) # sets xc and yc
            self.roi_x_edit.setText(str(self.image_handler.xc)) # update loaded value
            self.roi_y_edit.setText(str(self.image_handler.yc)) 
            self.roi_l_edit.setText(str(self.image_handler.roi_size))
            self.xc_label.setText(str(self.image_handler.xc))
            self.yc_label.setText(str(self.image_handler.yc))
            self.l_label.setText(str(self.image_handler.roi_size))
            self.roi.setPos(self.image_handler.xc - self.image_handler.roi_size//2, 
            self.image_handler.yc - self.image_handler.roi_size//2) # set ROI in image display
            self.roi.setSize(self.image_handler.roi_size, self.image_handler.roi_size)

        except OSError:
            pass # user cancelled - file not found


    def save_hist_data(self, trigger=None):
        """Prompt the user to give a directory to save the histogram data, then save"""
        default_path = self.get_default_path()
        try:
            if 'PyQt4' in sys.modules:
                save_file_name = QFileDialog.getSaveFileName(self, 'Save File', default_path, 'csv(*.csv);;all (*)')
            elif 'PyQt5' in sys.modules:
                save_file_name, _ = QFileDialog.getSaveFileName(self, 'Save File', default_path, 'csv(*.csv);;all (*)')
            
            # don't update the threshold  - trust the user to have already set it
            # if not self.thresh_toggle.isChecked(): # update the threshold unless it's set manually
            #     self.plot_current_hist(self.image_handler.hist_and_thresh)
            self.add_stats_to_plot()

            # include most recent histogram stats as the top two lines of the header
            self.image_handler.save_state(save_file_name,
                         hist_header=list(self.histo_handler.temp_vals.keys()),
                         hist_stats=list(self.histo_handler.temp_vals.values())) # save histogram
            
            try: 
                hist_num = self.histo_handler.stats_dict['Hist ID'][-1]
            except IndexError: # if there are no values in the stats_dict yet
                hist_num = -1

            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("File saved to "+save_file_name+"\n"+
                    "and appended histogram %s to log file."%hist_num)
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

        except OSError:
            pass # user cancelled - file not found

    def check_reset(self):
        """Ask the user if they would like to reset the current data stored"""
        reply = QMessageBox.question(self, 'Confirm Data Replacement',
            "Do you want to discard the current data?", 
            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)
        
        if reply == QMessageBox.Cancel:
            return 0

        elif reply == QMessageBox.Yes:
            self.image_handler.reset_arrays() # gets rid of old data

        return 1

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
            self.image_handler.reset_arrays() # get rid of old data
            self.hist_canvas.clear() # remove old histogram from display

        elif reply == QMessageBox.No:
            self.image_handler.reset_arrays() # get rid of old data
            self.hist_canvas.clear() # remove old histogram from display
        

    def load_from_files(self, trigger=None):
        """Prompt the user to select image files to process, then sequentially process
        them and update the histogram"""
        default_path = self.get_default_path(option='im')
            
        if self.check_reset():
            try:
                self.recent_label.setText('Processing files...') # comes first otherwise not executed
                if 'PyQt4' in sys.modules:
                    file_list = QFileDialog.getOpenFileNames(self, 'Select Files', default_path, 'Images(*.asc);;all (*)')
                elif 'PyQt5' in sys.modules:
                    file_list, _ = QFileDialog.getOpenFileNames(self, 'Select Files', default_path, 'Images(*.asc);;all (*)')
                for file_name in file_list:
                    try:
                        self.image_handler.process(file_name)
                        self.recent_label.setText('Just processed: '+os.path.basename(file_name)) # only updates at end of loop
                    except:
                        print("\n WARNING: failed to load "+file_name) # probably file size was wrong
            
                self.update_stats()
                if self.recent_label.text == 'Processing files...':
                    self.recent_label.setText('Finished Processing')

            except OSError:
                pass # user cancelled - file not found

    def load_from_csv(self, trigger=None):
        """Prompt the user to select a csv file to load histogram data from.
        It must have the specific layout that the image_handler saves in."""
        default_path = self.get_default_path()
            
        if self.check_reset():
            try:
                # the implementation of QFileDialog changed...
                if 'PyQt4' in sys.modules: 
                    file_name = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'csv(*.csv);;all (*)')
                elif 'PyQt5' in sys.modules:
                    file_name, _ = QFileDialog.getOpenFileName(self, 'Select A File', default_path, 'csv(*.csv);;all (*)')
                self.image_handler.load_from_csv(file_name)
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
            success = self.histo_handler.load_from_log(file_name)
            if not success:
                print('Data was not loaded from the log file.')
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
            print("\nFile processing event duration: %.4g "%(self.dir_watcher.event_handler.event_t*scale)+unit)
            print("Most recent idle time between events: %.4g "%(self.dir_watcher.event_handler.idle_t*scale)+unit)
            print("Image processing duration: %.4g "%(self.int_time*scale)+unit)
            print("Image plotting duration: %.4g "%(self.plot_time*scale)+unit)
            print("File writing duration: %.4g "%(self.dir_watcher.event_handler.write_t*scale)+unit)
            print("File copying duration: %.4g "%(self.dir_watcher.event_handler.copy_t*scale)+unit)
        else: 
            print("Initiate the directory watcher before testing timings")

    #### #### UI management functions #### #### 

    def closeEvent(self, event):
        """Prompt user to save data on closing"""
        reply = QMessageBox.question(self, 'Confirm Action',
            "Save before closing?", QMessageBox.Yes |
            QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)

        if reply == QMessageBox.Yes:
            self.save_hist_data()         # save current state
            
            if self.dir_watcher:          # make sure that the directory watcher stops
                self.dir_watcher.observer.stop()
                
            event.accept()
        
        elif reply == QMessageBox.No:
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
    main_win.show()
    if standalone: # if an app instance was made, execute it
        sys.exit(app.exec_()) # when the window is closed, the python code also stops

            
if __name__ == "__main__":
    run()