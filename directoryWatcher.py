"""Single Atom Image Analysis
Stefan Spence 12/04/19

 - watch the image_read_path directory for new images
 - save the new image with label into a dated subdirectory under image_storage_path
 - delete the original file so that a new file with the same name can be created
 
Assuming that image files are ASCII

watchdog creates an observer that waits for file creation events
the observer must be initiated and shut down properly to ensure that there isn't
one running behind the scenes which might overwrite previously saved files.
"""
import numpy as np
import os
import time
import shutil
try:
    from PyQt4.QtCore import QThread, pyqtSignal, QEvent
except ModuleNotFoundError:
    from PyQt5.QtCore import QThread, pyqtSignal, QEvent
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

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
    
    def on_created(self, event):
        """On a new image being written, save the file with a synced label into the 
        image storage dir"""
        t0 = time.time()
        self.idle_t = self.end_t - t0   # duration between end of last event and start of current event
        # if event.event_type == 'created' or event.event_type == 'modified': # ignoring deletion or move events
        # print("--- --- --- ---")
        # print(event.event_type)
        # print(event.src_path)
        
        time.sleep(0.3) # need to wait sufficient time for Dexter to update the file number
        # get Dexter file number  
        with open(self.dexter_sync_file_name, 'r') as sync_file:
                self.dfn = str(int(sync_file.read()))
                
        # copy file with labeling: [species]_[date]_[Dexter file #] ---- this will overwrite if file already exists
        new_file_name = self.image_storage_path+r'\Cs-133_'+self.date+'_'+self.dfn+'.'+event.src_path.split(".")[-1]
        try:
            shutil.copyfile(event.src_path, new_file_name)
        except PermissionError:
            print("WARNING: added a pause because python tried to access the file before the other program had let go")
            time.sleep(0.2)
            shutil.copyfile(event.src_path, new_file_name)
        time.sleep(0.1)  # add in latency in case it tries to delete the file before it's been saved.
        
        try:
            os.remove(event.src_path)  # delete the old file so that we can see a new created file event
        except PermissionError:
            print("WARNING: added a pause because python tried to access the file before the other program had let go")
            time.sleep(1)
            os.remove(event.src_path)
        
        # time.sleep(0.5)
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
        path_label_text = ['Image Storage Path: ', 'Log File Path: ', 
                'Dexter Sync File: ', 'Image Read Path: ', 'Results Path: ']
        self.dirs_dict = {key:value for (key, value) in 
                        np.array([path_label_text, self.get_dirs()]).T}  # handy list contains them all
        (self.image_storage_path, self.log_file_path, self.dexter_sync_file_name, 
                    self.image_read_path, self.results_path) = self.get_dirs()
        
        # create the watchdog object
        self.observer = Observer()
        
        # get the date to be used for file labeling
        self.date = time.strftime("%d %b %B %Y", time.localtime()).split(" ") # day short_month long_month year
        self.image_storage_path += r'\%s\%s\%s'%(self.date[3],self.date[2],self.date[0])
        
        self.event_handler = system_event_handler(self.image_storage_path, 
                                self.dexter_sync_file_name, self.date[0]+self.date[1]+self.date[3])
        
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
            elif "results path" in row:
                results_path = row.split('--')[-1]       # default folder to save histogram csv files to
                
        if os.path.split(dexter_sync_file_name)[0] == image_read_path:
            print("WARNING: The directory watcher acts on all file change events, so the Dexter sync file path and image read path must be different.")
        return [image_storage_path, log_file_path, dexter_sync_file_name, image_read_path, results_path]
        
    @staticmethod
    def print_dirs(image_storage_path, log_file_path, dexter_sync_file_name, image_read_path, results_path):
        """Return a string containing information on the paths used"""
        outstr = '// list of required directories for SAIA\n'
        outstr += 'image storage path\t--'+image_storage_path+'\n'
        outstr += 'log file path\t\t--'+log_file_path+'\n'
        outstr += 'dexter sync file\t\t--'+dexter_sync_file_name+'\n'
        outstr += 'image read path\t--'+image_read_path+'\n'
        outstr += 'results path\t\t--'+results_path+'\n'
        return outstr
    
    def run(self):
        pass
        
    def save_config(self):
        with open('./config.dat', 'w+') as config_file:
            config_file.write(self.print_dirs(*[d for d in self.dirs_dict.values()]))
  