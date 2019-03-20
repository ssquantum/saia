"""Single Atom Image Analysis
Stefan Spence 14/03/19

Separate out the image_handler class for processing single atom images from the
director watcher and Qt GUI. This allows it to be imported for other purposes.

Assume that there are two peaks in the histogram which are separated by a 
region of zeros.

Assuming that image files are ASCII.
"""
import os
import sys
import numpy as np


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
                              usecols=range(1,self.pic_size+1))
        
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
            
