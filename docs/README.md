**** Version 1.1 ****
Produces histograms for an image containing a either single atom or none.
**** ****

How to run Single Atom Image Analysis (SAIA):
	• Start the file: 
		○ Execute run_with_enthought.bat   --- a windows batch file with a hardcoded link to the Enthought python executable
		○ Execute run_with_conda.bat       --- activate the Anaconda environment (you must first create the saiaenvironment, which can be done using create_environment.bat) and run using Anaconda.
		○ Or run from a python distribution (e.g.  python main.py)
		
	• A window pops up showing the loaded file config and asking to start the directory watcher
		○ 'Yes' will start the directory watcher to process file creation events from a directory.
		○ 'No' starts up the program without the directory watcher (it can be initiated later)
Image storage path	Where SAIA will save images to (in subdirectories by date)
Log file path		Where SAIA will save log files to (in subdirectories by date, collects histogram statistics)
Dexter sync file	Absolute path to the file where Dexter stores the current file number
Image read path		Absolute path to the folder where Andor will save new image files to (note that no other file creation events should occur in this folder, or they will be processed by the directory watcher as well)
Results path		The default location to open the file browser for saving csv files
	
	• Note that the image size in pixels must be set before any images are processed. 
		○ If the image size is known, type it into the 'Image size in pixels:' text edit
		○ The image size can also be taken from an image file by clicking 'Load size from image'
		○ The 'Get ROI from image' will set the image size and then centre the ROI on the pixel with maximum intensity
		
	• There are several running modes:
		○ Active directory watcher (real time processing of images straight after the file is saved to the image read path. Copies then deletes images)
		○ Passive directory watcher (real time processing of images straight after the file is saved to the image read path. Doesn't alter the file)
		○ Load data from csv (the format is: file#, counts, atom detected?, max count, pixel x position, pixel y position, mean count, standard deviation)
		○ Load data from a selection of image files
		
	• Note that when loading in new data it will use the current ROI settings on display. It will ask whether you want to clear the current array, which will prevent mixing of data with different ROI settings.
	
	• For the current histogram there are several binning options:
Automatic Binning (Default)	Numpy decides the binning automatically
Manual	When the Max, Min, and #Bins text edits are populated, the histogram will be set with those limits
No Update	The directory watcher will still run and files will still be processed, but the histogram will not be replotted (speeds up processing)
	
	• Selecting 'Auto-Display Last Image' plots a 2D colourmap of the image file last processed.
		○ This can take up to 1s for 512x512 images and so causes lag if file events occur faster than this.
		○ The user can set an ROI by clicking 'ROI' and then dragging the box:
			§ Dragging from the box area translates the box
			§ The top left circle can be used to rotate the box
			§ The bottom right square can be used to resize the box
		○ The ROI can also be set by the text inputs in the settings tab (all must be filled in before there is any change)
		○ The ROI can be centred on the max pixel in an image by clicking 'Get ROI from image'
		○ Changing the ROI by either of these ways sets the region of the image that is processed, and this will be retained until the settings are next changed. 

We will take fluorescence images of atoms in the tweezer from the Andor camera (512x512, well-depth 180,000 e-)
We want a program that will real-time readout the integrated counts from the image and display a graph which identifies images with or without single atoms in.

Results structure
	• An image is processed by saving a copy to the image storage path then calculating:
File #				Taken from currentfile.txt
Integrated counts in ROI	User sets ROI, sum the counts in all of the pixels
Atom detected			Counts // threshold. This is greater than zero if an atom is detected
Mid count			the count in the pixel at the centre of the ROI
xc				x-position of max count (in pixels)
yc				y-position of max count (in pixels)
Mean count			Take the mean of the image outside of the ROI to estimate background
Standard deviation		Take the standard deviation of the image outside of the ROI
The integrated count is added to the current histogram and the rest of the information is stored in an array
These are also the column headings when the histogram is saved.
	
	• The histogram statistics are analysed and displayed in the 'Histogram Statistics', they will be appended to a log file when the histogram csv is saved. The log file contains the following columns:
Histogram # 			Increments by one every time a line is appended to the log file
Variable			Variable set by the user, must be a float
Loading probability		Ratio of images with counts above threshold to total images processed
Error in loading probability	The statistical confidence in the loading probability from the binomial distribution
Background peak count		Calculated position of the background peak in counts
Background peak width		Calculated width of the background peak in counts
Signal peak count		Calculated position of the atom peak in counts
Signal peak width		Calculated width of the atom peak in counts
Fidelity			The probability of correctly identifying atom presence
Error in fidelity		The range of possible fidelities based on peak position uncertainty
Separation			Signal peak count - background peak count
Threshold			Mean of background and signal peak counts
Images processed	Number of images processed in the current histogram

	• Each time a line is appended to the log file, the data will also be added to arrays which can be accessed in the 'Plotting' tab.
This happens when a histogram is saved, or when the 'Add to plot' button is pressed. 
	
Peak calculations
There are several different ways to estimate the background and signal peak centres and widths and therefore calculate the threshold:
	1) Scipy.signal.find_peaks - quite good at finding peaks but it isn't clear what the width is
	2) Use a threshold to split the sorted list of counts from each image. Take the mean and standard deviation of the upper (atom) and lower (background) arrays then set a new threshold at 5σ above background (not yet implemented)
	3) Use a threshold to split the histogram into background and single atom peaks, then fit Gaussian curves to get the mean and standard deviation.


	• Any of the plots can be saved by right clicking on the plot area and selecting 'Export…'
Under 'Item to export' make sure to select 'Entire Scene', otherwise it will not save.


Settings Tab

Histogram Tab

Histogram Statistics Tab

Image Tab

Plotting Tab

-----------------------------------------------------------------------------------
Possible Architectures:

Original method:
	Experiment run by Dexter in loop -> Andor running in loop receives trigger (saves images) -> Python runs a loop checking for file changes -> python processes the new file, then plots a histogram
	Needs several threads to run in parallel – a directory watcher to notice created files, and a graph plotter to process the data. Working version using pyqt

Andor Trigger Method 1:
	Experiment run by Dexter in loop -> Andor running in loop receives trigger (saves images) -> triggers Python with TTL to process new image file
	Maybe Python only updates the histogram every [x] number of runs? Python still needs to run continuously waiting for a trigger from Andor if we want it to retain data, but is told when files are saved rather than having to watch and wait.
	Python could send a trigger to start the next Dexter run
	
Andor Trigger Method 2:
	Experiment run by Dexter in a loop -> Andor running in loop receives trigger (saves images with Dexter file #) 
	Python script runs independently analysing files. Since Andor would save separate files with their Dexter label, then if python lags behind it will not lose sync. This requires Andor staying in sync with the Dexter file #
	
Dexter Saves Method:
	Write a .vi in LabView that gets Dexter to control the camera instead of Andor SDK. Then Dexter can save the files with its Dexter file # and we have no problem with syncing. Analysing the data can then be done separately with no time pressure. The Strontium project already do this with their MPD SPC3 SPAD.
	This is the ideal method that will be implemented at some point.

Input:
Andor saves an image to a set directory as it's running. 
Background subtraction is probably not necessary since it would just shift the histogram along the x axis.
The Andor camera will be set to take an image of a Region of Interest (ROI) around the atom.
The Andor camera can also bin pixels to combine their counts (eg. 4x4 pixel array -> 1 binned pixel) but we might do this binning in post-processing.

Original Method:
1) scan Andor output directory for presence of a new file (dir watcher)
2) save the image with a new label (dir watcher emits signal to dir watcher: on_created)
	[species]_[date]_[dexter file #]
	^requires synchronisation with dexter for file #. Might be slow over the network, or if the network drops then it goes out of sync.
	For multiruns we can get the file numbers from Dexter's outputted measure file
	Then delete the old image file so that a new one can be saved with the same name.
3) load image (dir watcher emits signal to image handler)
	• Takes 10 - 50 ms but will miss files and crash with an exception is the rate is too fast (works for 0.5s delays but crashes for 0 delay)
4) Find the position of the atom. Not really necessary if we've already selected an ROI.
	• Fit a Gaussian (should only be over a couple of pixels, fitting makes it slower...)
	• Locate the brightest pixel
5) Integrate the atom signal into a count
6) determine whether a single atom is present through reference to a threshold count
	Need some way to estimate parameters of peaks to separate them 
	The separation between peaks (1 atom, 2 atoms, etc) should be the same
	Could estimate threshold as:   
		○ midway between peaks 
		○ The middle of the gap between peaks
		○ Where there is maximum curvature (not robust)
		○ A set number of standard deviations above background (this method is used since it fixes our statistical confidence in the atom signal)
7) plot the integrated count in a histogram (image handler emits signal to histogram plotter)
	Preferably real-time but if that's too slow then every 100 shots or so.
	Would also be good to display the image but again might be too slow for each shot.
	Takes 10 - 50 ms (including loading image)
8) use the histogram to update threshold values (contained in the image handler)
	might need upper bound for two atoms trapped as well as lower bound for no atoms trapped
9) save the histogram with references to the image files so that the data can be re-analysed later as well
	Preferably also storing the input parameters, where available (dexter measure file saves input parameters for a multirun)
	10) Collect histogram statistics in a log file and add them to a plot each time a histogram is saved

Output:
	• A directory with subdirectories ordered by date storing all of the labelled image files
	• A real-time display of the histogram 
	• A summary (log) file with histogram data:
		○ Essential: files included, so that they can be re-analysed at a later date. Columns are:
			§ File label (dexter #)
			§ Integrated counts
			§ Variable (if possible - taken from dexter multirun)


Side tasks and notes:
	• Test the camera computing speed – how fast can it output files? It's faster with just an ROI
	• Will need to test timings: how fast can you process an image? How fast can you plot? How fast does the experiment run? How fast does the andor camera output the file? How fast are files saved? How fast can they be accessed over the network?
	• Long term stability - will it go out of sync if the network goes down?
	• Simulate acquiring an image using the red guide beam
		○ randomly fire/don't fire the AOM (requires an analog channel) – then we can test the success rate of identifying an 'atom'
		○ Set the exposure to the anticipated experimental value (20ms) and include MOT beams etc. To simulate experimental background.
	• Different ways to trigger; in order to time the imaging we need a trigger sent from the andor software (TTL)
		○ We can connect a TTL to python from Andor and from python to dexter via USB.
		○ After Andor camera script runs -> TTL to python 
		○ After python script runs -> TTL to dexter
		○ Then the next dexter experimental run could be started without losing sync (but has to wait)
	• We could measure one of the TTL channels on an oscilloscope, get the time between its triggers and subtract off the duration of the python script and the experiment duration. This would give us the time from camera trigger to python script starting (which is taking the images then saving the files).
		○ However this might not include the time that windows processes files after Andor tells them to be saved
		○ Further, subtracting similar numbers increases the relative numerical error
	• we could send an auxout TTL from the camera after the images have been taken to give readout speed. Then after the files have been saved to test that duration.


Timing for the image analysis program:
	• Time to update plot (2s between image creation): 3 – 4 ms
	• Time to make histogram (2s between image creation): 2 – 230 ms
	• Time to copy file (2s between image creation): 8 – 12 ms (sometimes tries to access an empty file because the copying isn't finished)
Updated timings:
	• Plotting: 10-15ms
	• Make histogram: 
		○ 250-350ms while live plotting (512x512) (probably because of lag from interface)
		○ 10ms while live plotting (64x64)
	• File copying event: 
		○ Occasionally seen to take 200-300ms when there is lag from live plotting
		○ 1-4ms while live plotting (64x64) and (512x512)
	

Timing for the camera:
	• Time to take images (readout - change with ROI?): 512x512 – 800ms, 64x64 - 300ms
	• Time to save images: 4ms


Debugging:
	• Care had to be taken to ensure that the directory watcher was stopped when the program was closed, otherwise its thread would keep running in the background which could cause overwriting of files.
	• It was found that when the directory watcher triggered on file modified events, it would recognise several events for the same image being saved (i.e. for one experimental run it would try and save the image up to 8 times). Particularly problematic was that the first of these would be before the Dexter sync file had been updated, and so it would overwrite the file with the previous Dexter number as well as writing to the current one.
		○ This was solved by making adding some delays and making the directory watcher trigger on file creation events only, then delete the file after a copy has been saved
	• Sometimes when the experiment is running fast the dir watcher processes the file event before Dexter has changed the current file number. This causes python to go out of sync and could possibly overwrite the previous file.

Future Developments:
	• Set several ROIs and make one of a set of TTLs high if an atom is detected in a particular ROI. Could be used as a trigger for the next Dexter run
	• Decluttering the display of a single image (removing histogram at the side)
	• Fix the intensity and zoom of an image for comparison when new images come in

Python script can read in/out TTL through USB
Use to check timings like an oscilloscope
Dexter can wait for a TTL signal from python