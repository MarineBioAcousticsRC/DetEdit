# DetEdit

Visual user-interface for acoustic detection annotation designed for use with HARP data (x.wav files).

----------
Author: John A. Hildebrand, based on initial version by Sean M. Wiggins.

Copyright: J. A. Hildebrand 2016  
Date: 1/26/2016


### Workflow

#### 1. Detect

##### Edetect.m  
  `Edetect` is a basic energy detector designed for use with xwav files.
  
  - **Input:**  
    User is prompted to supply inputs including:  
     1. A text file or spreadsheet containing detection parameters (see `GOM_BW_pfile_320.xlsx` for example).  
     2. A directory containing xwavs.
    
  - **Output:**  
   A `*_TPWS.mat` file containing matrices of detected signal parameters.

     MPP = peak to peak amplitudes  
     MTT = time  
     MSN = timeseries (bandpassed)  
     MSP = spectrum (bandpassed)  

  *Optional*:  
     MUSN = unfiltered timeseries  
     MUSP = unfiltered spectrum  

#### 2. Make LTSA snippets

##### mkLTSAsessions.m  
  `mkLTSAsessions` prepares small LTSAs for each detection bout.  

  - **Input:**  
   1. User is prompted to input a code to identify the detected species.  
     Current options include:  
     `Ko` or `k`: Kogia  
     `Zc` or `z`: Cuvier's beaked whale  
     `Me` or `m`: Gervais' beaked whale  
     `BWG` or `g`: Unidentified beaked whale BWG  
     `Md` or `d`: Beaked whale BW31  
     `De` or `de`: Delphinid  
     `Pm` or `pm`: Sperm whale

   2. User is prompted to select transfer function file  

  - **Output:**  
   A `*_LTSA.mat` file is produced.

#### 3. Edit detections

##### detEdit.m

   - **Input:**  
   User is prompted to:  
     * input a code to identify the detected species (as in step 2).  
     * input an interval for looking at false detections.  
     (To estimate a false positive rate, a good number might = total # of detections/N, where N is 300 or more.)  

     * select a tranfer function.   
     * select directory containing `_TPWS.mat` and `_LTSA.mat` files.  
     * Starting session (use **1** to start with first bout).  

   - **Editing tools & shortcuts:**  
  	Brushing - Matlab's paintbrush tool is used to label detections.  
  	Colors have different meanings:  
        	Red:  False positive  
		Black:  True detection  
		Yellow: Temporarily shows details of brushed clicks in black outline but does not change their designation 
		Bright Green:  Misidentified detection  
		10 other colors available for ID labeling include:  

			[255, 153, 200] = type 1 pink  
    			[218, 179, 255] = type 2 purple  
    			179,  200, 255] = type 3 light-blue  
    			[174, 235, 255] = type 4 pale-blue  
    			[0,   255, 255] = type 5 cyan  
    			[255, 177, 100] = type 6 peach  
    			[255,   0, 255] = type 7 magenta  
    			[122,  15, 227] = type 8 purple  
    			[20,   43, 140] = type 9 dark blue  
    			[221, 125,   0] = type 10 orange  


   Keyboard Shortcuts:  
	'r' Label currently selected clicks as false  
	'f' Label all clicks in current window as false  
	'i' Label currently selected clicks as true (does not work if brush color is red,yellow, or bright green)  
	't' Label all clicks in current window as true  
	'm' Label all clicks in current window as misidentified  
	'y' Display summary info for only currently selected clicks (mean spectrum and mean timeseries will be shown in black)  
	'u' Update window contents according to current brush selection and color  
	'j' Jump to a non-consecutive session (prompt in matlab command window will ask you for a session number  
	'b' Go back one bout  
	'a' Adjust LTSA contrast and brightness  
	's' Update maximum ICI scale  
	'<' Change RMS threshold in plot 51  
	':' Change PP threshold in plot 51 
	'^' Change high frequency threshold in plot 51  
	'!' Change frequency scale in plot 53
	'^' Change high frequency threshold in plot 53  
     	'd' Change recieved level scale on top subplot of fig 201  
    	'x' or 'z' Test a random subset of detections to estimate false positive cue rate  
    	'w' Test a random subset of time bins to estimate false positive bin rate  

   - **Output:**  
     `*ID.mat` - This file contains a 2xn matrix in which the first column contains the detection time, and the second column contains an ID number associated with that detection. The ID number is given by the brush color used by the analyste to identify that detection type.  

     `*FD.mat` - This file contains a 1Xn vector consisting of the times of all detections marked as false positives by the user.  

     `*MD.mat` - This file contains a 1Xn vector consisting of the times of all detections marked as misidentifications by the user.  

     NOTE: `ID.mat`,`FD.mat`, and `MD.mat` files are updated each time a bout is edited.

#### 4. Update detections

##### modDet.m

   - **Input:**
  ToDo

   - **Output:**
  ToDo
