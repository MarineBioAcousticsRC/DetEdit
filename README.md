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
     Red:
     Green:
     
    

   - **Output:**
     *ID.mat - This file contains a 2xn matrix in which the first column contains the detection time, and the second column contains an ID number associated with that detection. The ID number is given by the brush color used by the analyste to identify that detection type.
     *FD.mat - This file contains a 1Xn vector consisting of the times of all detections marked as false positives by the user.
     NOTE: ID.mat and FD.mat files are updated each time a bout is edited.

#### 4. Update detections

##### modDet.m

   - **Input:**
  ToDo

   - **Output:**
  ToDo
