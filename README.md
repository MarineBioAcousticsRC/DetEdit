# DetEdit
A graphical user interface for annotating and editing events detected in long-term acoustic monitoring data

----------
## Introduction
*DetEdit* is a tool to visualize and annotate events from large acoustic datasets. It is highly-configurable and can be used for multipurposes ranging from annotation and classification of individual or batches of signals and evaluation of signal properties, to removing false detections and obtaining false positive rates.

The tool has been used mainly to process extensive acoustic datasets of odontocete echolocation clicks and human impulsive noise. It can be used to process any stereotyped impulsive signals such as species-specific acoustic signals or calls from fish, crustaceans, bats, birds, insects, etc.


## Software and data requirments
*DetEdit* is a MATLAB-based graphical user interface (GUI). It requires: 
* MATLAB R2014b or newer versions (www.mathworks.com). A version of the repository for older versions of MATLAB (R2013b and older) is also provided.
* Audio files in WAV or XWAV format

## Runtime and setup

### Obtaining DetEdit
- Download repository at https://github.com/ScrippsWhaleAcoustics/DetEdit
- Examples of DetEdit files with different odontocete species: https://drive.google.com/drive/folders/0B1N3RJM5Uw4ha25lVEhLRF9FMUk

Set up *DetEdit* repository folder in MATLAB's path by using the MATLAB File pull-down and *Set Path...*. Click the *Add with Subfolders...* button, browse and select the folder containing the version of *DetEdit* that you want to run. Remove other versions of *DetEdit* from the path with the *Remove* button. Click *Save* and then *Close*.


## Typical Workflow

### Initial data preparation
Create parameter files required to use the interface:

#### 1. Make Long-Term Spectral Average (LTSA) files
LTSA files (`.ltsa`) are generated from a collection of WAV or XWAV files by averaging spectra over a long time periods and arranging these spectra sequentially as frequency-time spectrogram plots. LTSA files provide a quick overview of long-term recordings.

Using the Command Window in MATLAB, the `.ltsa` file is created as follows:   
```bash
> mkLtsa
```
and follow the prompt windows to specify audio file format, file directories and parameters (e.g. time average length and frequency bin size to average the data).

Only 5 filename formats of WAV and XWAV files are supported:
- `yymmdd-HHMMSS`: SiteName_190326-123700.wav
- `yymmdd_HHMMSS`: SiteName_190326_123700.wav
- `yyyymmdd_HHMMSS`: SiteName_20190326_123700.wav
- `yymmddHHMMSS`: SiteName_190326123700.wav
- `yyyymmddTHHMMS`: SiteName_20190326T123700.wav

#### 2. Make `TPWS` files (start Time, Peak-to-peak amplitude, Waveform and Spectra parameters) 
A `*_TPWS.mat` file contains the following matrices of detected signal parameters:

|#| Variable | Description                               |
|-|----------|-------------------------------------------|
|1|   `MTT`  | Vector of start times of detections       |
|2|   `MPP`  | Vector of received level amplitudes (dB<sub>pp</sub>)|
|3|   `MSP`  | Matrix of detection spectra               |
|4|   `MSN`  | Matrix of waveforms                       |

Provide `*_TPWS.mat` file with these matrices is required:
- user creates matrices manually
- if user has start times of detected acoustic signals, can create the `*_TPWS.mat` file as follows:

  ```bash
  > edit make_TPWS
  ```
  In editor window, modify input/output locations and detection parameters to run script.
- if user has no detections, a generic detector can be applied,
	- for use with XWAV files:
	  ```bash
          > Edetect
  	  ```
	  User is prompted to supply inputs including a text file or spreadsheet containing detection parameters (see `GOM_BW_pfile_320.xlsx` for example), and directory containing XWAV files.
		
	- for use with WAV files:
	  ```bash
          > Edetect_wav ('paramFile','E:\DetEdit\GOM_BW_pfile_320.xlsx',...
          'tfFile','E:\MyTransferFunctionFiles\example_sig1_invSensit.tf',...
          'timeFile','E:\DetEdit\BW_Effort.xlsx',...
          'wavDir','E:\MyWAVFiles',...
          'channel',1)
  	  ```

#### 2. Make LTSA snippets
A `*_LTSA.mat` file contains the following matrices of detected signal parameters:

|#| Variable | Description                               |
|-|----------|-------------------------------------------|
|1|   `pt`   | Vector of start times of spectral averages|
|2|   `pwr`  | Matrix of power spectral densities        |


+++++++++++++++ continue editing+++++++++++++++

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

	'p' Activate brush (need it for brushing in MATLAB2014b or higher). It requires to click the figure that you want to brush before pressing the key.

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
