# DetEdit

Visual user-interface for acoustic detection annotation designed for use with HARP data (x.wav files).

----------
Author: John A. Hildebrand, based on initial version by Sean M. Wiggins.

Copyright: J. A. Hildebrand 2016

Date: 1/26/2016


### Workflow

#### 1) Detection

- `Edetect.m`

	`Edetect` is a basic energy detector designed for use with xwav files.

	#####Input:

	User is prompted to supply inputs including:

	a) A text file or spreadsheet containing detection parameters (see `GOM_BW_pfile_320.xlsx` for example).
		
	b) A directory containing xwavs.
	

	#####Output:

	A `*_TPWS.mat` file containing matrices of detected signal parameters.

		- MPP = peak to peak amplitudes 		

		- MTT = time 

		- MSN = timeseries (bandpassed)

		- MSP = spectrum (bandpassed)

		_Optional_:

		- MUSN = unfiltered timeseries

		- MUSP = unfiltered spectrum


#### 2) Make LTSA snippets

- `mkLTSAsessions.m`

	`mkLTSAsessions` prepares small LTSAs for each detection bout.

	##### Input:

	a) User is prompted to input a code to identify the detected species. 

	Current options include:

	`Ko` or `k`: Kogia

	`Zc` or `z`: Cuvier's beaked whale

	`Me` or `m`: Gervais' beaked whale

	`BWG` or `g`: Unidentified beaked whale BWG

	`Md` or `d`: Beaked whale BW31

	`De` or `de`: Delphinid

	b) User is prompted to select transfer function file
	

	##### Output:

	A `*_LTSA.mat` file is produced.


#### 3) Edit detections

- `detEdit.m`

	##### Input:

	User is prompted to:

	a) input a code to identify the detected species (as in step 2).

	b) input an interval for looking at false detections (To estimate a false positive rate, a good number might = total # of detections/N, where N is 300 or more.) 

	c) select a tranfer function

	d) select directory containing `_TPWS.mat` and `_LTSA.mat` files.

	e) Starting session (use `1` to start with first bout)

	##### Editing tools & codes:

	TODO

	##### Output

	TODO

#### 4) Update detections

- `modDet.m`

	##### Input

	TODO

	##### Output

	TODO
