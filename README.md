Radboud Memory Dynamics - Calcium Imaging Analysis
=====

This is the main calcium imaging analysis tool for the Radboud Memory Dynamics group.

The project structure is based on Cookie cutter Data Science by Driven Data. The analysis tools are provided by CaImAn by Flatiron Institute. 

This project is a joint work done by Melisa Maidana Capitan (m.maidanacapitan@donders.ru.nl),
Casper ten Dam, Sebastian Tiesmeyer, Morgane Audrain under the supervision of Francesco Battaglia. 

Project Organization
 -----
    
    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    
General 
-----

The calcium imaging program have two ways of working, you can choose to run:

    •	Pipeline

    •	Pipeline_parameters_settings

The pipeline code run the all pipeline with define parameters that you can find on the SQL database,
I will come back to the workability of this pipeline later.

The pipeline parameters settings is a pipeline where you can choose to run one or two trial of a specific mouse 
or a specific session for set up the more accurate parameters for this session or entire mouse.

The parameters settings pipeline plot for your different picture of each for you to choose the cropping region 
or the filter that you want to use for example.
You can run with this pipeline every step that you want but is not necessary and useful to run all the steps in row with this pipeline.
Once the parameters choose and upgrade on the database you can return to use the global pipeline.
Now that the difference between the two pipelines are explained, I will go more in details about the general use of the pipeline.

For a normal use of the pipeline you are not obligate to go through the organization of each steps, you can just run this pipeline. 

When you start the pipeline, you will be asked to complete some settings:
    
    •	the mouse number:  Identification number of the mouse that you want to analyze    
                → you can analyze just one mouse by one mouse, so you need to enter here just one number. For example, 32364
            
    •	Sessions: Sessions that you want to analyze for this specific mouse, you can choose to analyze one session or more. For multiple sessions you need to enter an array. For example, one session 1, multiple session [1 2 3]
    
    •	Begin of trial: the first trial where you want to start the analysis
    
    •	Last trial: Last trial you want to analyze

After setting those settings the pipeline will ask you which steps you want to run. The choice that you have it between:

    •	Decoding
    
    •	Cropping
    
    •	Motion correction
    
    •	Alignment
    
    •	Equalization
    
    •	Source extraction
    
    •	Component evaluation
    
    •	Registration
    
    •	All

 I will come back later the use and option of this different steps. After running the steps that you select the final question of the pipeline is if you want to run other steps with the same parameters meaning the same mouse, sessions and trial. If you decide to don’t run another step the pipeline stop every process that can still run and display the goodbye message.


---
1. Decoding
---

Decodes the raw files with Inscopix.

This function only needs the mouse number, the session, and trials that you want to decode. Create a tif file. 

---
2. Cropping
---

Crop the region of interest (ROI) for later analysis. 

This step takes a decoded file and crop it according to specified cropping points.

This input is a tif file and the output is a tif file as well.
You do not need to specify here the decoding version because is by default version one. 
You will be ask here which cropping region you want to applicate, for visualization purpose you can choose to run first pipeline_parameters_settings, choose the parameters that you want and after applicate those parameters for the all session.
The spatial cropping points are referred to the selected interval. y1: start point in y axes, y2: end point in y axes. z1: start point in z axes, z2: end point in z axes. 
TEST CROPPING: parameter_setting scripts can be run in a small section of the filed of view (FOV) in order to be able to run faster multiple motion correction and source extraction parameters. But, it needs to be verified that the convergence of an extended area of the FOV will return ~ the same neural footprints. For that, may use of script in parameters_setting_cropping_impact. 
In this script we select 5 different cropping with different sizes, and run the pipeline (with out component evaluation). That output is a figure that is saved in  directory = '/home/sebastian/Documents/Melisa/calcium_imaging_analysis/data/interim/cropping/meta/figures/cropping_initialization/'. This figure contains tha contour plots for different cropping, but always using the same motion correction and source extraction parameters. This step is important because the algorithm is really sensitive to changes in the seed or initial condition. Changes in the size of the FOV of course will lead to changes in the initial condition, that can lead to huge differences in the extraction. Do not forget to verify before generalizing your parameters to a bigger FOV. 



3. Motion correction
Motion correction can be run both in the server and it the local machine, but it is supposed to be faster in the server. 

Motion correction parameters as specified in a dictionary as follows : 

	parameters_motion_correction = {'motion_correct': True, 'pw_rigid': True, 'save_movie_rig': False,
              'gSig_filt': (7, 7), 'max_shifts': (25, 25), 'niter_rig': 1, 'strides': (96, 96),
              'overlaps': (48, 48), 'upsample_factor_grid': 2, 'num_frames_split': 80, 'max_deviation_rigid': 15,
              'shifts_opencv': True, 'use_cuda': False, 'nonneg_movie': True, 'border_nan': 'copy'}

	# motion correction parameters
	motion_correct = True    # flag for performing motion correction
	pw_rigid = True          # flag for performing piecewise-rigid motion correction (otherwise just rigid)
	gSig_filt = (3, 3)       # size of high pass spatial filtering, used in 1p data
	max_shifts = (5, 5)      # maximum allowed rigid shift
	strides = (48, 48)       # start a new patch for pw-rigid motion correction every x pixels
	overlaps = (24, 24)      # overlap between pathes (size of patch strides+overlaps)
	max_deviation_rigid = 3  # maximum deviation allowed for patch with respect to rigid shifts
	border_nan = 'copy'      # replicate values along the boundaries


This script runs a few examples of gSig_filt size and saves the figure in folder: 

'/home/sebastian/Documents/Melisa/calcium_imaging_analysis/data/interim/motion_correction/meta/figures'

for each mouse,session,trial,is_rest parameters and version analysis. 

Once several parameters of motion_correction had been run, it is necessary to compare performance. Easier quality measurement of motion correction is crispness. And the end the script goes over all the analysed motion corrected versions (using function compare_crispness that is defined in analysis.metrics), and with plot_crispness_for_parameters (analysis.figures) creates a plot for cripness for the mean summary image and for the corr summary image. 


Also the script can run several rigid mode, pw_rigid modes, or other variations of the parameters. 

Strides size and overlap are relevant parameters if 'pw_rigid' mode is use. 

Optimal parameters correspond to --- values of crispness. The script prints in the screen the full dictionary corresponding to the optimal ones. 

(Ideally this should be save in the parameters data base, to be read later and directly implemented in the pipeline...for now, let's do it manually)

Ready for next step.

2. Equalizer


Equalization is the last added step in the pipeline. This step should be run with to main objectives. 

1) Reducing the photobleching effect (as mentioned in https://imagej.net/Bleach_Correction )
2) Have an stationary signal over different days of recoding. 

The equalization procedure consists in mapping the histogram of pixel values of one image into the histogram of the other image (see '/home/sebastian/Documents/Melisa/calcium_imaging_analysis/data/interim/equalizer/meta/figures/MonaLisa_Vincent_equalization.png' for exemplification). 

It is known that bleaching diminishes the quality of the image turning it darker. The idea here is to use the histogram of the first recording videos as a template to match the following videos.


26/11/2019 This function is under construction but will take parameters: 

	parameters_equalizer = {'make_template_from_trial': '6_R', 'equalizer': 'histogram_matching', 'histogram_step': h_step}



	#make_template_from_trial < - Refers to where the reference histogram is taking from
	#equalizer <- Sets the method to histogram matching. There are other methods like exponential fitting, or mean equalization that are not yet implemented. 
	#a part of the code also produces the histogram to verify equalization. h_step sets the sets the bin size


4. Source extraction
---

Source extraction can be run both in the server and it the local machine, but again it is supposed to be faster in the server. 

Source extraction parameters as specified in a dictionary as follows : 

	parameters_source_extraction ={'session_wise': False, 'fr': 10, 'decay_time': 0.1, 'min_corr': 0.77, 'min_pnr': 6.6,
                               'p': 1, 'K': None, 'gSig': (5, 5), 'gSiz': (20, 20), 'merge_thr': 0.7, 'rf': 60,
                               'stride': 30, 'tsub': 1, 'ssub': 2, 'p_tsub': 1, 'p_ssub': 2, 'low_rank_background': None,
                               'nb': 0, 'nb_patch': 0, 'ssub_B': 2, 'init_iter': 2, 'ring_size_factor': 1.4,
                               'method_init': 'corr_pnr', 'method_deconvolution': 'oasis',
                               'update_background_components': True,
                               'center_psf': True, 'border_pix': 0, 'normalize_init': False,
                               'del_duplicates': True, 'only_init': True}

	# parameters for source extraction and deconvolution
	p = 1               # order of the autoregressive system
	K = None            # upper bound on number of components per patch, in general None
	gSig = (3, 3)       # gaussian width of a 2D gaussian kernel, which approximates a neuron
	gSiz = (13, 13)     # average diameter of a neuron, in general 4*gSig+1
	Ain = None          # possibility to seed with predetermined binary masks
	merge_thr = .7      # merging threshold, max correlation allowed
	rf = 40             # half-size of the patches in pixels. e.g., if rf=40, patches are 80x80
	stride_cnmf = 20    # amount of overlap between the patches in pixels
	#                     (keep it at least large as gSiz, i.e 4 times the neuron size gSig)
	tsub = 2            # downsampling factor in time for initialization,
	#                     increase if you have memory problems
	ssub = 1            # downsampling factor in space for initialization,
	#                     increase if you have memory problems
	#                     you can pass them here as boolean vectors
	low_rank_background = None  # None leaves background of each patch intact,
	#                     True performs global low-rank approximation if gnb>0
	gnb = 0             # number of background components (rank) if positive,
	#                     else exact ring model with following settings
	#                         gnb= 0: Return background as b and W
	#                         gnb=-1: Return full rank background B
	#                         gnb<-1: Don't return background
	nb_patch = 0        # number of background components (rank) per patch if gnb>0,
	#                     else it is set automatically
	min_corr = .8       # min peak value from correlation image
	min_pnr = 10        # min peak to noise ration from PNR image
	ssub_B = 2          # additional downsampling factor in space for background
	ring_size_factor = 1.4  # radius of ring is gSiz*ring_size_factor


The script runs different selections of gSig (gSiz = 4 * gSig + 1) and saves the resulting corr and pnr summary image (as well as the combination) in '/home/sebastian/Documents/Melisa/calcium_imaging_analysis/data/interim/source_extraction/trial_wise/meta/figures/corr_pnr/'. From here exploration of the effect of different gSig on the summary images can be done.

Later exploration of min_corr and min_pnr can be done. The scripts creates a bunch of histograms to get an intuitive idea of the values of corr and pnr of a video. 

Selection a range or corr and pnr values, there is a script that computes source extraction and plots a figure with all the contour plots for the analyzed parameters. This plotting helps to understand how initial seeds change the final result of the extraction.  

Visual inspection of this figure can help to decide which values of pnr and corr are adequate (selects the higher number of neurons that are 'real neurons'). Parameter selection of source extraction can be done in combination with parameter selection of component evaluation in order to get a better solution.

5. Component Evaluation

For running component evaluation, a version of all previous states should be chosen. This step can be run in the local machine.

Component evaluation parameters are 3: minimal value of signal to noise ration in the extracted calcium traces, minimal pearson correlation coefficient between a template for the 'neuron footprint' extracted as a mean over all the images and the source extracted footprint, and a boolean variable that specifies whether the assessment will or will not use a CNN for classification as good/bad component.

Componente evaluation parameters as specified in a dictionary as follows : 

	parameters_component_evaluation = {'min_SNR': 3,
                                   'rval_thr': 0.85,
                                   'use_cnn': False}


The script proposed runs for all source extraction versions selected different selection of component evaluation parameters and makes some plot where accepted and rejected components can be seen in the contour plot and also in  the traces ( plot_contours_evaluated and plot_traces_multiple_evaluated)



PARAMETER SELECTION FOR SESSION (or day) WISE ANALYSIS
==================================================================
Bleaching effect ===> Because of continuous exposure to the light of the microscope, the image gets bleached. The effect of the bleaching can be seen in figures save in the folders: 

'/home/sebastian/Documents/Melisa/calcium_imaging_analysis/data/interim/source_extraction/session_wise/meta/figures/contours/'

and

'/home/sebastian/Documents/Melisa/calcium_imaging_analysis/data/interim/source_extraction/trial_wise/meta/figures/fig:corrpnrphotobleaching56165.png'

In the later the mean of the pnr and correlation image is taken, and the figure shows the evaluation of the mean value over different days (Implement this as a separate part from source extraction to make it faster, because only the corr and pnr image are required).

In the first folder, there are different figures that shows the bleaching effect in the source extraction by showing the contours of different trials within a day and using different parameters for the extraction. Visual inspection can show that during the first trial more neurons are detected, and during late trials that corr image gets blurier and less amount of neurons are detected. 

Next problem is then how to select source extraction parameters (and component evaluation parameters) that are useful for all the days, and that are chosen with a 'good enough' criteria. Three different paths are being explore here. 


1) Looking for optimality ===> Use the same parameter for every day. For doing so, source extraction and component evaluation is performed in a small cropped region with a intensive parameter exploration. First idea is to select the source extraction and component evaluation parameters that near to a border transition from accepting every neuron in the last day maximize the cell counting in the last days while minimizing the false positive counting (using same source extraction and component evaluation params) in the first day.

2) Template matching ===> Based on the assumption that if a neuron is detected during the first trial, it should also be there for the later trials, once parameters for the first trial has been selected, use the neural footprint as a template for the other trials. For this, make an exploration of source extraction parameters during every other trial and select the one that maximizes overlapping or matching between new selected neurons and the ones selected during the first trial.

3) Source extraction of all the trials together. For this, it is important first to do the alignment between different trials. Use motion correction for alignment and run source extraction and component evaluation for everything together. NEW SUGGESTION BT F.STELLA (IN DEVELOPMENT): Verify source extraction during resting periods. 	 



Finally, compare the results ==> compare cell counting, final footprints and calcium traces ! 

----
Running the complete pipeline step by step 
----

Once parameters are selected...

Run the pipeline with all the steps with the best parameters, all the trials and resting/non resting conditions, all the sessions, all the mice! HAPPY :) 


