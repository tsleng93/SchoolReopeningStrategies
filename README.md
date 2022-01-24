


# SchoolReopeningStrategies



This repository contains the underlying code and functions used in the paper "Assessing the impact of secondary school reopening strategies on within-school Covid-19 transmission and absences: a modelling study", by Trystan Leng*, Edward M. Hill, Robin N. Thompson, Michael J. Tildesley, Matt J. Keeling and Louise Dyson.

Corresponding author: Trystan Leng, email: trystan.leng@warwick.ac.uk

Underlying code, functions and visualisations are written in matlab version 2019b.


Also contained in this The underlying model and model fitting code used in the paper "Quantifying pupil-to-pupil SARS-CoV-2 transmission and the impact of lateral flow testing in English secondary schools" by Trystan Leng*, Edward M. Hill, Alex Holmes, Emma Southall, Robin N. Thompson, Michael J. Tildesley, Matt J. Keeling, and Louise Dyson is contained in the subfolder 'FittedModel'. The read me file for this is included below, in the section 'ReadMe - Quantifying within-school SARS-CoV-2 transmission'.

# ReadMe -  Assessing the impact of secondary school reopening strategies 


## Model functions
The function 'interactingyeargroups.m' is the main function of the model, simulating the spread of infection in a secondary school over the course of a half term. 

The function 'Interactingyeargroups_externalinfection.m' is a function used within the main function to simulate external infections to pupils.

The functions 'Modeloutputs.m' and 'Moremodeloutputs.m' extract key quantities from the 'history' output of the schools model, such as prevalence, total infected over a half term, absences, and testing volume. The specific outputs are described in comments within the .m files.

## Regenerating results
To generate the underlying data and plots for Figures 1 and 2, users should run the code 'main_plots.m'. Figures 1 (b), (c), and (d) require [Violin plots for Matlab](https://github.com/bastibe/Violinplot-Matlab#:~:text=A%20violin%20plot%20is%20an,overlays%20the%20data%20points%20itself.&text=violinplot%20is%20meant%20as%20a,boxplot%20(excluding%20named%20arguments).). Figure 1 (e) requires [Break Y Axis by MikeCF (2021)](https://www.mathworks.com/matlabcentral/fileexchange/45760-break-y-axis)

To generate the plots for the PCR and LFT test probability profiles (Supplementary Figure S1), users should run the code 'test_probability_profiles.m'. Test probability profiles for symptomatic individuals were obtained directly from [Hellewell et al. (2020)](https://cmmid.github.io/topics/covid19/pcr-positivity-over-time.html). Test probability profiles for asymptomatic individuals were obtained by assuming that profiles were equal to that of symptomatic individuals until the peak of infection, but then decay more rapidly. 

To generate the underlying data and plots for the sensitivity analyses (Supplementary Figures S3-S6), users should run the code 'sensitivity_analyses.m'

To generate the underlying data and plots for Figure 3, and Supplementary Figures S7-S11) users should run the code 'additional_plots.m'.

# ReadMe - Quantifying pupil-to-pupil SARS-CoV-2 transmission 

The model code is written in Matlab, and has been tested on Matlab 2019a and Matlab 2021a. No non-standard hardware is required to run the model code. 

Relevant scripts and functions are contained in the subfolder 'FittedModel'. By downloading this subfolder, a user will be able to run the relevant model functions. This download should be completed in a few minutes.  We extend the model described in 'interactingyeargroups.m' to incorporate greater aspects of realism.

The main analyses in the paper require data not publicly available. However, we provide comparable data in the file 'exampleworkspace.mat'. In particular:

1) We used modelled secondary school sizes and close contact group sizes, rather than those directly obtained from data
2) We include sigmoidal fits of the relative frequency of B.1.1.7 variant for each LTLA, rather than the raw data regarding S-gene negatives
3) Smoothed community PCR positive testing rates for each LTLA are provided, obtained using the Matlab function 'smooth' with the setting 'rloess, rather than the raw data obtained from Pillar 2 data.
4) We assume 36% uptake LFT uptake across all LTLAs, rather than obtaining uptake levels from Pillar 2 data that vary through time and by LTLA.

Comparable analyses can be obtained by running the script 'LTLAComparisonMain.m'.  Comparable figures to figures 2 and 3 (omitting data) can be obtained through running the script 'ExampleFigures.m'.  This runs the baseline strategy of twice weekly mass testing and the isolation of close contacts. Different school control strategies can be obtained by adjusting the 'Strategy' structure:

 - to run the 'isolating close contacts' strategy, set Strategy.masstesting = 0 , Strategy.isolation = 1, Strategy.SCT = 0
 - to run the 'mass testing' strategy (i.e. twice weekly mass testing alone), set Strategy.masstesting = 2, Strategy.isolation = 0, Strategy.SCT = 0
 - to run the 'mass testing + serial contact testing' strategy, set Strategy.masstesting = 2, Strategy.isolation = 0, Strategy.SCT = 1

LTLAComparisonMain.m should take ~45 minutes to run, using a parallel pool of 10 workers.

## Model functions
The function 'interactingyeargroupsextended.m' is the main function of the model, simulating the spread of infection in a secondary school from 31st August 2020 until the 23rd May 2021.

The function 'interactingyeargroupsquicker.m' is a quicker version of the above model, owing to a  quicker within-model infection process, but does not track within-school R explcitly. This function is used for the model fitting.

The functions 'Interactingyeargroups_infection.m' and 'Interactingyeargroups_externalinfection.m' are functions that used to simulate parts of the infection process in interacting year groups. These have been converted to .mex files to speed up model code. 

The function 'SchoolPopulation.m' is a function creating the matrix structure of school contacts for the longer version of the code, while 'SchoolPopulationQuicker.m' is a function creating the matrix structure of school contacts for the quicker version of the code.

The function 'LTLAComparisonMain.m' is the function used to obtain the main model results. 

The function 'FittingModelFunction.m' is a function containing the fitting procedure used in the model.

The functions 'Modeloutputscondensed.m' extracts key quantities from the 'history' output of the schools model. The specific outputs are described in comments within the .m files.


