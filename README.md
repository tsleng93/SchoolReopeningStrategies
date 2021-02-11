# SchoolReopeningStrategies

This repository contains the underlying code and functions used in the paper "Assessing the impact of secondary school reopening strategies on within-school Covid-19 transmission and absences: a modelling study", by Trystan Leng*, Edward M. Hill, Robin N. Thompson, Michael J. Tildesley, Matt J. Keeling and Louise Dyson.

Corresponding author: Trystan Leng, email: trystan.leng@warwick.ac.uk

Underlying code, functions and visualisations are written in matlab version 2019b.


# Model functions
The function 'interactingyeargroups.m' is the main function of the model, simulating the spread of infection in a secondary school over the course of a half term.

The functions 'Modeloutputs.m' and 'Moremodeloutputs.m' extract key quantities from the 'history' output of the schools model, such as prevalence, total infected over a half term, absences, and testing volume. The specific outputs are described in comments within the .m files.



# Regenerating results
To generate the underlying data and plots for Figures 1 and 2, users should run the code 'main_plots.m'. Figures 1 (b), (c), and (d) require [Violin plots for Matlab](https://github.com/bastibe/Violinplot-Matlab#:~:text=A%20violin%20plot%20is%20an,overlays%20the%20data%20points%20itself.&text=violinplot%20is%20meant%20as%20a,boxplot%20(excluding%20named%20arguments).). Figure 1 (e) requires [Break Y Axis by MikeCF (2021)](https://www.mathworks.com/matlabcentral/fileexchange/45760-break-y-axis)

To generate the plots for the PCR and LFT test probability profiles (Supplementary Figure S1), users should run the code 'test_probability_profiles.m'. Test probability profiles for symptomatic individuals were obtained directly from [Hellewell et al. (2020)](https://cmmid.github.io/topics/covid19/pcr-positivity-over-time.html). Test probability profiles for asymptomatic individuals were obtained by assuming that profiles were equal to that of symptomatic individuals until the peak of infection, but then decay more rapidly. 

To generate the underlying data and plots for the sensitivity analyses (Supplementary Figures S2-S4), users should run the code 'sensitivity_analyses.m'

To generate the underlying data and plots for Supplementary Figuure 5, users should run the code 'suppfigure_K.m'.
