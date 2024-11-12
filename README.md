# restore4cs ghg flux expert
This repository contains the scripts to track expert manual assessment of CO2 and CH4 incubation timeseries.
Experts are prompted to select what look like safe data to calculate CO2 flux. For CH4, experts are asked to manually select what looks to them like diffusion only (or no ebullition).

Then, fluxes are being processed with different approaches inspired from the R package GoFlux (https://github.com/Qepanna/goFlux).
Data and fluxes estimates are saved on Dropbox, not on the remote.

<img src="https://github.com/camilleminaudo/restore4cs-scripts/blob/main/RESTORE4Cs_LOGO_DEF.jpg" width=50% height=50%>


## Install
First, you'll need to make sure you have R and R Studio correctly installed.
You need to make sure you have access to the RESTORE4Cs Dropbox.

Then, you will need to install the GoFlux package, which is stored on GitHub.
To install a package from GitHub, one must first install the package
`devtools` from the CRAN:

``` r
if (!require("devtools", quietly = TRUE))
install.packages("devtools")
```

Then, install the `goFlux` package from GitHub:

``` r
try(detach("package:goFlux", unload = TRUE), silent = TRUE)
devtools::install_github("Qepanna/goFlux")
```

Then, download or ideally clone the ghg-flux-expert repository, and it should be enough to have the scripts working.
The main 

The functioning of the package depends on many other packages
(`data.table`, `dplyr`, `ggnewscale`, `ggplot2`, `graphics`,
`grDevices`, `grid`, `gridExtra`, `lubridate`, `minpack.lm`, `msm`,
`pbapply`, `plyr`, `purrr`, `rjson`, `rlist`, `SimDesign`, `stats`,
`tibble`, `tidyr`, `utils`), which will be installed when installing
`GoFluxYourself`. If prompted, it is recommended to update any
pre-installed packages.

# Run the scripts
When ready, open scripts/expert_vs_automatic_processing.R
Enter your username (Firstname + Lastname) and the correct path to the Dropbox folder.
You can also change the number of incubation timeseries you want to work on in this session. I recommend to keep nb_draw at a value of 10.
Then run the script.

1.	The script loads all the fieldsheets + the temperature logger info. We are here only working on measurements performed with the floating chamber, i.e. the open water strata. No vegetation, no transparent chambers.
2.	10 incubations are randomly chosen.
3.	You’ll be prompted in a new R window to choose for each of these incubations what are to you the safe start and end time to compute the CO2 flux. To do so, you simply need to click on the graphic. ONE click for the start, ONE click for the end. When you’ve done so, the same incubation will be shown but showing this time the selected data.
4.	Then the script calculates CO2 and CH4 fluxes according to your timeframe selection
5.	The script calculates CO2 and CH4 fluxes blindly, without any expert assessment.
A bunch of plots will appear, showing the automatized computation of CH4 diffusion vs bubbling.
6.	To estimate the diffusive CH4 flux, you’ll be prompted to chose on a pCH4 plot what looks to you like the start and end times of a purely diffusive behaviour. You can only select one part of the plot, multiple selection is at the moment impossible.
7.	The script computes CH4 ebullitive and diffusive fluxes based on your selection.
8.	Everything is saved in a few csv files and a figure is displayed just for fun.

IMPORTANT: if for a given incubation you think none of the timeseries shown can be used (erratic or anomalous data), select a very small section of it (a few seconds), this will disable the flux calculation and help me track suspicious measurements.

You can play again and again and again. Going through 10 incubations takes about 2 minutes. You can chose to process more than 10 incubations at once, but be aware that this is a boring game, sorry!.
I recommend keeping the number of incubations to process (nb_draw) at 10, but you can repeat your experience as much as you want by re-running the entire program.
