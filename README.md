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

# ghg-flux-expert