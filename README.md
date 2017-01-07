
### Updated R code and data from Gasparrini & Armstrong BMCmrm 2013

--------------------------------------------------------------------------------

An illustration on methods for reducing estimates of bi-dimensional exposure-lag-response associations obtained by DLNMs  from multiple studies, and then pooling them. The example reproduces the example included in the paper:

Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from distributed lag non-linear models. *BMC Medical Research Methodology*. 2013;**13**(1):1. [[freely available here](http://www.ag-myresearch.com/2013_gasparrini_bmcmrm.html)]

The code uses functions in the R packages [R package dlnm](https://github.com/gasparrini/dlnm) and [R package mvmeta](https://github.com/gasparrini/mvmeta).

--------------------------------------------------------------------------------

The code:

  * *regEngWales.csv* stores the daily time series data from 10 locations corresponding to regions of England and Wales in the period 1993–2006
  * the numbered files from *01.prep.R* to *06.metareg.R* reproduce the results of the illustrative example

Download as a ZIP file using the green button *Clone or download* above