### Reducing and meta-analysing estimates from distributed lag non-linear models

------------------------------------------------------------------------

An illustration of the method for reducing estimates of bi-dimensional exposure-lag-response associations obtained by distributed lag non-linear models (DLNMs) from multiple studies, and then pooling them. The example reproduces the example included in the article:

Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from distributed lag non-linear models. *BMC Medical Research Methodology*. 2013;**13**(1):1. DOI: 10.1186/1471-2288-13-1. PMID: 23297754. [[freely available here](http://www.ag-myresearch.com/2013_gasparrini_bmcmrm.html)]

The code uses functions in the R packages [dlnm](https://github.com/gasparrini/dlnm) and [mvmeta](https://github.com/gasparrini/mvmeta).

------------------------------------------------------------------------

The code:

-   *regEngWales.csv* stores the daily time series data from 10 locations corresponding to regions of England and Wales in the period 1993–2006
-   the numbered files from *01.prep.R* to *06.metareg.R* reproduce the results of the illustrative example

Download as a ZIP file using the green button *Clone or download* above
