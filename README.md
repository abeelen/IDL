
Several small [IDL][IDL] that I wrote for my own use but that may be usefull to other. They are all under the [GPL][GPL], so everyone is permitted to redistribute and/or modify them under the terms of the [GPL][GPL]


* [GILDAS][GILDAS] Line fitting (gildas):
> Deriving line parameter from a UVFIT in the [Grenoble Image and Line Data Analysis System][GILDAS] is not an easy task. I developped a few program to be able to do it easily. First you should patch the gio library of [GILDAS][GILDAS] in order to be able to save uvfit result in fits format. The task gildas_fits or fits will then be able to save UVFITS result table in fits format (use standard fits format and number of bits=-32). The resulting fits file can then be used to derive line parameters, for example with the IDL script included, using the Markwardt IDL Library to make the fit. The very usefull TexToIDL library is also used to legend the plot. Uncertainties on the fitted line parameters are derived from spectra uncertainties, computed from the visibilities and thus does not required the complex process of deconvolution/reconvolution involved in the mapping of interferometric data.

    - ```tofits.diff``` : small patch for the tofits.f file of the gio library of GILDAS
    - ```fit_line.pro``` : fit a gaussian line to the spectrum
    - ```mp_gauss.pro``` : define the fitted function, a gaussian with a continnum


* [IDL][IDL] [ASURV][ASURV] (asurv):
> To easily use and test the [Astronomical SURVival Analysis package][ASURV], I wrote a simple C wrapper which can be directly be called, as a library, from [IDL][IDL] or other language/program. This work has been highly simplified by the splitting of the original monolitical code of [ASURV][ASURV] by the Starlink project. The univariate two sample tests, bivariate correlation tests and linear regression have been wrapped.
    - ```asurv.tar.gz``` : The, now gone, starlink version of [ASURV][ASURV]
    - ```lib_bivar.f  lib_univar.f  lunivar.c``` : actual wrapper
    -  ```Makefile``` : makefile to build the wrapper

* [IDL][IDL] Statistics (stats)
> Two very usefull routines to plot the histogram and the partition function of a dataset. The first use the IDL  histogram function but also compute the proper abscissa.

    - ```plot_repart.pro``` : draw an histogram of a variable
    - ```plot_histo.pro``` : plot the repartition function of a variable


* [IDL][IDL] [SDSS][SDSS] (sdss)
> I wrote some [IDL][IDL] routine to easily manipulate [SDSS][SDSS] spectra. It has been tested on the [DR3QSO][DR3QSO] catalog and spectra.

    - ```filename_sdss.pro``` : construct the filename of a spectra from information in the [DR3QSO][DR3QSO] catalog
    - ```read_sdss_spectrum.pro``` : read an [SDSS][SDSS] spectrum into a [IDL][IDL] structure
    - ```plot_sdss_spectrum.pro``` : plot an [SDSS][SDSS] spectrum the [SDSS][SDSS] way

* [IDL][IDL] Misc. (mis)
> Here are a few [IDL][IDL] scripts that I wrote for every day use.
    - ```zeta.pro``` compute the zeta function
    - ```barre.pro``` draw an horizontal bar
    - ```const.pro``` define some physical constant (SI)
    - ```copyright.pro``` add copyright to your plot
    - ```draft.pro``` add a draft comment to your plot

[GPL]: http://www.gnu.org/licenses/gpl-3.0.txt  "The GPL v3 License"
[IDL]: http://www.exelisvis.com/docs/using_idl_home.htmlscripts "The Interactive Data Language"
[SDSS]: http://www.sdss.org/ "The Sloan Digital Sky Survey"
[DR3QSO]: http://www.sdss.org/dr4/products/value_added/qsocat_dr3.html "The SDSS QSO Data Release 3"
[GILDAS]: https://www.iram.fr/IRAMFR/GILDAS/https://www.iram.fr/IRAMFR/GILDAS/ "The IRAM GILDAS Software"
[ASURV]: http://astrostatistics.psu.edu/statcodes/sc_censor.html "The ASURV software"
