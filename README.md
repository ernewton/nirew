# Equivalent widths from NIR spectra

This code measures equivalent widths (EWs) from spectra, with an emphasis on NIR absorptiion features. It also can be used to infer stellar parameters, including metallicity, spectral type, temperature, radius, and luminosity based on established relationships with spectral features.

## Important notes

The H-band feature definitions in Table 1 in Newton et al. (2015) are from a previous version of this code. The feature definitions included here are correct, and are consistent with the rest of the paper.

Minor modifications from the original version of the radial code improve robustness, but result in small differences in the measured EWs and the inferred stellar properties. These differences are within estimated errors.

## Routines

### example
Demonstrates use of `ern_rv` and `measure_ew` on spectra, showing the effect of oversampling and spectral resolution.

### jl_example
Demonstrates use of `nir_rv` and `measure_ew` to calculate radial velocity and various EWs.

### measure_na

Measure the Na EW from a spectrum and use it to estimate metallicity using relation from Newton et al. (2014). Errors (1-sigma) can also be calculated. See Newton et al. (2014, 2015) for notes on where this relation is applicable. Also calculates H2O-K2 index and uses it to estimate NIR spectral type (Rojas-Ayala et al. 2012, Newton et al. 2014).

### measure_ew

Measures an EW using a pseudo-continuum and the trapezoidal rule using the routine `tsum`. It is recommended that you oversample your spectrum.

## Dependencies

idlutils is required.

Code to calculate EWs directly has no other dependencies. However, spectra must be shifted to rest wavelengths before EWs can be calculated. Higher-level code and the examples use `tellrv` to do this, which is [available on github](https://github.com/ernewton/nirew).

## Reference

If you use this code in published research, please cite [Newton et al. (2014)](http://adslabs.org/adsabs/abs/2014AJ....147...20N/)

If temperature, radius, or luminosity is used, please cite [Newton et al. (2015)](http://adslabs.org/adsabs/abs/2015ApJ...800...85N/)

## License

Copyright 2015 Elisabeth R. Newton. Licensed under the terms of the MIT License (see LICENSE).