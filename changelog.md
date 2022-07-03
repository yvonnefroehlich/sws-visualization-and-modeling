# Changelog


## Release v1.0 (2022/07/XX)

Changes and differences compared to [sws_tools](https://github.com/michaelgrund/sws_tools) by [Michael Grund](https://github.com/michaelgrund)
- [001_stereoplot]()
   - support of _rotation-correlation method_ ([**_Bowman & Ando 1987_**]( https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.))
   - support of _eigenvalue method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
   - support of splits of _simultaneous inversion of multiple waveforms_ ([**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899))
   - fix colormap set up for _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862))
   - support of cyclic colormaps of _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862))
   - remove diverging and multi-sequential colormaps of _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862))
   - support of sequential and cyclic colormaps of _cmocean colormaps_ ([**_Thyng et al. 2016_**](http://dx.doi.org/10.5670/oceanog.2016.66))
   - add check whether files / data is from one single recording station
- [003_modeling]()
   - support of input files of _SplitLab_ version 1.2.1 (Porritt 2014) (additionally to version 1.0.5 (original))
   - sort models separately for the different model types (based on rote mean square error)
   - no double-counting of stacks and singles (within 5° BAZ bins)
   - allow for color-coding with respect to the fast polarization direction (support of _phase_ of _cmocean colormaps_ ([**_Thyng et al. 2016_**](http://dx.doi.org/10.5670/oceanog.2016.66)))
   - allow for color-coding with respect to the root mean square error (support of _grayC_ of _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862)))
   - read multi-event analysis results from txt files instead of from structure
   - support of splits of _simultaneous inversion of multiple waveforms_ ([**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899)) (nulls are only plotted)
   - generate model parameter plots for the different model types
   - generate synthetic stereoplots separately for the different model types
   - add expected BAZ directions for nulls for the two-layer scenario (under development)
   - output txt files with model parameters for the different model types

**Contributors**: [Yvonne Fröhlich](https://github.com/yvonnefroehlich)
