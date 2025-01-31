# Changelog

## Release v2.0 (2024/MM/DD) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13799086.svg)](https://doi.org/10.5281/zenodo.13799086)

- General
   - Improve checks for input SWSMs data ([PR #11](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/11))
   - Improve coding style ([PR #27](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/27))
   - Add section "Related topics" to main README ([PR #25](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/25))
   - Move "Contents" section to separate READMEs & introduce "_image folder ([PR #31](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/31))

- [001_stereoplot](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/001_stereoplot)
   - Add selection / querries for phase and observation type ([PR #23](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/23))
   - Allow plotting recording station in the center of the stereoplot ([PR #33](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/33))
   - Update to R2024a ([PR #21](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/21))
   - Update regarding crameri v1.09 ([PR #15](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/15))
   - Use function "crameri.m" instead of struct "CrameriColourMaps*.mat" to check existence of SCM ([PR #32](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/32))

- [002_visualization](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization)
   - Add loading function for datasets on RADAR4KIT in MATLAB, Python, and R ([PR #24](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/24))

- [003_modeling](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling)
    - Fixes regarding colormap and colorbar for phi ([PR #13](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/13))
    - Fix query for plotting not considered BAZ range in gray ([PR #14](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/14))
    - Improve visualization of stats regarding the model type ([PR #22](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/22))
    - Save structs for model types separately ([PR #16](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/16))
    - Update table for required disc space in modeling README ([PR #17](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pull/17))



**Contributors**: [Yvonne Fröhlich](https://github.com/yvonnefroehlich)

For changes and differences compared to [sws_tools](https://github.com/michaelgrund/sws_tools) by [Michael Grund](https://github.com/michaelgrund) see
the [changelog entry for release v1.0](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/blob/main/changelog.md#release-v10-20221016-)


## Release v1.0 (2022/10/16) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7213157.svg)](https://doi.org/10.5281/zenodo.7213157)

Changes and differences compared to [sws_tools](https://github.com/michaelgrund/sws_tools) by [Michael Grund](https://github.com/michaelgrund):
- [001_stereoplot](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/001_stereoplot)
   - Support of _rotation-correlation method_ ([**_Bowman & Ando 1987_**](https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.))
   - Support of _eigenvalue method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
   - Support of _splits_ of _simultaneous inversion of multiple waveforms_ ([**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899))
   - Fix colormap set up for _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862)) ([PR](https://github.com/michaelgrund/sws_tools/pull/4) opened, will be fixed)
   - Support of cyclic colormaps of _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862))
   - Remove diverging and multi-sequential colormaps of _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862))
   - Support of sequential and cyclic colormaps of _cmocean colormaps_ ([**_Thyng et al. 2016_**](http://dx.doi.org/10.5670/oceanog.2016.66))
   - Add check whether files / data is from one single recording station
- [003_modeling](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling)
   - Support of input files of _SplitLab_ version 1.2.1 (Porritt 2014) (additionally to version 1.0.5 (original))
   - Sort models separately for the different model types (based on rote mean square error)
   - No double-counting of stacks and singles (within 5° backazimuth bins)
   - Allow for color-coding with respect to the fast polarization direction (support of _phase_ of _cmocean colormaps_ ([**_Thyng et al. 2016_**](http://dx.doi.org/10.5670/oceanog.2016.66)))
   - Allow for color-coding with respect to the root mean square error (support of _grayC_ of _Scientific Colour maps_ ([**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862)))
   - Read multi-event analysis results from txt files instead of from structure
   - Support of _splits_ of _simultaneous inversion of multiple waveforms_ ([**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899)) (_nulls_ are only plotted)
   - Generate model parameter plots for the different model types
   - Generate synthetic stereoplots separately for the different model types
   - Add expected backazimuth directions for _nulls_ for the two-layer scenario (under development)
   - Output txt files with model parameters for the different model types

**Contributors**: [Yvonne Fröhlich](https://github.com/yvonnefroehlich)
