# Visualization and Modeling of Shear Wave Splitting [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7213156.svg)](https://doi.org/10.5281/zenodo.7213156)

_MATLAB_ functions for visualization and modeling of shear wave splitting observations:

- Optimized for the output of _SplitLab_ ([**_Wüstefeld et al. 2008_**](https://doi.org/10.1016/j.cageo.2007.08.002)) and
  [_StackSplit_](https://github.com/michaelgrund/stacksplit) ([**_Grund 2017_**](https://doi.org/10.1016/j.cageo.2017.04.015)).
- The modeling routine is applicable for the _energy minimization_ method ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899)).
- Extended and strongly modified from [sws_tools](https://github.com/michaelgrund/sws_tools) by [Michael Grund](https://github.com/michaelgrund).



## Citation

_Zenodo DOI for all / latest version_: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7213156.svg)](https://doi.org/10.5281/zenodo.7213156)

If you make use of this material, please acknowledge the relating publications in which framework these functions were implemented:

**Peer-reviewed journal article**
- [**_Fröhlich Y., Grund M., Ritter J. R. R. (2024)_**](https://doi.org/10.1093/gji/ggae245).
  Lateral and vertical variations of seismic anisotropy in the lithosphere-asthenosphere system underneath Central Europe from long-term splitting measurements.
  *Geophysical Journal International*, 239(1), 112-135.
  https://doi.org/10.1093/gji/ggae245.
- [**_Ritter J. R. R., Fröhlich Y., Sanz Alonso Y. & Grund M. (2022)_**](https://doi.org/10.1007/s10950-022-10112-w).
  Short-scale laterally varying SK(K)S shear wave splitting at BFO, Germany – implications for the determination of anisotropic structures.
  *Journal of Seismology*, 26, 1137-1156.
  https://doi.org/10.1007/s10950-022-10112-w, correction https://doi.org/10.1007/s10950-023-10136-w.
- [**_Grund M., Ritter J. R. R. (2020)_**](https://doi.org/10.1093/gji/ggaa388).
  Shear-wave splitting beneath Fennoscandia – evidence for dipping structures and laterally varying multilayer anisotropy.
  *Geophysical Journal International*, 223(3), 1525-1547.
  https://doi.org/10.1093/gji/ggaa388.

**Presentations**
- [**_Fröhlich Y., Ritter J. R. R. (2024)_**](http://dx.doi.org/10.5281/zenodo.14510993).
  Vertical and Small-scale Lateral Varying Seismic Anisotropy in the Upper Mantle Underneath the Upper Rhine Graben, Central Europe.
  *Annual Meeting of the American Geophysical Union*, Washington D.C..
  Division Session Exploring Innovations and New Directions in Seismic Anisotropy and Attenuation: Observations, Models, and Experiments I Oral, DI21A-02.
  [Abstract ID 1578275](https://agu.confex.com/agu/agu24/meetingapp.cgi/Paper/1578275).
  http://dx.doi.org/10.5281/zenodo.14510993.



## Content

_Example figures_: Generated with the provided test data

- **[001_stereoplot](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/001_stereoplot)**
- **[002_visualization](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization)**
- **[003_modeling](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling)**

![](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/raw/main/_images/000_repo_readme_image.png)



## Requirements

_Tested with_: R2022a, R2021a,b under Linux and Windows

- **_MATLAB_**: Forward calculation
  - Deep Learning Toolbox
  - Mapping Toolbox
  - [_MATLAB Seismic Anisotropy Toolkit_ (MSAT)](https://www1.gly.bris.ac.uk/MSAT/) ([**_Walker & Wookey 2012_**](https://doi.org/10.1016/j.cageo.2012.05.031))
- **Data**: Shear wave splitting observations
  - Output *.txt files (_nulls_, _splits_) of _SplitLab_ version 1.5.0 ([**_Wüstefeld et al. 2008_**](https://doi.org/10.1016/j.cageo.2007.08.002)) or 1.2.1 (**_Porritt 2014_**)
  - Output *.mat structure and *.txt files (STACK, SIMW) of _StackSplit_ ([**_Grund 2017_**](https://doi.org/10.1016/j.cageo.2017.04.015))
- **Colormaps** (optional): Color-coding of the fast polarization direction and the root mean square error
  - [MatPlotLib Perceptually Uniform Colormaps](https://de.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps)\
    version v2.1.3, MATLAB File Exchange, last access 2022 June 26
  - [crameri perceptually uniform scientific colormaps](https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)\
    version v1.09, MATLAB File Exchange, last access 2023 April 10; based on [**_Crameri (2021)_**](https://zenodo.org/record/5501399)
  - [cmocean perceptually-uniform colormaps](https://de.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)\
    version v2.02, MATLAB File Exchange, last access 2022 June 18; based on [**_Thyng et al. (2016)_**](http://dx.doi.org/10.5670/oceanog.2016.66)



## Releases

| Release | Zenodo DOI | Publication | RADAR4KIT dataset |
| --- | --- | --- | --- |
| [dev](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main) |  |  |  |
| [v2.0](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/releases/tag/v2.0) | | [**_Fröhlich et al. (2024)_**](https://doi.org/10.1093/gji/ggae245) | https://dx.doi.org/10.35097/685 |
| [v1.0](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/releases/tag/v1.0) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7213157.svg)](https://doi.org/10.5281/zenodo.7213157) | [**_Ritter et al. (2022)_**](https://doi.org/10.1007/s10950-022-10112-w) | https://dx.doi.org/10.35097/684 |

For details of the individual releases as well as for changes and differences compared to [sws_tools](https://github.com/michaelgrund/sws_tools)
by [Michael Grund](https://github.com/michaelgrund) see the [changelog](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/blob/main/changelog.md).



## Known Issues

- Modeling of multi-event analysis: Only using either STACK or SIMW results is supported
- Model parameter distribution for T1: Under development, not fully tested
- Synthetic stereoplot for T1 and H2: Backazimuths of predicted nulls are partly wrong
- Synthetic stereoplot for T1: Gray arrow is partly not exactly placed in the center



## Contributing

For bug reports, suggestions or recommendations feel free to [open an issue](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/issues) or
[submit a pull request](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/pulls) directly here on
[GitHub](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main).



## Related topics

| Software | Language | Description | Author |
| --- | --- | --- | --- |
| [MSAT](https://github.com/andreww/MSAT) | MATLAB | A new toolkit for the analysis of elastic and seismic anisotropy | Andrew M. Walker & James Wookey |
| [PyDRex](https://github.com/seismic-anisotropy/PyDRex) | Python | Simulate crystallographic preferred orientation evolution in polycrystals | Leon Bilton & Thomas Duvernay |
| [PyRaysum](https://github.com/paudetseis/PyRaysum) | Python | Software for modeling ray-theoretical body-wave propagation | Wasja Bloch & Pascal Audet |
| [sws_tools](https://github.com/michaelgrund/sws_tools) | MATLAB | Tools for modeling and plotting of shear-wave splitting data | Michael Grund |



## References

- [**_Bowman, J. R. & Ando, M. (1987)_**](https://doi.org/10.1111/j.1365-246X.1987.tb01367.x).
  Shear-wave splitting in the upper-mantle wedge above the Tonga subduction zone.
  *Geophysical Journal International*, volume 88, issue 1, pages 25-41.
  https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.
- [**_Crameri, F. (2021)_**](https://zenodo.org/record/5501399).
  Scientific colour maps, version 7.0.1. *Zenodo*. http://www.fabiocrameri.ch/colourmaps.php. https://zenodo.org/record/5501399.
- [**_Grund, M. (2017)_**](https://doi.org/10.1016/j.cageo.2017.04.015).
  StackSplit - a plugin for multi-event shear wave splitting analyses in SplitLab.
  *Computers & Geosciences*, volume 105, pages 43-50.
  https://doi.org/10.1016/j.cageo.2017.04.015.
  versions [1.0](https://doi.org/10.5281/zenodo.464385), [2.0](https://doi.org/10.5281/zenodo.7118716), and [3.0](https://doi.org/10.5281/zenodo.5802051)
  available at https://github.com/michaelgrund/stacksplit.
- **_Porritt, R. W. (2014)_**. SplitLab version 1.2.1.
  available at https://robporritt.wordpress.com/software/.
- [**_Restivo, A. & Helffrich, G. (1999)_**](https://doi.org/10.1046/j.1365-246x.1999.00845.x).
  Teleseismic shear wave splitting measurements in noisy environments.
  *Geophysical Journal International*, volume 137, pages 821-830.
  https://doi.org/10.1046/j.1365-246x.1999.00845.x.
- [**_Roy, C., Winter, A., Ritter, J. R. R. & Schweitzer, J. (2017)_**](https://doi.org/10.1093/gji/ggw470).
  On the improvement of SKS splitting measurements by the Simultaneous Inversion of Multiple Waveforms (SIMW).
  *Geophysical Journal International*, volume 208, pages 1508-1523.
  https://doi.org/10.1093/gji/ggw470.
- [**_Silver, P. G. & Chan, W. W. (1991)_**](https://doi.org/10.1029/91JB00899).
  Shear wave splitting and subcontinental mantle deformation.
  *Journal of Geophysical Research*, volume 96, issue B10, pages 16429-16454.
  https://doi.org/10.1029/91JB00899.
- [**_Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M. & DiMarco, S. F. (2016)_**](http://dx.doi.org/10.5670/oceanog.2016.66).
  True colors of oceanography: Guidelines for effective and accurate colormap selection.
  *Oceanography*, volume 29, issue 3, pages 9-13.
  http://dx.doi.org/10.5670/oceanog.2016.66.
- [**_Walker, A. M. & Wookey, J. (2012)_**](https://doi.org/10.1016/j.cageo.2012.05.031).
  MSAT — A new toolkit for the analysis of elastic and seismic anisotropy.
  *Computer & Geosciences*, volume 49, pages 81-90.
  https://doi.org/10.1016/j.cageo.2012.05.031.
  available at https://www1.gly.bris.ac.uk/MSAT/, https://github.com/andreww/MSAT.
- [**_Wolfe, C. J. & Silver, P. G. (1998)_**](https://doi.org/10.1029/97JB02023).
  Seismic anisotropy of oceanic upper mantle: Shear wave splitting methodologies and observations.
  *Journal of Geophysical Research: Solid Earth*, volume 103, issue B1, pages 749-771.
  https://doi.org/10.1029/97JB02023.
- [**_Wüstefeld, A., Bokelmann, G., Zaroli, C. & Barruol, G. (2008)_**](https://doi.org/10.1016/j.cageo.2007.08.002).
  SplitLab: A shear-wave splitting environment in Matlab.
  *Computers & Geosciences*, volume 34, issue 5, pages 515-528.
  https://doi.org/10.1016/j.cageo.2007.08.002.
  version 1.0.5 **was** available at http://splitting.gm.univ-montp2.fr, version 1.9.0 is available at https://github.com/IPGP/splitlab.



## Funding

The presented research and YF received support from various sources:

- Scholarship of the [Graduate Funding from the German States](https://www.khys.kit.edu/english/graduate_funding.php)
- [DFG grant RI1133/14-1](https://gepris.dfg.de/gepris/projekt/521545943?language=en) within the [DFG Priority Program 2404 DeepDyn](https://www.geo.lmu.de/deepdyn/en/)
