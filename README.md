# Visualizing and Modeling of Shear Wave Splitting

_MATLAB_ functions for visualizing and modeling of shear wave splitting observation:
- Optimized for the output of _SplitLab_ ([**_Wüstefeld et al. 2008_**](https://doi.org/10.1016/j.cageo.2007.08.002)) and [_StackSplit_](https://github.com/michaelgrund/stacksplit) ([**_Grund 2017_**](https://doi.org/10.1016/j.cageo.2017.04.015)).
- The modeling routine is applicable for the energy-minimum method ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899)).
- Strongly modified and extended from https://github.com/michaelgrund/sws_tools ([**_Grund & Ritter 2020_**](https://doi.org/10.1093/gji/ggaa388)).


## Citation
If you make use of this material, please acknowledge the relating publications in which framework these functions were implemented:

- **_Fröhlich, Yvonne, Grund, Michael & Ritter, Joachim R. R. (2022)_**. Laterally and vertically varying seismic anisotropy in the lithosphere-asthenosphere system revealed from SK(K)S splitting at neighboring sites in the Upper Rhine Graben area, Central Europe. in preparation for *Geophysical Journal International*.
- **_Ritter, Joachim R. R., Fröhlich, Yvonne, Sanz Alonso, Yasmin & Grund, Michael (2022)_**. under review by *Journal of Seismology*.
- [**_Grund, Michael & Ritter, Joachim R. R. (2020)_**](https://doi.org/10.1093/gji/ggaa388). Shear-wave splitting beneath Fennoscandia - Evidence for dipping structures and laterally varying multilayer anisotropy. *Geophysical Journal International*, 223, 1525-1547. https://doi.org/10.1093/gji/ggaa388.
- [**_Grund, Michael (2019)_**](https://doi.org/10.5445/IR/1000091425). Exploring geodynamics at different depths with shear wave splitting. *Dissertation*, Karlsruhe Institute of Technology (KIT). https://doi.org/10.5445/IR/1000091425.


Furthermore you can cite the [Zenodo Doi]() given above.


## Content
- [`001_stereoplot`](): Stereoplot representation
  - different shear wave splitting measurement methods
  - singel-event and multi-event anlysis results
  - color-coding of fast polarization direction (see Requirements/Colormaps)
  - shaded background or backzimuth sector
- [`002_xxx`](): General visualization (under development)
- [`003_modeling`](): Modeling of one layer with horizontal and tilted symmetry axis, two layers with horizontal symmetry axes
  - _Forward calculation_: Synthetic splitting parameters for the energy-minimum method in a ray theory reference frame
  - _Observation Comparison_: Minimizing the root mean square error of the splitting parameters
  - _Result visualization_: Backazimuthal variation, model type distribution, model parameter distribution, synthetic stereoplot
  - _Result data_: txt files with model parameter


## Requirements
- **Software**: _MATLAB_ , tested with R2022a, R2021a,b (under Linux and Windows)
- **Forward calculation**: [_MATLAB Seismic Anisotropy Toolbox_ (MSAT)](https://www1.gly.bris.ac.uk/MSAT/) ([**_Walker & Wookey 2012_**](https://doi.org/10.1016/j.cageo.2012.05.031))
- **Data**: Shear wave splitting observations
  - Output txt files (nulls, splits) of _SplitLab_ versions 1.5.0 (original) or 1.2.1 (**_Porritt 2014_**)
  - Output txt files (stack, simw) of _StackSplit_
- **Colormaps** (optional): color-coding of fast polarization direction
  - [MatPlotLib Perceptually Uniform Colormaps](https://de.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps)
    (v2.1.3, last access 2022 June 26)
  - [crameri perceptually uniform scientific colormaps](https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)
    (v1.08, last access 2022 June 25; based on [**_Crameri 2021_**](http://doi.org/10.5281/zenodo.1243862))
  - [cmocean perceptually-uniform colormaps](https://de.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)
    (v2.02, last access 2022 June 18; based on [**_Thyng et al. 2016_**](http://dx.doi.org/10.5670/oceanog.2016.66))


## How to use
- `001_stereoplot`
  - Header of function `SWS_Analysis_BASICS_stereoplot.m`
- `002_xxx`
  - xxx
- `003_modeling`
  - Seperat [README]()
  - Supporting Information of **_Fröhlich et al. (2022)_**
  - [Supporting Information](https://academic.oup.com/gji/article/223/3/1525/5893297#supplementary-data) of **_Grund & Ritter (2020)_**


## Releases
- dev ([main branch]())
- [v1.0]()

For details on the individual releases see the [changelog]().


## Contributing

For bug reports, suggestions or recommendations feel free to open an issue or submit a pull request directly here on GitHub.


## References

[**_Bowman, J. R. & Ando, M. (1987)_**](https://doi.org/10.1111/j.1365-246X.1987.tb01367.x).
Shear-wave splitting in the upper-mantle wedge above the Tonga subduction zone.
*Geophysical Journal International*, volume 88, issue 1, pages 25-41.
https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.

[**_Crameri, F. (2021)_**](http://doi.org/10.5281/zenodo.1243862).
Scientific colour maps. *Zenodo*. http://www.fabiocrameri.ch/colourmaps.php. http://doi.org/10.5281/zenodo.1243862.

[**_Grund, M. (2017)_**](https://doi.org/10.1016/j.cageo.2017.04.015).
StackSplit - a plugin for multi-event shear wave splitting analyses in SplitLab.
*Computers & Geosciences*, volume 105, pages 43-50.
https://doi.org/10.1016/j.cageo.2017.04.015.
versions [1.0](https://doi.org/10.5281/zenodo.464385), 2.0, and [3.0](https://doi.org/10.5281/zenodo.5802051)
available at https://github.com/michaelgrund/stacksplit.

[**_Grund, M. (2019)_**](https://doi.org/10.5445/IR/1000091425).
Exploring geodynamics at different depths with shear wave splitting.
*Dissertation*, Karlsruhe Institute of Technology (KIT). https://doi.org/10.5445/IR/1000091425.

[**_Grund, M. & Ritter, J. R. R. (2020)_**](https://doi.org/10.1093/gji/ggaa388).
Shear-wave splitting beneath Fennoscandia - Evidence for dipping structures and laterally varying multilayer anisotropy.
*Geophysical Journal International*, volume 223, pages 1525-1547.
https://doi.org/10.1093/gji/ggaa388.

**_Porritt, R. W. (2014)_**. SplitLab version 1.2.1.
available at https://robporritt.wordpress.com/software/.

[**_Silver, P. G. & Chan, W. W. (1991)_**](https://doi.org/10.1029/91JB00899).
Shear wave splitting and subcontinental mantle deformation.
*Journal of Geophysical Research*, volume 96, issue B10, pages 16429-16454.
https://doi.org/10.1029/91JB00899.

[**_Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M. & DiMarco, S. F. (2016)_**](http://dx.doi.org/10.5670/oceanog.2016.66).
True colors of oceanography: Guidelines for effective and accurate colormap selection.
*Oceanography*, volume 29, issue 3, pages 9–13.
http://dx.doi.org/10.5670/oceanog.2016.66.

[**_Walker & Wookey (2012)_**](https://doi.org/10.1016/j.cageo.2012.05.031).
MSAT—A new toolkit for the analysis of elastic and seismic anisotropy.
*Computer & Geosciences*, volume 49, pages 81-90.
https://doi.org/10.1016/j.cageo.2012.05.031.
available at https://www1.gly.bris.ac.uk/MSAT/, https://github.com/andreww/MSAT.

[**_Walsh, E., Arnold, R. & Savage, M. K. (2013)_**](https://doi.org/10.1002/jgrb.50386).
Silver and Chan revisited.
*Journal of Geophysical Research: Solid Earth*, volume 118, issue 10, pages 5500-5515.
https://doi.org/10.1002/jgrb.50386.

[**_Wüstefeld, A., Bokelmann, G., Zaroli, C. & Barruol, G. (2008)_**](https://doi.org/10.1016/j.cageo.2007.08.002).
SplitLab: A shear-wave splitting environment in Matlab.
*Computers & Geosciences*, volume 34, issue 5, pages 515-528.
https://doi.org/10.1016/j.cageo.2007.08.002.
version 1.0.5 available at http://splitting.gm.univ-montp2.fr/ and version 1.9.0 available at https://github.com/IPGP/.
