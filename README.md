# Visualization and Modeling of Shear Wave Splitting

_MATLAB_ functions for visualization and modeling of shear wave splitting observations:
- Optimized for the output of _SplitLab_ ([**_Wüstefeld et al. 2008_**](https://doi.org/10.1016/j.cageo.2007.08.002)) and [_StackSplit_](https://github.com/michaelgrund/stacksplit) ([**_Grund 2017_**](https://doi.org/10.1016/j.cageo.2017.04.015)).
- The modeling routine is applicable for the _energy minimization method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899)).
- Extended and strongly modified from [sws_tools](https://github.com/michaelgrund/sws_tools) by [Michael Grund](https://github.com/michaelgrund).


## Citation

If you make use of this material, please acknowledge the relating publications in which framework these functions were implemented:

<!---
- **_Fröhlich, Yvonne, Grund, Michael & Ritter, Joachim R. R. (2022)_**. Laterally and vertically varying seismic anisotropy in the lithosphere-asthenosphere system revealed from SK(K)S splitting at neighboring sites in the Upper Rhine Graben area, Central Europe. in preparation for *Geophysical Journal International*.
-->
- **_Ritter, Joachim R. R., Fröhlich, Yvonne, Sanz Alonso, Yasmin & Grund, Michael (2022)_**. Short-scale laterally varying SK(K)S shear wave splitting at BFO, Germany – implications for the determination of anisotropic structures. accepted by *Journal of Seismology*. DOI: 10.1007/s10950-022-10112-w.
- [**_Grund, Michael & Ritter, Joachim R. R. (2020)_**](https://doi.org/10.1093/gji/ggaa388). Shear-wave splitting beneath Fennoscandia - Evidence for dipping structures and laterally varying multilayer anisotropy. *Geophysical Journal International*, 223, 1525-1547. https://doi.org/10.1093/gji/ggaa388.
- [**_Grund, Michael (2019)_**](https://doi.org/10.5445/IR/1000091425). Exploring geodynamics at different depths with shear wave splitting. *Dissertation*, Karlsruhe Institute of Technology (KIT). https://doi.org/10.5445/IR/1000091425.

<!---
Furthermore you can cite the [Zenodo Doi]() given above.
-->


## Content

### **[001_stereoplot](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/001_stereoplot)**

 _How to use_: Header of function [`SWS_Analysis_BASICS_stereoplot.m`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/blob/main/001_stereoplot/SWS_Analysis_BASICS_stereoplot.m)

- Plot **single-event analysis results** (output of _SplitLab_)
  - _Rotation-correlation method_ ([**_Bowman & Ando 1987_**]( https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.))
  - _Energy minimization method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
  - _Eigenvalue method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
- Plot **multi-event analysis results** (output of _StackSplit_)
  - Stacking of error surfaces ([**_Wolfe & Silver 1998_**](https://doi.org/10.1029/97JB02023), [**_Restivo & Helffrich 1999_**](https://doi.org/10.1046/j.1365-246x.1999.00845.x))
  - _Splits_ of _simultaneous inversion of multiple waveforms_ ([**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899))
- Color-code bars with respect to the fast polarization direction (see Requirements/Colormaps)
- Shade background or backazimuth sector

_Example figures produced with the provided [TEST_data_stereoplot](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/001_stereoplot/TEST_data_stereoplot)_:

![figures_SWS_stereo_README_BFO_orange_TEST](https://user-images.githubusercontent.com/94163266/191190219-8570c195-045f-4ad2-9c57-79e140f8e11d.png)

### **[002_visualization](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization)**

- under development

<!---
_How to use_:

- xxx
- xxx

_Example figures produced with the provided [Test_data_visualization]()_:
-->

### **[003_modeling](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling)**

_How to use_: [README](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/blob/main/003_modeling/README.md); header of function [`SWS_modeling_calc_misfit.m`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/blob/main/003_modeling/SWS_modeling_calc_misfit.m)\
_Details_: [Supporting Information](https://academic.oup.com/gji/article/223/3/1525/5893297#supplementary-data) of [**_Grund & Ritter (2020)_**](https://doi.org/10.1093/gji/ggaa388)

- **Model types**
  - Transverse isotropy
  - Structural anisotropy: one horizontal or dipping layer, two horizontal layers
- **Forward calculation**
  - Synthetic splitting parameters based on the _energy minimization method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
  - Ray theory reference frame
- **Observation** (output of _SplitLab_ and _StackSplit_)
  - Single-event analysis: _energy minimization method_ ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
  - Multi-event analysis:
    stacking of error surfaces ([**_Wolfe & Silver 1998_**](https://doi.org/10.1029/97JB02023), [**_Restivo & Helffrich 1999_**](https://doi.org/10.1046/j.1365-246x.1999.00845.x)),
    _simultaneous inversion of multiple waveforms_ ([**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899))
- **Comparison of forward calculation and observation**
  - Calculation and minimizing the root mean square error regarding the splitting parameters
  - Joint fitting of fast polarization direction and delay time, separate fitting of the fast polarization direction
- **Result visualization**
  - Backazimuthal variation of the splitting parameters (forward calculation and observation)
  - Model type distribution (histogram)
  - Model parameter distribution (scatter plot)
  - Stereoplot of synthetic splitting parameters (polar plot)
- **Result data**
  - Model parameters for each model type as separate *.txt files

_Example figures produced with the provided [TEST_data_modeling](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/003_modeling/TEST_data_modeling)_:

![figure_SWS_modeling_README_BFO_orange_TEST](https://user-images.githubusercontent.com/94163266/183926434-f510331c-0ded-4fb9-9867-b30727568432.png)


## Requirements

_Tested with_: R2022a, R2021a,b under Linux and Windows

- **_MATLAB_**: Forward calculation
  - Deep Learning Toolbox
  - Mapping Toolbox
  - [_MATLAB Seismic Anisotropy Toolkit_ (MSAT)](https://www1.gly.bris.ac.uk/MSAT/) ([**_Walker & Wookey 2012_**](https://doi.org/10.1016/j.cageo.2012.05.031))
- **Data**: Shear wave splitting observations
  - Output *.txt files (nulls, splits) of _SplitLab_ versions 1.5.0 (original) or 1.2.1 (**_Porritt 2014_**)
  - Output *.mat structure and *.txt files (stack, simw) of _StackSplit_
- **Colormaps** (optional): Color-coding of to the fast polarization direction and the root mean square error
  - [MatPlotLib Perceptually Uniform Colormaps](https://de.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps)\
    version v2.1.3, MATLAB File Exchange, last access 2022 June 26
  - [crameri perceptually uniform scientific colormaps](https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)\
    version v1.08, MATLAB File Exchange, last access 2022 June 25; based on [**_Crameri (2021)_**](http://doi.org/10.5281/zenodo.1243862)
  - [cmocean perceptually-uniform colormaps](https://de.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)\
    version v2.02, MATLAB File Exchange, last access 2022 June 18; based on [**_Thyng et al. (2016)_**](http://dx.doi.org/10.5670/oceanog.2016.66)


## Releases

- dev ([main branch](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main))
<!---
- [v1.0]()
-->

For details of the individual releases as well as for changes and differences compared to [sws_tools](https://github.com/michaelgrund/sws_tools) by [Michael Grund](https://github.com/michaelgrund) see the [changelog]().


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

[**_Restivo, A. & Helffrich, G. (1999)_**](https://doi.org/10.1046/j.1365-246x.1999.00845.x).
Teleseismic shear wave splitting measurements in noisy environments.
*Geophysical Journal International*, volume 137, pages 821-830.
https://doi.org/10.1046/j.1365-246x.1999.00845.x.

[**_Roy, C., Winter, A., Ritter, J. R. R. & Schweitzer, J. (2017)_**](https://doi.org/10.1093/gji/ggw470).
On the improvement of SKS splitting measurements by the Simultaneous Inversion of Multiple Waveforms (SIMW).
*Geophysical Journal International*, volume 208, pages 1508-1523.
https://doi.org/10.1093/gji/ggw470.

[**_Silver, P. G. & Chan, W. W. (1991)_**](https://doi.org/10.1029/91JB00899).
Shear wave splitting and subcontinental mantle deformation.
*Journal of Geophysical Research*, volume 96, issue B10, pages 16429-16454.
https://doi.org/10.1029/91JB00899.

[**_Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M. & DiMarco, S. F. (2016)_**](http://dx.doi.org/10.5670/oceanog.2016.66).
True colors of oceanography: Guidelines for effective and accurate colormap selection.
*Oceanography*, volume 29, issue 3, pages 9-13.
http://dx.doi.org/10.5670/oceanog.2016.66.

[**_Walker, A. M. & Wookey, J. (2012)_**](https://doi.org/10.1016/j.cageo.2012.05.031).
MSAT — A new toolkit for the analysis of elastic and seismic anisotropy.
*Computer & Geosciences*, volume 49, pages 81-90.
https://doi.org/10.1016/j.cageo.2012.05.031.
available at https://www1.gly.bris.ac.uk/MSAT/, https://github.com/andreww/MSAT.

[**_Wolfe, C. J. & Silver, P. G. (1998)_**](https://doi.org/10.1029/97JB02023).
Seismic anisotropy of oceanic upper mantle: Shear wave splitting methodologies and observations.
*Journal of Geophysical Research: Solid Earth*, volume 103, issue B1, pages 749-771.
https://doi.org/10.1029/97JB02023.

[**_Wüstefeld, A., Bokelmann, G., Zaroli, C. & Barruol, G. (2008)_**](https://doi.org/10.1016/j.cageo.2007.08.002).
SplitLab: A shear-wave splitting environment in Matlab.
*Computers & Geosciences*, volume 34, issue 5, pages 515-528.
https://doi.org/10.1016/j.cageo.2007.08.002.
version 1.0.5 available at http://splitting.gm.univ-montp2.fr/ and version 1.9.0 available at https://github.com/IPGP/.
