# Modelling of Shear Wave Splitting Observations

_MATLAB_ functions for visualization and modelling of shear wave splitting observation:
- Optimized for the output of _SplitLab_ ([**_Wüstelfeld et al. 2008_**](https://doi.org/10.1016/j.cageo.2007.08.002)) and [_StackSplit_](https://github.com/michaelgrund/stacksplit) ([**_Grund 2017_**](https://doi.org/10.1016/j.cageo.2017.04.015)).
- The modelling procedure is applicable for the energy-minimum method ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899)).
- Strongly modified and extended from https://github.com/michaelgrund/sws_tools ([**_Grund & Ritter 2020_**](https://doi.org/10.1093/gji/ggaa388)).


## Citation
If you make use of this material, please acknowledge the relating publications in which framework these functions were implemented:

- **_Fröhlich, Yvonne, Grund, Michael & Ritter, Joachim R. R. (2022)_**. Laterally and vertically varying seismic anisotropy in the lithosphere-asthenosphere system revealed from SK(K)S splitting at neighboring sites in the Upper Rhine Graben area, Central Europe. in preperation for *Geophysical Journal International*.
- **_Ritter, Joachim R. R., Fröhlich, Yvonne, Sanz Alonso, Yasmin & Grund, Michael (2022)_**. under review by *Journal of Seismology*.
- [**_Grund, Michael & Ritter, Joachim R. R. (2020)_**](https://doi.org/10.1093/gji/ggaa388). Shear-wave splitting beneath Fennoscandia - Evidence for dipping structures and laterally varying multilayer anisotropy. *Geophysical Journal International*, 223, 1525–1547. https://doi.org/10.1093/gji/ggaa388.
- [**_Grund, M. (2019)_**](https://doi.org/10.5445/IR/1000091425). Exploring geodynamics at different depths with shear wave splitting. *Dissertation*, Karlsruhe Institute of Technology (KIT). https://doi.org/10.5445/IR/1000091425.


Furthermore you can cite the [Zenodo Doi]() given above.


## Content
- folder `001_stereoplots`: Stereoplot representation of shear wave splitting observations
- folder `002_modelling`: Modelling of shear wave splitting observations
  - Forward calculation of synthetic splitting parameters for the energy-minimum method in a ray theory reference frame
  - Comparison with the observations by minimizing the root mean square error
  - Testing for different scenarios: one layer with horizontal and tilted symmetry axis, two layers with horizontal symmetry axes
  - Visualization: backazimuthal variation, model parameter distribution, model type distribution 


## Requirements
- **Software**: _MATLAB_, tested with R2022a, R2021a,b (under Linux and Windows)
- **Forward calculation**: _MATLAB Seismic Anisotropy Toolbox_ (MSAT) ([**_Walker & Wookey 2012_**](https://doi.org/10.1016/j.cageo.2012.05.031))
- **Data**: Shear wave splitting observations made with the energy-minimum method ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
- **Format**: Output (txt files) of _SplitLab_ ([**_Wüstefeld et al., 2018_**](https://doi.org/10.1016/j.cageo.2007.08.002)) version 1.2.1 (**_Poritt 2014_**) and of _StackSplit_ ([**_Grund 2017_**](https://doi.org/10.1016/j.cageo.2017.04.015))


## How to use
-
-
-

see also [Supporting Information](https://academic.oup.com/gji/article/223/3/1525/5893297#supplementary-data) of **_Grund & Ritter (2020)_**.


## Releases
- dev ([main branch]())
- [v1.0]()

For details on the individual releases see the [changelog]().


## Contributing
 
For bug reports, suggestions or recommondations feel free to open an issue or submitt a pull request directly on GitHub.


## References

[**_Bowman, J. R. & Ando, M. (1987)_**](https://doi.org/10.1111/j.1365-246X.1987.tb01367.x).
Shear-wave splitting in the upper-mantle wedge above the Tonga subduction zone.
*Geophysical Journal International*, volume 88, issue 1, pages 25-41.
https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.

[**_Grund, M. & Ritter, J.R.R. (2020)_**](https://doi.org/10.1093/gji/ggaa388).
Shear-wave splitting beneath Fennoscandia - Evidence for dipping structures and laterally varying multilayer anisotropy.
*Geophysical Journal International*, volume 223, pages 1525–1547.
https://doi.org/10.1093/gji/ggaa388.

[**_Grund, M. (2019)_**](https://doi.org/10.5445/IR/1000091425).
Exploring geodynamics at different depths with shear wave splitting.
*Dissertation*, Karlsruhe Institute of Technology (KIT). https://doi.org/10.5445/IR/1000091425.
 
[**_Grund, M. (2017)_**](https://doi.org/10.1016/j.cageo.2017.04.015).
StackSplit - a plugin for multi-event shear wave splitting analyses in SplitLab.
*Computers & Geosciences*, volume 105, pages 43-50.
https://doi.org/10.1016/j.cageo.2017.04.015.
versions [1.0](https://doi.org/10.5281/zenodo.464385), 2.0, and [3.0](https://doi.org/10.5281/zenodo.5802051)
available at https://github.com/michaelgrund/stacksplit.

**_Porritt, R. W. (2014)_**. SplitLab version 1.2.1.
available at https://robporritt.wordpress.com/software/.

[**_Silver, P. G. & Chan, W. W. (1991)_**](https://doi.org/10.1029/91JB00899).
Shear wave splitting and subcontinental mantle deformation.
*Journal of Geophysical Research*, volume 96, issue B10, pages 16429-16454.
https://doi.org/10.1029/91JB00899.

[**_Walker & Wookey (2012)_**](https://doi.org/10.1016/j.cageo.2012.05.031).
MSAT—A new toolkit for the analysis of elastic and seismic anisotropy.
*Computer & Geosciences*, volume 49, pages 81-90.
https://doi.org/10.1016/j.cageo.2012.05.031.
available at https://github.com/andreww/MSAT.

[**_Walsh, E., Arnold, R. & Savage, M. K. (2013)_**](https://doi.org/10.1002/jgrb.50386).
Silver and Chan revisited.
*Journal of Geophysical Research: Solid Earth*, volume 118, issue 10, pages 5500-5515.
https://doi.org/10.1002/jgrb.50386.

[**_Wüstefeld, A., Bokelmann, G., Zaroli, C. & Barruol, G.  (2008)_**](https://doi.org/10.1016/j.cageo.2007.08.002).
SplitLab: A shear-wave splitting environment in Matlab.
*Computers & Geosciences*, volume 34, issue 5, pages 515-528.
https://doi.org/10.1016/j.cageo.2007.08.002.
version 1.0.5 available at http://splitting.gm.univ-montp2.fr/ and version 1.9.0 available at https://github.com/IPGP/.
