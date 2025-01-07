# 001_stereoplot

 _How to use_: [Header of function `SWS_Analysis_BASICS_stereoplot.m`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/blob/main/001_stereoplot/SWS_Analysis_BASICS_stereoplot.m)



## Overview

- Plot **single-event analysis results** (output of _SplitLab_)
  - _Rotation-correlation_ method ([**_Bowman & Ando 1987_**]( https://doi.org/10.1111/j.1365-246X.1987.tb01367.x.))
  - _Energy minimization_ method ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
  - _Eigenvalue_ method ([**_Silver & Chan 1991_**](https://doi.org/10.1029/91JB00899))
- Plot **multi-event analysis results** (output of _StackSplit_)
  - Stacking of error surfaces (STACK; [**_Wolfe & Silver 1998_**](https://doi.org/10.1029/97JB02023), [**_Restivo & Helffrich 1999_**](https://doi.org/10.1046/j.1365-246x.1999.00845.x))
  - _Splits_ of _simultaneous inversion of multiple waveforms_ (SIMW; [**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899))
- Color-code bars with respect to the fast polarization direction (see Requirements/Colormaps)
- Shade background or backazimuth sector

![figures_SWS_stereo_README_BFO_orange_TEST](https://user-images.githubusercontent.com/94163266/191190219-8570c195-045f-4ad2-9c57-79e140f8e11d.png)
