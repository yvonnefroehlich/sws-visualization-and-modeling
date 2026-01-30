# 003_modeling

_See also_: [Supporting Information](https://academic.oup.com/gji/article/223/3/1525/5893297#supplementary-data) of [**_Grund, Ritter (2020)_**](https://doi.org/10.1093/gji/ggaa388)



## Overview

- **Model types**
  - Transverse isotropy
  - Structural anisotropy: one horizontal layer (H1), one dipping layer (T1), two horizontal layers (H2)
- **Forward calculation**
  - Synthetic splitting parameters based on the _energy minimization_ method ([**_Silver, Chan 1991_**](https://doi.org/10.1029/91JB00899))
  - Ray theory reference frame
- **Observation** (output of _SplitLab_ and _StackSplit_)
  - Single-event analysis: _energy minimization_ method ([**_Silver, Chan 1991_**](https://doi.org/10.1029/91JB00899))
  - Multi-event analysis:
    stacking of error surfaces (STACK; [**_Wolfe, Silver 1998_**](https://doi.org/10.1029/97JB02023), [**_Restivo, Helffrich 1999_**](https://doi.org/10.1046/j.1365-246x.1999.00845.x)),
    _simultaneous inversion of multiple waveforms_ (SIMW; [**_Roy et al. 2017_**](https://doi.org/10.1029/91JB00899))
- **Comparison of forward calculation and observation**
  - Calculation and minimizing the root mean square error regarding the splitting parameters
  - Joint fitting of fast polarization direction and delay time, separate fitting of the fast polarization direction
- **Result visualization**
  - Backazimuthal variation of the splitting parameters (forward calculation and observation)
  - Model type distribution (bar plot)
  - Model parameter distribution (scatter plot)
  - Stereoplot of synthetic splitting parameters (polar plot)
- **Result data**
  - Model parameters for each model type as separate *.txt files

![](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/raw/main/_images/003_modeling_readme_image.png)



## Instruction: How to perform a modeling run for a single seismological recording station


### Requirements

- Official _MATLAB_ Toolboxes
  - Deep Learning Toolbox
  - Mapping Toolbox

- _MATLAB Seismic Anisotropy Toolkit_ (MSAT) by [**_Walker, Wookey (2012)_**](https://doi.org/10.1016/j.cageo.2012.05.031)
  - Download MSAT from https://www1.gly.bris.ac.uk/MSAT/ or https://github.com/andreww/MSAT (last access 2022/07/05).
  - Add the whole MSAT package to your _MATLAB_ path.
  - If you want to get familiar with the modeling in MSAT and the behavior of different settings change to the directory
    `MSAT/examples/splitting_model/` in which the script `split_model.m` is located. This script is described in the
    publication mentioned above. You can play around with different parameters and see how these affect the splitting
    parameters.


### Forward Calculation

_Please note_: For small step sizes, computation time and structure size increase significantly

- Pre-compute splitting parameters for the different model types: **`SWS_modeling_precomp_models_main.m`**

| phi / deg | dt / s | dip angle / deg | down-dip direction / deg | thickness / km | size / GB | publication |
|---|---|---|---|---|---|---|
| [-45:45:90] | [1:1:4]       | [15:15:75] | [0:45:315] | [250:250:500] | 0.008 | [sws_modout_domper8s.mat](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/000_test_data) |
| [-85:5:90]  | [0.2:0.2:4]   | [5:5:75]   | [0:5:355]  | [5:5:250]     | 5.9   | [Fröhlich et al. (2024)](https://doi.org/10.1093/gji/ggae245) |
| [-85:5:90]  | [0.25:0.25:4] | [5:5:75]   | [0:5:355]  | [5:5:250]     | 4.2   | [Ritter et al. (2022)](https://doi.org/10.1007/s10950-022-10112-w) |
| [-90:5:90]  | [0.2:0.2:4]   | [5:5:75]   | [0:5:360]  | [5:25:500]    | 5.7   | [Grund & Ritter (2020)](https://doi.org/10.1093/gji/ggaa388) |

- All models of all model types are merged in a single nested _MATLAB_ structure with fields:
  - 1 | `modout.phi_eff`: effective or apparent phi values over backazimuth
  - 2 | `modout.dt_eff`: effective or apparent dt values over backazimuth
  - 3 | `modout.mod_paras`: model parameters depending on the model type
    - One horizontal layer see `SWS_modeling_precomp_singlelayer.m` _or_
    - Two horizontal layers see `SWS_modeling_precomp_twolayers.m` _or_
    - One dippping layer see `SWS_modeling_precomp_dippinglayer.m`, `SWS_modeling_calc_dipping.m`
  - 4 | `modout.type`: string corresponding to the model type
    - One horizontal layer as 'single_layer' _or_
    - Two horizontal layers as 'two_layers' _or_
    - One dipping layer as 'dipping'


### Comparison & Visualization

_Please note_: SWS data input is expected to be in standard _SplitLab_ and _StackSplit_ output formats (`SWS_modeling_read_data.m`)

- The observed SWS data is compared against all pre-computed models of all model types: **`SWS_modeling_calc_misfit.m`**
- The root mean square error (RMSE) between the observed and precomputed splitting parameters is calculated
- Visualization of the modeling results (`SWS_modeling_plot_results.m`, `SWS_modeling_plot_stereo_synthetic.m`)
- Take the minimum RMSE or something else metric to get the model which best describes the observed SWS data
