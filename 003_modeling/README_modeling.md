# Instruction: How to perform a modeling run for a single seismological recording station

- created: Michael Grund (ORCID 0000-0001-8759-2018)\
  https://github.com/michaelgrund/sws_tools
- re-written: Yvonne Fröhlich (ORCID 0000-0002-8566-0619)\
  https://github.com/yvonnefroehlich/sws-plotting-and-modeling

_Details_: [Supporting Information](https://academic.oup.com/gji/article/223/3/1525/5893297#supplementary-data) of [Grund & Ritter (2020)](https://doi.org/10.1093/gji/ggaa388.)


## Requirements

- Official _MATLAB_ Toolboxes
  - Deep Learning Toolbox
  - Mapping Toolbox

- _MATLAB Seismic Anisotropy Toolbox_ (MSAT) by [Walker & Wookey (2012)](https://doi.org/10.1016/j.cageo.2012.05.031)
  - Download it from (last access 2022/07/05)
    - https://www1.gly.bris.ac.uk/MSAT/
    - https://github.com/andreww/MSAT
  - Install MSAT
  - Add the whole MSAT package to your _MATLAB_ path
  - If you want to get familiar with the modeling in MSAT and the behavior of different settings
  change to the directory `MSAT/examples/splitting_model/` in which the script `split_model.m` is located.
  This script is described in the publication Walker & Wookey (2012).
  You can play around with different parameters and see how these affect the splitting parameters.


## Forward Calculation

_Please note_: For small stepsizes computation time and structure size increases significantly.

- Precompute splitting parameters for the different model types: **`SWS_modeling_precomp_models_main.m`**

|step_phi / deg|step_dt / s|step_dips / deg|step_dddir / deg|size / GB|publication|
|---|---|---|---|---|---|
|45 |1.00|15|45|0.008|provided TEST_data|
|5  |0.25|5 |5 |4.7  |Ritter, Fröhlich, Sanz Alonso, Grund (2022)|
|5  |0.20|5 |5 |6.6  |Grund & Ritter (2020)|

- All models of all model types are merged in a single nested _MATLAB_ structure with fields:
  - 1 | `modout.phi_eff`: effective or apparent phi values over backazimuth
  - 2 | `modout.dt_eff`: effective or apparent dt values over backazimuth
  - 3 | `modout.mod_paras`: model parameters depending on model type
    - one horizontal layer see `SWS_modeling_precomp_singlelayer.m` _or_
    - two horizontal layers see `SWS_modeling_precomp_twolayers.m` _or_
    - one dippping layer see `SWS_modeling_precomp_dippinglayer.m`, `SWS_modeling_calc_dipping.m`
  - 4 | `modout.type`: string corresponding to the model type
    - one horizontal layer as 'single_layer' _or_
    - two horizontal layers as 'two_layers' _or_
    - one dipping layer as 'dipping'

## Comparison & Visualization

_Please note_: SWS data input is expected to be in standard _SplitLab_ and/or _StackSplit_ output format (`SWS_modeling_read_data.m`).

- The observed SWS data is compared against all pre-computed models of all model types: **`SWS_modeling_calc_misfit.m`**
- The root mean square error (RMSE) between the observed and pre-computed splitting parameters is calculated
- Vizualization of the modeling results (`SWS_modeling_plot_results.m`, `SWS_modeling_plot_stereo_synthetic.m`)
- Take the minimum RMSE or something else metric to get the model which best describes the observed SWS data
