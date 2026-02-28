# 002_visualization

 _How to use_: See [`load_models.py`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization//02_explore_syn_split_para/load_models.py) and [`explore_syn_split_para.py`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization//02_explore_syn_split_para/explore_syn_split_para.py)

_Test data_: Use files from [000_test data](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/000_test_data)


## Overview

- **[01_load_radar4kit_data](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization/01_load_radar4kit_data)**: Scripts to load datasets provided via RADAR4KIT in MATLAB, Python, and R
  - Upper Rhine Graben Area (URG): https://dx.doi.org/10.35097/685; related to [**_Fröhlich et al. (2024)_**](https://doi.org/10.1093/gji/ggae245)
  - Blackforest Observatory (BFO): https://dx.doi.org/10.35097/684; related to [**_Ritter et al. (2022)_**](https://doi.org/10.1007/s10950-022-10112-w)
- **[02_explore_syn_split_para](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization/02_explore_syn_split_para)**: Explore forward calculated splitting parameters in Python via PyGMT
  - [`separate_modout_struct.m`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization/02_explore_syn_split_para/separate_modout_struct.m): Split MATLAB struct into separate structs for the different model types
  - [`load_models.py`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization//02_explore_syn_split_para/load_models.py): Function to load forward calculated splitting parameters for structural anisotropy models
  - [`explore_syn_split_para.py`](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/tree/main/002_visualization//02_explore_syn_split_para/explore_syn_split_para.py): Plot forward calculated splitting parameters for structural anisotropy models

![](https://github.com/yvonnefroehlich/sws-visualization-and-modeling/raw/main/_images/002_visualization_readme_image.png)
