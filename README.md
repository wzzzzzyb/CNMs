# Utilizing Causal Network Markers to Identify Tipping Points ahead of Critical Transition
*Early-warnin signal, Identification of clinical disease, Causal network marker*

[![DOI](https://img.shields.io/badge/DOI-advs.202415732-yellow)](https://doi.org/10.1002/advs.202415732)c
![MATLAB](https://img.shields.io/badge/MATLAB-R2023a-%23FF9500?logo=mathworks&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)
![Last Commit](https://img.shields.io/github/last-commit/wzzzzzyb/CNMs?color=blue)

## Overview
This repository contains two complementary MATLAB modules for implementing CNMs:
1. 🧬 ​​Validations on Benchmark Models;
2. 🧠 Validations on Real-world epileptic iEEG Datasets.

The code requires MATLAB R2023a (or other compatible version) with the Parallel Computing Toolbox

## 1. Validations on Benchmark Models
### Description
The MATLAB script 📊`Five_genetic_validation.m`, contained within the 📁`Validations on Benchmark models` directory, generates the results for five-genetic network data by Langevin simulation and gives the CNMs results. The implemented code structure is generalizable to both ecological networks and Turing networks.

### How to Run
Run 📊`Five_genetic_validation.m` directly.

### Customization​
The dynamic simulation parameters are fully customizable. Users may modify:
```matlab
dt = 0.008;                     % Time steps
steps = 1000;                   % Simulation steps
par_num = 2000;                 % Parallel trajectories
D = 1e-5;                       % Noise strength
P = [0.01:0.01:0.1,0.2:0.1:1];  % Selected parameters
```

## 2. Validations on Real-world Epileptic iEEG Datasets
### Description
The MATLAB script 📊`iEEG_data_validation.m`, contained within the 📁`Validations on iEEG Datasets` directory, is responsible for generating all subplot visualizations of the iEEG results that appear in the 📰`Supporting Information` of our publication.

### How to Run
1. Read the 📄`iEEG Dataset Download Instruction.pdf` to download the necessary datasets;
2. Run the 📊`iEEG_data_validation.m` directly to generate the MATLAB visualization plots.

### Customization​
To analyze different iEEG datasets, modify the following parameters: 
```matlab
ID = 1; Sz = 4;         % Select your interested dataset
```
(Please note that the referenced datasets must be pre-downloaded.)

## Citation
The APA format:
``` markdown
​​Bian, S., Wang, Z., Leng, S., Lin, W., & Shi, J.​​ (2025). Utilizing causal network markers to identify tipping points ahead of critical transition. Advanced Science, Article e15732. https://doi.org/10.1002/advs.202415732
```
The BibTEX format:
```bibtex
@article{Bian2025,
  author = {Bian, Shirui and Wang, Zezhou and Leng, Siyang and Lin, Wei and Shi, Jifan},
  title = {Utilizing Causal Network Markers to Identify Tipping Points ahead of Critical Transition},
  journal = {Advanced Science},
  pages = {e15732},
  year = {2025},
  doi = {https://doi.org/10.1002/advs.202415732},
}
```

## License
Open-source under MIT License.

## Additional Acknowledgments
We thank **Liuqian Guo** from *Dalian University of Technology* for her contribution to drawing the elements for Figure 2 (d) and (e) : )
<div align="center">
  <img src="https://github.com/wzzzzzyb/CNMs/raw/main/Figure/Figure2%20(d-e).png" width="95%">
</div>
