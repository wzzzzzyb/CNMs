# Utilizing Causal Network Markers to Identify Tipping Points ahead of Critical Transition
*Early-warning signal, Identification of clinical disease, Causal network marker*

[![DOI](https://img.shields.io/badge/DOI-arXiv.2412.16235-yellow)](https://doi.org/10.48550/arXiv.2412.16235)
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
The MATLAB script 📊`Five_genetic_validation.m`, contained within the 📁`Validations on Benchmark models` directory, generates the results for five-genetic network data by Langevin simulation and gives the CNMs results.

### How to Run
Run directly.

## 2. Validations on Real-world Epileptic iEEG Datasets
### Description
The MATLAB script 📊`iEEG_data_validation.m`, contained within the 📁`Validations on iEEG Datasets` directory, is responsible for generating all subplot visualizations of the iEEG results that appear in the 📰`Supporting Information` of our publication.

### How to Run
1. Read the 📄`iEEG Dataset Download Instruction.pdf` to download the necessary datasets;
2. Run the 📊`iEEG_data_validation.m` directly to generate the MATLAB visualization plots.

### Customization​
To analyze different iEEG datasets, modify the following parameters. Please note that the referenced datasets must be pre-downloaded.
```bibtex
ID = 1; Sz = 4;         % Select your interested dataset
```

## Citation
The BibTEX format is as given:
```bibtex
@article{bian2024cnms,
  title={Utilizing Causal Network Markers to Identify Tipping Points ahead of Critical Transition},
  author={Bian, Shirui and Wang, Zezhou and Leng, Siyang and Lin, Wei and Shi, Jifan},
  journal={arXiv preprint arXiv:2412.16235},
  year={2024}
}
```

## License
Open-source under MIT License.

## Additional Acknowledgments
We thank **Liuqian Guo** from *Dalian University of Technology* for her contribution to drawing the elements for Figure 2 (d) and (e) : )
