# lcc - Life Compiler Collection
Required programs: **Python3**\
Required Python packages: **NumPy**, **TensorFlow**\
Required dataset formats: **txt**, **tsv**, **csv**

This repository contains a **lcc**, or **Life Compiler Collection**, 
which compiles a code for **Function of Life** of the [LDP team](https://kim-lab.org/). 
Function of Life is a mathematical function that describes cellular processes,
using a series of matrix operations on GPU.

## Setup

For the best experience, we recommend using [Conda](https://docs.conda.io/projects/conda/en/latest/#) 
as a package, dependency and environment management system 
that is compatible with Windows, macOS, and Linux. 

### Setup Python Environment

    conda env create -f environment.yml

    conda activate lcc
    conda deactivate


#### GPU Support

	conda env create -f environment-gpu.yml

	conda activate lcc-gpu
	conda deactivate


GPU environment requires NVIDIA GPU Driver 450.x or higher

## Input

## Output
