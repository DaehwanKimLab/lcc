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

## Program Design

###**lcc.py**
lcc.py parse the genome using compiler database and generates a whole cell simulation code. 

#### Input
- Genome sequence file(s) (.fasta)

#### Output
- cell.py: Whole cell simulation execution code, a.k.a. **Function of Life**
- CompilerData save files (.npy)

###**cell.py**
cell.py 

    # Instantiate all data components.
    Cst = FConstant()
    Env = FEnvironment()
    Cel = FCellState()
    
    # Instantiate all reaction components.
    Exe = FRateGaugeModelOnly()
    Bch = FBiochemicalReactionRateFunction()
    Pol = FPolymerizationRateFunction()
    
    # Instantiate cell process objects.
    Replication = FReplication(Bch, Cel, Cst, Env, Exe, Pol)
    
    # Generate a dictionary of cell process object names
    Dict_CellProcesses = dict()
    Dict_CellProcesses['Replication'] = Replication
    
    # Instantiate simulation object.
    Sim = FSimulation(Bch, Cel, Cst, Env, Exe, Pol, Dict_CellProcesses)
        
    # Declare temporary parameters
    Cel.Vol = tf.constant([7e-16])
    
    # Run simulation.
    Sim.Initialize()
    Sim.Run()

#### Input
- CompilerData save files (.npy)
- 

#### Output
- Molecule count readout
- 
