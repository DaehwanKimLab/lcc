# lcc-project


# Downloads

# Build
## build lcc

1. Create a build directory in the project and change to it, and run cmake
   ```sh
   mkdir build
   cd build
   cmake ../
   ```

2. Run make command
   ```sh
   make
   ```

## compile lpp source cfile and run simulation

To run lcc command, activate lcc conda environment and run the command in the project directory.
```sh
conda activate lcc
./build/lcc Models/Ingalls2013/Ingalls2013_Model2_18_NumericalSimulation.lpp
python SimExecutor.py
```

# Build on Windows
Required packages:  
- Visual Studio 2017(or newer version, community version) [Download](https://my.visualstudio.com/Downloads?q=visual%20studio%202017&wt.mc_id=o~msft~vscom~older-downloads)  
- conda (3.17 or higher)

## Installation
- Install Visual Studio 2017
- Create lcc conda envrionment
- Activate lcc conda environment
- Install m2-flex, m2-bison, cmake package in the lcc conda environment  
  `conda install m2-flex m2-bison cmake`

## Build
    cd lcc/src
    mkdir build
    cd build
    cmake ../
    
You can find `lcc.sln` solution file in the build directory. Open the solution file with Visual Studio. Build `lcc` target. You can also build from a command line.

    cmake --build . --target lcc


## Vim syntax highlight
Copy etc/lpp.vim to $HOME/.vim/syntax directory and add below line to .vimrc or create a file to $HOME/.vim/ftdetect/lpp.vim
```
au BufRead,BufNewFile *.lpp set filetype=lpp
```
