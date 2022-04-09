# lcc-project


# Downloads

# Build
## build lcc-prep
    
    mkdir build
    cd build
    cmake ../
    make
    cd ..


## compile TCA.lpp (or any) and simulate

    ./build/lcc -i . TCA.lpp RNASynthesis.lpp ProteinSynthesis.lpp Ecoli.lpp
    python SimModule.py
    
    
    To see a nice curve:
    ./build/lcc -i . TCA_Reduced_2.lpp
    python SimModule.py
    
## plotting

    ./plot.py



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
