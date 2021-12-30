# lcc-project


# Downloads

# Build
## build lcc-prep
    
    mkdir build
    cd build
    cmake ../
    make


## compile TCA.lpp (or any) and simulate

    ./build/lcc -i . TCA.lpp RNASynthesis.lpp ProteinSynthesis.lpp Ecoli.lpp

    python SimModule.py

## plotting

    ./plot.py



# Build on Windows
Required packges:  
- Visual Studio 2017(or newer version, community version) [Download](https://my.visualstudio.com/Downloads?q=visual%20studio%202017&wt.mc_id=o~msft~vscom~older-downloads)  
- win_flex_bison package [Download](https://sourceforge.net/projects/winflexbison/files/latest/download)  
- cmake (>=3.17) [Download](https://github.com/Kitware/CMake/releases/download/v3.22.1/cmake-3.22.1-windows-x86_64.msi)  
- python

## Installation
- Install Visual Studio 2017
- Download win_flex_bison package. Uncompress it. Make sure there are no spaces in the directory name. You can find the win_bison.exe and win_flex.exe files in the directory. Rename the files to bison.exe and flex.exe. Add the directory name to the PATH environment. Open Command Prompt and run `bison.exe -h`, `flex.exe -h`
- Install cmake package
- Install python package. I recommend using conda environment. (use same environment file for lcc)

## Build
    cd lcc/src
    mkdir build
    cd build
    cmake ../

You can find `lcc.sln` solution file in the build directory. Open the solution file with Visual Studio. Build `lcc-prep` target.


## How to change PATH environment variable
`Start Menu` --> `Settings` --> Search 'Edit the system environment variables' on the `Find a setting` box. --> Click `Environment Variables` --> Select `Path` and click `Edit` button. --> Click `New` button and put the directory name --> Ok...

## Vim syntax highlight
Copy etc/lpp.vim to $HOME/.vim/syntax directory and add below line to .vimrc or create a file to $HOME/.vim/ftdetect/lpp.vim
```
au BufRead,BufNewFile *.lpp set filetype=lpp
```
