# lcc-project


# Downloads

# Build
## build lcc-prep
    
	mkdir build
	cd build
	cmake ../
	make


## compile TCA.lpp and simulate

        ./build/lcc-prep -i . TCA.lpp --simmodule Sim.py

## plotting

	./plot.py TCA_simout.tsv


## Vim syntax highlight
Copy etc/lpp.vim to $HOME/.vim/syntax directory and add below line to .vimrc or create a file to $HOME/.vim/ftdetect/lpp.vim
```
au BufRead,BufNewFile *.lpp set filetype=lpp
```


# Deprecated
## llvm
https://llvm.org/docs/GettingStarted.html  

Installing llvm development environment:  

    git clone https://github.com/llvm/llvm-project.git
    cd llvm-project
    mkdir build; cd build
    #cmake -G "Unix Makefiles" ../llvm -DCMAKE_INSTALL_PREFIX="/opt/llvm"
    cmake -G "Unix Makefiles" ../llvm -DLLVM_ENABLE_PROJECTS=clang -DCMAKE_INSTALL_PREFIX="/opt/llvm"
    make -j5
    make install 
    (or) DESTDIR=/home/somewhere make install

