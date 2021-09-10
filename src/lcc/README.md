# lcc-project






# Downloads
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
