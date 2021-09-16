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



## make token list from the metabolites.tsv

	./generate_mol_keyword.py -m ../../data/metabolites.tsv -o etc/mol_list
	
    # make vim syntax file
	./generate_mol_keyword.py -k etc/mol_list -v etc/lpp.vim -o ~/.vim/syntax/lpp.vim

	./generate_mol_keyword.py -k etc/mol_list -l lpp.l.template -o lpp.l
	./generate_mol_keyword.py -k etc/mol_list -y lpp.y.template -o lpp.y
	make


