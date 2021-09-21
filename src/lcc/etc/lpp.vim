" Vim Syntax file
" Language: LPP
"
" Place this file to $HOME/.vim/syntax/lpp.vim
"
" Add below line to .vimrc or create a file to $HOME/.vim/ftdetect/lpp.vim
" au BufRead,BufNewFile *.lpp set filetype=lpp
"

syn keyword lppTodo contained TODO FIXME XXX NOTE
syn keyword lppTopKeyword organism experiment protein pathway
syn keyword lppTopKeyword description reaction_id reaction
syn keyword lppTopKeyword property using module 
syn keyword lppMolecule oxaloacetate citrate
syn keyword lppMolecule isocitrate
syn match lppMolecule display "keto-glutarate"
syn match lppMolecule display "acetyl-CoA"
syn match lppOperators display "+"
syn match lppOperators display "-->"
syn match lppOperators display "|"

"<!-- LPP -->

syn region lppString start='"' skip='\\"' end='"'

syn region lppComment start="//" end="$" contains=lppTodo
syn region lppBlockComment start="/\*" end="\*/" contains=lppTodo


hi def link lppTopKeyword Type
hi def link lppMolecule Special
hi def link lppComment Comment
hi def link lppBlockComment Comment
hi def link lppString Constant
hi def link lppOperators Operator
hi def link lppTodo Todo
