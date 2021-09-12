%{
#include "node.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

NBlock *ProgramBlock;

extern int yylex();
extern int yylineno;
void yyerror(const char *s) { std::printf("Error(line %d): %s\n", yylineno, s); std::exit(1); }
%}

%union {
	NNode *Node;
	NBlock *Block;
	NExpression *Expr;
	NStatement *Stmt;
	NIdentifier *Ident;
	NReaction *Reaction;
	NMoleculeIdentifer *MolIdent;
    NPathwayExpression *PathwayExpr;

	std::vector<NMoleculeIdentifer*> *MolVec;
	std::string *String;
	int Token;
}

%token <Token> T_PROTEIN T_PATHWAY

%token <String> T_ATP T_CO2 T_OXALOACETATE T_ACETYL_COA T_COA
%token <String> T_CITRATE T_ISOCITRATE T_KETO_GLUTARATE
%token <String> T_SUCCINATE T_FUMARATE T_MALATE

%token <String> T_INTEGER T_DOUBLE T_IDENTIFIER

%token <Token> T_LPAREN T_RPAREN T_LBRACE T_RBRACE T_COMMA T_DOT
%token <Token> T_PLUS T_MINUS T_ARROW T_BIARROW
%token <Token> T_EQUAL T_OR T_SEMIC

%type <Ident> ident
%type <Block> program stmts
%type <Stmt> stmt protein_decl pathway_decl
%type <MolVec> mol_expr
%type <MolIdent> mol_ident
%type <Reaction> protein_decl_args
%type <PathwayExpr> pathway_expr pathway_decl_args

%right T_ARROW T_BIARROW
%left T_PLUS
%left T_OR
%start program

%%

program        : stmts { ProgramBlock = $1; }
               ;

stmts          : stmt { $$ = new NBlock(); $$->Statements.push_back($<Stmt>1); }
               | stmts stmt { $1->Statements.push_back($<Stmt>2); } 
               ;

stmt           : protein_decl T_SEMIC
               | pathway_decl T_SEMIC
               ;

protein_decl   : T_PROTEIN ident T_LPAREN protein_decl_args T_RPAREN 
                 { $$ = new NProteinDeclaration(*$2, *$4); }
               ;

protein_decl_args : /* blank */ { $$ = new NReaction(); }
                  | mol_expr T_ARROW mol_expr { $$ = new NReaction(*$1, *$3, false); }
                  | mol_expr T_BIARROW mol_expr { $$ = new NReaction(*$1, *$3, true); }
                  ;

pathway_decl   : T_PATHWAY ident T_EQUAL pathway_decl_args { $$ = new NPathwayDeclaration(*$2, *$4); } 
               ;

pathway_decl_args : /* blank */ { $<Ident>$ = new NIdentifier(); }
                  | pathway_expr 
                  ;

pathway_expr   : ident { $<Ident>$ = new NIdentifier(*$1); delete $1; }
               | pathway_expr T_PLUS pathway_expr  { $$ = new NPathwayExpression(*$1, *$3, $2); }
               | pathway_expr T_OR pathway_expr    { $$ = new NPathwayExpression(*$1, *$3, $2); }
               | pathway_expr T_ARROW pathway_expr { $$ = new NPathwayExpression(*$1, *$3, $2); }
               ;

mol_expr       : mol_ident { $$ = new MoleculeList(); $$->push_back($<MolIdent>1); }
               | mol_expr T_PLUS mol_ident { $1->push_back($<MolIdent>3); }
               ;
	
mol_ident      : T_ATP                  { $$ = new NMoleculeIdentifer(T_ATP, *$1); delete $1; }
               | T_CO2                  { $$ = new NMoleculeIdentifer(T_CO2, *$1); delete $1; } 
               | T_COA                  { $$ = new NMoleculeIdentifer(T_COA, *$1); delete $1; }
               | T_OXALOACETATE         { $$ = new NMoleculeIdentifer(T_OXALOACETATE, *$1); delete $1; }
               | T_ACETYL_COA           { $$ = new NMoleculeIdentifer(T_ACETYL_COA, *$1); delete $1; }
               | T_CITRATE              { $$ = new NMoleculeIdentifer(T_CITRATE, *$1); delete $1; }
               | T_ISOCITRATE           { $$ = new NMoleculeIdentifer(T_ISOCITRATE, *$1); delete $1; }
               | T_KETO_GLUTARATE       { $$ = new NMoleculeIdentifer(T_KETO_GLUTARATE, *$1); delete $1; }
               | T_SUCCINATE            { $$ = new NMoleculeIdentifer(T_SUCCINATE, *$1); delete $1; }
               | T_FUMARATE             { $$ = new NMoleculeIdentifer(T_FUMARATE, *$1); delete $1; }
               | T_MALATE               { $$ = new NMoleculeIdentifer(T_MALATE, *$1); delete $1; }
               ;

ident          : T_IDENTIFIER { $$ = new NIdentifier(*$1); delete $1; } 
               ;
                 
%%