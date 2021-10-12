%{
#include "node.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <memory>

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
	NMoleculeReaction *MolReaction;
	NReaction *Reaction;
	NMoleculeIdentifier *MolIdent;
	NPathwayExpression *PathwayExpr;

	std::vector<std::shared_ptr<NMoleculeIdentifier>> *MolVec;
	std::vector<std::shared_ptr<NIdentifier>> *IdentVec;
	std::string *String;
	int Token;
}

%token <Token> T_PROTEIN T_PATHWAY T_EXPERIMENT T_ORGANISM
%token <Token> T_DESCRIPTION T_REACTION T_REACTION_ID
%token <Token> T_PROPERTY T_USING T_MODULE
%token <Token> T_COFACTOR T_DOMAIN T_STEP T_SEQUENCE

%token <String> T_STRING_LITERAL

%token <String> T_ATP T_CO2 T_OXALOACETATE T_ACETYL_COA T_COA
%token <String> T_CITRATE T_ISOCITRATE T_KETO_GLUTARATE
%token <String> T_SUCCINATE T_FUMARATE T_MALATE

%token <String> T_NUMBER T_INTEGER T_DOUBLE T_IDENTIFIER

%token <Token> T_LPAREN T_RPAREN T_LBRACE T_RBRACE T_COMMA T_DOT
%token <Token> T_PLUS T_MINUS T_ARROW T_BIARROW
%token <Token> T_EQUAL T_OR T_SEMIC

%type <Ident> ident gen_ident
%type <Block> program stmts block
%type <Block> pathway_block pathway_stmts
%type <Block> protein_block protein_stmts
%type <Stmt> stmt protein_decl pathway_decl
%type <MolVec> mol_expr
%type <MolIdent> mol_ident
%type <Reaction> protein_decl_args
%type <Reaction> gen_reaction_decl_args protein_step_decl_args
%type <PathwayExpr> pathway_expr pathway_decl_args
%type <Stmt> organism_decl experiment_decl pathway_stmt protein_stmt
%type <Stmt> pathway_description_stmt pathway_reaction_id_stmt pathway_reaction_stmt
%type <Stmt> protein_cofactor_stmt protein_domain_stmt protein_step_stmt protein_sequence_stmt
%type <IdentVec> ident_list protein_cofactor_decl_args protein_domain_decl_args gen_expr protein_sequence_decl_args
%type <Stmt> using_stmt

%type <Block> experiment_block experiment_stmts
%type <Stmt> experiment_stmt property_stmt

%right T_ARROW T_BIARROW
%left T_PLUS
%left T_OR
%start program

%%

program        : stmts { ProgramBlock = $1; }
               ;

stmts          : stmt { $$ = new NBlock(); $$->Statements.emplace_back($<Stmt>1); }
               | stmts stmt { $1->Statements.emplace_back($<Stmt>2); }
               ;

stmt           : protein_decl T_SEMIC
               | pathway_decl T_SEMIC
               | organism_decl T_SEMIC
               | experiment_decl T_SEMIC
               | using_stmt T_SEMIC
               ;

block          : T_LBRACE stmts T_RBRACE { $$ = $2; }
               | T_LBRACE T_RBRACE { $$ = new NBlock(); }
               ;

organism_decl  : T_ORGANISM ident T_STRING_LITERAL { $$ = new NOrganismDeclaration(*$2, *$3); delete $2; delete $3; }
               ;

experiment_decl : T_EXPERIMENT ident T_STRING_LITERAL { $$ = new NExperimentDeclaration(*$2, *$3); delete $2; delete $3; }
                | T_EXPERIMENT ident experiment_block  { $$ = new NExperimentDeclaration(*$2, $3); delete $2; }
                ;

protein_decl   : T_PROTEIN ident T_LPAREN protein_decl_args T_RPAREN 
                 { $$ = new NProteinDeclaration(*$2, *$4); delete $2; delete $4; }
               | T_PROTEIN ident T_LPAREN protein_decl_args T_RPAREN protein_block
                 { $$ = new NProteinDeclaration(*$2, *$4, $6); delete $2; delete $4; }
               ;

protein_decl_args : gen_reaction_decl_args { $$ = $1; }
                  ;

protein_block  : T_LBRACE protein_stmts T_RBRACE { $$ = $2; }
               | T_LBRACE T_RBRACE { $$ = new NBlock(); }
               ;

protein_stmts  : protein_stmt { $$ = new NBlock(); $$->Statements.emplace_back($<Stmt>1); }
               | protein_stmts protein_stmt { $1->Statements.emplace_back($<Stmt>2); }
               ;

protein_stmt   : protein_cofactor_stmt T_SEMIC { $$ = $1; }
               | protein_domain_stmt T_SEMIC { $$ = $1; }
               | protein_step_stmt T_SEMIC { $$ = $1; }
               | protein_sequence_stmt T_SEMIC { $$ = $1; }
               ;

protein_cofactor_stmt : T_COFACTOR ident T_LPAREN protein_cofactor_decl_args T_RPAREN { $$ = new NProteinCofactorStatement(*$2, *$4); delete $2; delete $4; }
                      ;

protein_domain_stmt : T_DOMAIN ident T_LPAREN protein_domain_decl_args T_RPAREN { $$ = new NProteinDomainStatement(*$2, *$4); delete $2; delete $4; }
                    ;

protein_step_stmt : T_STEP ident T_LPAREN protein_step_decl_args T_RPAREN { $$ = new NProteinStepStatement(*$2, *$4); delete $2; delete $4; }
                  ;

protein_sequence_stmt : T_SEQUENCE ident T_LPAREN protein_sequence_decl_args T_RPAREN { $$ = new NProteinSequenceStatement(*$2, *$4); delete $2; delete $4; }
                      ;

protein_cofactor_decl_args : ident_list { $$ = $1; }
                           ;

protein_domain_decl_args : ident_list { $$ = $1; }
                         ;

protein_step_decl_args : gen_reaction_decl_args { $$ = $1; }
                       ;

protein_sequence_decl_args : ident_list { $$ = $1; }
                           ;

pathway_decl   : T_PATHWAY ident T_EQUAL pathway_decl_args { $$ = new NPathwayDeclaration(*$2, $4); delete $2; }
               | T_PATHWAY ident pathway_block { $$ = new NPathwayDeclaration(*$2, $3); delete $2; }
               ;

pathway_decl_args : /* blank */ { $<Ident>$ = new NIdentifier(); }
                  | pathway_expr 
                  ;

pathway_expr   : ident { $<Ident>$ = new NIdentifier(*$1); delete $1; }
               | pathway_expr T_PLUS pathway_expr  { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
               | pathway_expr T_OR pathway_expr    { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
               | pathway_expr T_ARROW pathway_expr { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
               ;

pathway_block  : T_LBRACE pathway_stmts T_RBRACE { $$ = $2; }
               | T_LBRACE T_RBRACE { $$ = new NBlock(); }
               ;

pathway_stmts  : pathway_stmt { $$ = new NBlock(); $$->Statements.emplace_back($<Stmt>1); }
               | pathway_stmts pathway_stmt { $1->Statements.emplace_back($<Stmt>2); }
               ;

pathway_stmt   : pathway_description_stmt T_SEMIC { $$ = $1; }
               | pathway_reaction_id_stmt T_SEMIC { $$ = $1; }
               | pathway_reaction_stmt T_SEMIC { $$ = $1; }
               ;

pathway_description_stmt : T_DESCRIPTION T_STRING_LITERAL { $$ = new NDescriptionStatement(*$2); delete $2; }
                         ;

pathway_reaction_id_stmt : T_REACTION_ID ident { $$ = new NPathwayReactionIDStatement(*$2); delete $2; }
                         ;

pathway_reaction_stmt    : T_REACTION pathway_decl_args { $$ = new NPathwayReactionStatement(*$2); delete $2; }
                         ;

experiment_block : T_LBRACE experiment_stmts T_RBRACE { $$ = $2; }
                 | T_LBRACE T_RBRACE { $$ = new NBlock(); }
                 ;

experiment_stmts : experiment_stmt { $$ = new NBlock(); $$->Statements.emplace_back($1); }
                 | experiment_stmts experiment_stmt { $1->Statements.emplace_back($2); }
                 ;

experiment_stmt : T_DESCRIPTION T_STRING_LITERAL T_SEMIC { $$ = new NDescriptionStatement(*$2); delete $2; }
                | property_stmt T_SEMIC { $$ = $1; }
                ;

property_stmt  : T_PROPERTY ident T_NUMBER { $$ = new NPropertyStatement($2->Name, *$3); delete $2; delete $3; }
               | T_PROPERTY ident T_STRING_LITERAL { $$ = new NPropertyStatement($2->Name, *$3); delete $2; delete $3; }
               | T_PROPERTY ident ident { $$ = new NPropertyStatement($2->Name, $3->Name); delete $2, delete $3; }
               ;

using_stmt     : T_USING T_MODULE ident { $$ = new NUsingStatement($2, *$3); delete $3; }
               ;

ident_list     : /* blank */ { $$ = new IdentifierList(); }
               | ident { $$ = new IdentifierList(); $$->emplace_back($<Ident>1); }
               | ident_list T_COMMA ident { $1->emplace_back($<Ident>3); }
               ;

gen_reaction_decl_args : /* blank */ { $$ = new NReaction(); }
                       | gen_expr T_ARROW gen_expr { $$ = new NReaction(*$1, *$3, false); delete $1; delete $3; }
                       | gen_expr T_BIARROW gen_expr { $$ = new NReaction(*$1, *$3, true); delete $1; delete $3; }
                       ;

gen_expr       : gen_ident { $$ = new IdentifierList(); $$->emplace_back($<Ident>1); }
               | gen_expr T_PLUS gen_ident { $1->emplace_back($<Ident>3); }
               ;

gen_ident      : mol_ident { $$ = $1; }
               | ident { $$ = $1; }
               ;

mol_expr       : mol_ident { $$ = new MoleculeList(); $$->emplace_back($<MolIdent>1); }
               | mol_expr T_PLUS mol_ident { $1->emplace_back($<MolIdent>3); }
               ;
	
mol_ident      : T_ATP                  { $$ = new NMoleculeIdentifier(T_ATP, *$1); delete $1; }
               | T_CO2                  { $$ = new NMoleculeIdentifier(T_CO2, *$1); delete $1; }
               | T_COA                  { $$ = new NMoleculeIdentifier(T_COA, *$1); delete $1; }
               | T_OXALOACETATE         { $$ = new NMoleculeIdentifier(T_OXALOACETATE, *$1); delete $1; }
               | T_ACETYL_COA           { $$ = new NMoleculeIdentifier(T_ACETYL_COA, *$1); delete $1; }
               | T_CITRATE              { $$ = new NMoleculeIdentifier(T_CITRATE, *$1); delete $1; }
               | T_ISOCITRATE           { $$ = new NMoleculeIdentifier(T_ISOCITRATE, *$1); delete $1; }
               | T_KETO_GLUTARATE       { $$ = new NMoleculeIdentifier(T_KETO_GLUTARATE, *$1); delete $1; }
               | T_SUCCINATE            { $$ = new NMoleculeIdentifier(T_SUCCINATE, *$1); delete $1; }
               | T_FUMARATE             { $$ = new NMoleculeIdentifier(T_FUMARATE, *$1); delete $1; }
               | T_MALATE               { $$ = new NMoleculeIdentifier(T_MALATE, *$1); delete $1; }
               ;

ident          : T_IDENTIFIER { $$ = new NIdentifier(*$1); delete $1; } 
               ;
                 
%%
