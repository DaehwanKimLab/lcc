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
	NChainReaction *ChainReaction;
	NChainReactionExpression *ChainReactionExpr;
    NSubstrate *Substrate;
    NDeclaraionStatement *DeclStmt;
    NInitializerExpression *InitExpr;

	std::vector<std::shared_ptr<NStatement>> *StmtVec;
	std::vector<std::shared_ptr<NMoleculeIdentifier>> *MolVec;
	std::vector<std::shared_ptr<NIdentifier>> *IdentVec;
	std::vector<std::shared_ptr<NProteinDeclaration>> *ProteinDeclVec;
    std::vector<std::shared_ptr<NPropertyStatement>> *PropertyVec;
    std::vector<std::shared_ptr<NSubstrate>> *SubstrateVec;
    std::vector<std::shared_ptr<NExpression>> *ExprVec;
	std::string *String;
	int Token;
}

%token <Token> T_REACTION T_PROTEIN T_PROTEIN_COMPLEX T_PATHWAY T_EXPERIMENT T_ORGANISM T_PROCESS
%token <Token> T_DESCRIPTION T_REACTION_ID
%token <Token> T_PROPERTY T_USING T_MODULE
%token <Token> T_COFACTOR T_DOMAIN T_STEP T_SEQUENCE T_PDB
%token <Token> T_POLYMERASE T_RIBOSOME
%token <Token> T_REPLICATION_ORIGIN T_REPLICATION_TERMINUS T_RIBOSOME_BINDING_SITE T_TRANSLATION_TERMINATOR
%token <Token> T_INITIATION T_ELONGATION T_TERMINATION
%token <Token> T_CONTAINER T_PETRIDISH
%token <Token> T_FOR T_WHILE T_IF T_ELSE
%token <Token> T_INT T_FLOAT T_ARRAY T_DICT T_AND T_L_OR T_NOT
%token <Token> T_GT T_LT T_GE T_LE T_EQ T_NE
%token <Token> T_STAR T_DIV
%token <Token> T_LBRACKET T_RBRACKET

%token <String> T_STRING_LITERAL

%token <String> T_MOLECULE

%token <String> T_NUMBER T_INTEGER T_DOUBLE T_IDENTIFIER

%token <Token> T_LPAREN T_RPAREN T_LBRACE T_RBRACE T_COMMA T_DOT
%token <Token> T_PLUS T_MINUS T_ARROW T_INARROW T_BIARROW
%token <Token> T_ASSIGN T_OR T_SEMIC T_COLON
%token <Token> T_INC T_DEC

%type <Ident> ident gen_ident
%type <Block> program stmts block
%type <Block> pathway_block pathway_stmts
%type <Block> protein_block protein_stmts
%type <Block> process_block process_stmts
%type <Stmt> reaction_decl
%type <Stmt> protein_decl
%type <Stmt> process_decl
%type <MolVec> mol_expr
%type <MolIdent> mol_ident
%type <Reaction> reaction_decl_args
%type <Reaction> protein_decl_args
%type <Reaction> process_decl_args
%type <Reaction> gen_reaction_decl_args step_decl_args
%type <Reaction> gen_reaction_decl_reaction
%type <ChainReaction> pathway_expr pathway_decl_args
%type <Stmt> pathway_stmt protein_stmt process_stmt
%type <Stmt> pathway_description_stmt pathway_reaction_id_stmt pathway_reaction_stmt
%type <Stmt> protein_cofactor_stmt protein_domain_stmt protein_step_stmt protein_sequence_stmt protein_pdb_stmt protein_reaction_stmt
%type <Stmt> process_step_stmt
%type <IdentVec> ident_list protein_cofactor_decl_args protein_domain_decl_args gen_expr protein_sequence_decl_args
%type <IdentVec> protein_complex_decl_args
%type <ChainReaction> chain_reaction_decl_args
%type <ChainReactionExpr> chain_expr
%type <PropertyVec> gen_property_args gen_reaction_decl_property_args
%type <Stmt> gen_property_arg
%type <SubstrateVec> gen_reaction_expr
%type <Substrate> gen_reaction_substrate

%type <Block> experiment_block experiment_stmts
%type <Stmt> property_stmt experiment_stmt

%type <StmtVec> for_stmt if_stmt while_stmt
%type <StmtVec> stmt reaction_decl_stmt protein_decl_stmt protein_complex_decl pathway_decl organism_decl using_stmt experiment_decl reaction_decls protein_decls ribosome_decl_stmt polymerase_decl_stmt process_decl_stmt
%type <StmtVec> ribosome_decl_args ribosome_args polymerase_decl_args polymerase_args
%type <StmtVec> container_stmt petridish_stmt
%type <StmtVec> p_expr_stmt p_expr_decl_stmt
%type <StmtVec> decl_stmt init_declarator_list
%type <Stmt> ribosome_arg polymerase_arg

%type <Stmt> gen_initiation_stmt gen_elongation_stmt gen_termination_stmt
%type <Reaction> gen_elongation_decl_arg
%type <StmtVec> ribosome_binding_site_stmt translation_terminator_stmt replication_origin_stmt replication_terminus_stmt

%type <Expr> p_expr variable p_range_expr p_const_expr
%type <ExprVec> p_expr_list
%type <DeclStmt> init_declarator
%type <Ident> declarator
%type <InitExpr> initializer initializer_list
%type <Token> declaration_specifiers type_specifier
%type <String> unit

%right T_ARROW T_INARROW T_BIARROW
%left T_PLUS
%left T_OR
%start program

%%

program        : stmts { ProgramBlock = $1; }
               ;

stmts          : stmt { $$ = new NBlock(); $$->AddStatment($<StmtVec>1); }
               | stmts stmt { $1->AddStatment($<StmtVec>2); }
               ;

stmt           : reaction_decl_stmt T_SEMIC
               | protein_decl_stmt T_SEMIC
               | protein_complex_decl T_SEMIC
               | pathway_decl T_SEMIC
               | process_decl_stmt T_SEMIC
               | organism_decl T_SEMIC
               | organism_decl
               | experiment_decl T_SEMIC
               | using_stmt T_SEMIC
			   | ribosome_decl_stmt T_SEMIC
			   | polymerase_decl_stmt T_SEMIC
			   | ribosome_binding_site_stmt T_SEMIC
			   | translation_terminator_stmt T_SEMIC
			   | replication_origin_stmt T_SEMIC
               | replication_terminus_stmt T_SEMIC
               | container_stmt T_SEMIC
               | container_stmt
               | petridish_stmt T_SEMIC
               | petridish_stmt
               | for_stmt T_SEMIC
               | for_stmt
               | while_stmt T_SEMIC
               | while_stmt
               | if_stmt T_SEMIC
               | if_stmt
               | p_expr_stmt T_SEMIC
               | decl_stmt T_SEMIC
               ;

block          : T_LBRACE stmts T_RBRACE { $$ = $2; }
               | T_LBRACE T_RBRACE { $$ = new NBlock(); }
               | /* empty */ { $$ = new NBlock(); }
               ;

for_stmt       : T_FOR T_LPAREN p_expr_decl_stmt T_SEMIC p_expr T_SEMIC p_expr T_RPAREN block { $$ = NNodeUtil::InitStatementList(new NLoopStatement(*$3, $5, $7, $9)); delete $3; }
               ;

if_stmt        : T_IF T_LPAREN p_expr T_RPAREN block { $$ = NNodeUtil::InitStatementList(new NIfStatement($3, $5)); }
               | T_IF T_LPAREN p_expr T_RPAREN block T_ELSE block { $$ = NNodeUtil::InitStatementList(new NIfStatement($3, $5, $7)); }
               ;

while_stmt     : T_WHILE T_LPAREN p_expr T_RPAREN block { $$ = NNodeUtil::InitStatementList(new NLoopStatement($3, $5)); }
               ;

container_stmt : T_CONTAINER ident block { $$ = NNodeUtil::InitStatementList(new NContainerStatement(*$2, $3)); delete $2; }
               ;

petridish_stmt : T_PETRIDISH ident block { $$ = NNodeUtil::InitStatementList(new NPetridishStatement(*$2, $3)); delete $2; }
               ;

p_expr_decl_stmt : p_expr_stmt
                 | decl_stmt
                 ;

p_expr_stmt    : p_expr { $$ = NNodeUtil::InitStatementList(new NExpressionStatement($1)); }
               | p_expr_stmt T_COMMA p_expr { $1->emplace_back(new NExpressionStatement($3)); }
               ;

organism_decl  : T_ORGANISM ident T_STRING_LITERAL { $$ = NNodeUtil::InitStatementList(new NOrganismDeclaration(*$2, *$3)); delete $2; delete $3; }
               | T_ORGANISM ident block { $$ = NNodeUtil::InitStatementList(new NOrganismDeclaration(*$2, $3)); delete $2; }
               ;

experiment_decl : T_EXPERIMENT ident T_STRING_LITERAL  { $$ = NNodeUtil::InitStatementList(new NExperimentDeclaration(*$2, *$3)); delete $2; delete $3; }
                | T_EXPERIMENT ident experiment_block  { $$ = NNodeUtil::InitStatementList(new NExperimentDeclaration(*$2, $3)); delete $2; }
                ;

protein_decl_stmt  : T_PROTEIN protein_decls { $$ = $2; }
                   ;

protein_decls  : protein_decl { $$ = NNodeUtil::InitStatementList($1); }
               | protein_decls T_COMMA protein_decl { $1->emplace_back($3); }
               ;

protein_decl   : ident { $$ = new NProteinDeclaration(*$1); delete $1; }
               | ident T_LPAREN protein_decl_args T_RPAREN { $$ = new NProteinDeclaration(*$1, *$3); delete $1; delete $3; }
               | ident T_LPAREN protein_decl_args T_RPAREN protein_block { $$ = new NProteinDeclaration(*$1, *$3, $5); delete $1; delete $3; }
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
               | protein_pdb_stmt T_SEMIC { $$ = $1; }
               | protein_reaction_stmt T_SEMIC { $$ = $1; }
               ;

protein_cofactor_stmt : T_COFACTOR ident T_LPAREN protein_cofactor_decl_args T_RPAREN { $$ = new NProteinCofactorStatement(*$2, *$4); delete $2; delete $4; }
                      ;

protein_domain_stmt : T_DOMAIN ident T_LPAREN protein_domain_decl_args T_RPAREN { $$ = new NProteinDomainStatement(*$2, *$4); delete $2; delete $4; }
                    ;

protein_step_stmt : T_STEP ident T_LPAREN step_decl_args T_RPAREN { $$ = new NStepStatement(*$2, *$4); delete $2; delete $4; }
                  ;

protein_sequence_stmt : T_SEQUENCE ident T_LPAREN protein_sequence_decl_args T_RPAREN { $$ = new NProteinSequenceStatement(*$2, *$4); delete $2; delete $4; }
                      ;

protein_pdb_stmt : T_PDB T_STRING_LITERAL { $$ = new NPropertyStatement("PDB", new NConstantExpression(*$2)); delete $2; }
                 ;

protein_cofactor_decl_args : ident_list { $$ = $1; }
                           ;

protein_domain_decl_args : ident_list { $$ = $1; }
                         ;

step_decl_args : gen_reaction_decl_args { $$ = $1; }
                       ;

protein_sequence_decl_args : ident_list { $$ = $1; }
                           ;

protein_reaction_stmt : T_REACTION ident T_LPAREN gen_reaction_decl_args T_RPAREN { $4->SetID(*$2); delete $2; $$ = $4; }
                      ;

protein_complex_decl : T_PROTEIN_COMPLEX ident T_ASSIGN protein_complex_decl_args { $$ = NNodeUtil::InitStatementList(new NProteinComplexDeclaration(*$2, *$4)); delete $2; delete $4; }
                     ;

protein_complex_decl_args : gen_expr { $$ = $1; }
                          ;

pathway_decl   : T_PATHWAY ident T_ASSIGN pathway_decl_args { $$ = NNodeUtil::InitStatementList(new NPathwayDeclaration(*$2, $4)); delete $2; }
               | T_PATHWAY ident pathway_block { $$ = NNodeUtil::InitStatementList(new NPathwayDeclaration(*$2, $3)); delete $2; }
               ;


pathway_decl_args : chain_reaction_decl_args { $$ = $1; }
                  ;

pathway_expr   : ident { $<Ident>$ = new NIdentifier(*$1); delete $1; }
               | pathway_expr T_PLUS pathway_expr    { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
               | pathway_expr T_OR pathway_expr      { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
               | pathway_expr T_ARROW pathway_expr   { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
               | pathway_expr T_INARROW pathway_expr { $$ = new NPathwayExpression($1, $3, $2); /* delete $1; delete $3; */}
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

process_decl_stmt : T_PROCESS process_decl { $$ = NNodeUtil::InitStatementList($2); }
                  ;

process_decl      : ident { $$ = new NProcessDeclaration(*$1); delete $1; }
                  | ident T_LPAREN process_decl_args T_RPAREN { $$ = new NProcessDeclaration(*$1, *$3); delete $1; delete $3; }
                  | ident T_LPAREN process_decl_args T_RPAREN process_block { $$ = new NProcessDeclaration(*$1, *$3, $5); delete $1; delete $3; }
                  ;

process_decl_args : /* blank */ { $$ = new NReaction(); } /* gen_reaction_decl_args */
                  ;

process_block     : T_LBRACE process_stmts T_RBRACE { $$ = $2; }
                  | T_LBRACE T_RBRACE { $$ = new NBlock(); }
                  ;

process_stmts     : process_stmt { $$ = new NBlock(); $$->AddStatment($<Stmt>1); }
                  | process_stmts process_stmt { $1->AddStatment($<Stmt>2); }
                  ;

process_stmt      : process_step_stmt T_SEMIC { $$ = $1; }
                  ;

process_step_stmt : T_STEP ident T_LPAREN step_decl_args T_RPAREN { $$ = new NStepStatement(*$2, *$4); delete $2; delete $4; }
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

property_stmt  : T_PROPERTY ident T_NUMBER { $$ = new NPropertyStatement($2->Name, new NConstantExpression(*$3)); delete $2; delete $3; }
               | T_PROPERTY ident T_NUMBER unit { $$ = new NPropertyStatement($2->Name, new NConstantExpression(*$3, *$4)); delete $2; delete $3; delete $4; }
               | T_PROPERTY ident T_STRING_LITERAL { $$ = new NPropertyStatement($2->Name, new NConstantExpression(*$3)); delete $2; delete $3; }
               | T_PROPERTY ident ident { $$ = new NPropertyStatement($2->Name, new NVariableExpression(*$3)); delete $2, delete $3; }
               ;

using_stmt     : T_USING T_MODULE ident { $$ = NNodeUtil::InitStatementList(new NUsingStatement($2, *$3)); delete $3; }
               ;


ribosome_binding_site_stmt : T_RIBOSOME_BINDING_SITE ident { $$ = NNodeUtil::InitStatementList(new NRibosomeBindingSite(*$2)); delete $2; }
                           ;

translation_terminator_stmt : T_TRANSLATION_TERMINATOR ident { $$ = NNodeUtil::InitStatementList(new NTranslationTerminator(*$2)); delete $2; }
						    ;

replication_origin_stmt : T_REPLICATION_ORIGIN ident { $$ = NNodeUtil::InitStatementList(new NReplicationOrigin(*$2)); delete $2; }
					    ;

replication_terminus_stmt : T_REPLICATION_TERMINUS ident { $$ = NNodeUtil::InitStatementList(new NReplicationTerminus(*$2)); delete $2; }
						  ;

ribosome_decl_stmt : T_RIBOSOME ident T_LPAREN ribosome_decl_args T_RPAREN { $$ = NNodeUtil::InitStatementList(new NRibosomeDeclaration(*$2, $4)); delete $2; }
                   ;

ribosome_decl_args : ribosome_args { $$ = $1; }
                   ;

ribosome_args  : /* blank */ { $$ = NNodeUtil::InitStatementList(); }
               | ribosome_arg { $$ = NNodeUtil::InitStatementList(); $$->emplace_back($1); }
               | ribosome_args T_COMMA ribosome_arg { $1->emplace_back($3); }
			   ;

ribosome_arg   : gen_initiation_stmt { $$ = $1; }
               | gen_elongation_stmt { $$ = $1; }
			   | gen_termination_stmt { $$ = $1; }
			   ;

polymerase_decl_stmt : T_POLYMERASE ident T_LPAREN polymerase_decl_args T_RPAREN { $$ = NNodeUtil::InitStatementList(new NPolymeraseDeclaration(*$2, $4)); delete $2; }
                     ;

polymerase_decl_args : polymerase_args { $$ = $1; }
                     ;

polymerase_args : /* blank */ { $$ = NNodeUtil::InitStatementList(); }
                | polymerase_arg { $$ = NNodeUtil::InitStatementList(); $$->emplace_back($1); }
                | polymerase_args T_COMMA polymerase_arg { $1->emplace_back($3); }
                ;

polymerase_arg  : gen_initiation_stmt { $$ = $1; }
                | gen_elongation_stmt { $$ = $1; }
		        | gen_termination_stmt { $$ = $1; }
			    ;

gen_initiation_stmt : T_INITIATION T_COLON ident { $$ = new NInitiationStatement(*$3); delete $3; }
                    ;

gen_elongation_stmt : T_ELONGATION T_COLON gen_elongation_decl_arg { $$ = new NElongationStatement(*$3); delete $3; }
                    ;

gen_termination_stmt : T_TERMINATION T_COLON ident { $$ = new NTerminationStatement(*$3); delete $3; }
                     ;

gen_elongation_decl_arg : gen_reaction_decl_reaction { $$ = $1; }
                        ;

ident_list     : /* blank */ { $$ = new IdentifierList(); }
               | ident { $$ = new IdentifierList(); $$->emplace_back($<Ident>1); }
               | ident_list T_COMMA ident { $1->emplace_back($<Ident>3); }
               ;

chain_reaction_decl_args : /* blank */ { $$ = new NChainReaction(); }
                         | chain_expr { $$ = new NChainReaction(); $$->Add($1); }
                         | chain_reaction_decl_args T_ARROW chain_expr { $1->Add($3, 1); }
                         | chain_reaction_decl_args T_INARROW chain_expr { $1->Add($3, 1); }
                         | chain_reaction_decl_args T_BIARROW chain_expr { $1->Add($3, 2); }
                         ;

chain_expr : gen_ident { $$ = new NChainReactionExpression(); $$->Add(*$1); delete $1; }
           | chain_expr T_PLUS gen_ident { $1->Add(*$3, 1); delete $3; }
           | chain_expr T_OR gen_ident { $1->Add(*$3, 2); delete $3; }
           ;

reaction_decl_stmt  : T_REACTION reaction_decls { $$ = $2; }
                    ;

reaction_decls  : reaction_decl { $$ = NNodeUtil::InitStatementList($1); }
                | reaction_decls T_COMMA reaction_decl { $1->emplace_back($3); }
                ;

reaction_decl   : ident { $$ = new NReactionDeclaration(*$1); delete $1; }
                | ident T_LPAREN reaction_decl_args T_RPAREN { $$ = new NReactionDeclaration(*$1, *$3); delete $1; delete $3; }
                ;

reaction_decl_args : gen_reaction_decl_args { $$ = $1; }
                   ;


gen_reaction_decl_args : /* blank */ { $$ = new NReaction(); }
                       | gen_reaction_decl_reaction { $$ = $1; }
					   | gen_reaction_decl_property_args { $$ = new NReaction(); $$->SetProperty(*$1); delete $1; }
					   | gen_reaction_decl_reaction T_COMMA gen_reaction_decl_property_args { $$ = $1; $$->SetProperty(*$3); delete $3; }
					   ;

gen_reaction_decl_reaction : gen_reaction_expr T_ARROW gen_reaction_expr { $$ = new NReaction(*$1, *$3, true); delete $1; delete $3; }
                           | gen_reaction_expr T_INARROW gen_reaction_expr { $$ = new NReaction(*$1, *$3, false); delete $1; delete $3; }
                           | gen_reaction_expr T_BIARROW gen_reaction_expr { $$ = new NReaction(*$1, *$3, true); delete $1; delete $3; }
                           ;

gen_reaction_decl_property_args : gen_property_args { $$ = $1; }
                           ;

gen_property_args : gen_property_arg { $$ = new PropertyList(); $$->emplace_back(static_cast<NPropertyStatement*>($1)); }
                  | gen_property_args T_COMMA gen_property_arg { $1->emplace_back(static_cast<NPropertyStatement*>($3)); }
                  ;

gen_property_arg  : /* */ { $$ = new NPropertyStatement(); }
                  | ident { $$ = new NPropertyStatement($1->Name); delete $1; }
                  | ident T_ASSIGN p_expr { $$ = new NPropertyStatement($1->Name, $3); delete $1; }
                  ;

gen_reaction_expr : gen_reaction_substrate { $$ = new SubstrateList(); $$->emplace_back($1); }
                  | gen_reaction_expr T_PLUS gen_reaction_substrate { $$->emplace_back($3); }
                  ;

gen_reaction_substrate : gen_ident { $$ = new NSubstrate(*$1); delete $1; }
                       | T_INTEGER gen_ident { $$ = new NSubstrate(*$2, *$1); delete $1; delete $2; }
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

mol_ident      : T_MOLECULE { $$ = new NMoleculeIdentifier(T_MOLECULE, *$1); delete $1; }
               ;

ident          : T_IDENTIFIER { $$ = new NIdentifier(*$1); delete $1; }
               ;

p_const_expr   : T_NUMBER { $$ = new NConstantExpression(*$1); delete $1; }
               | T_NUMBER unit { $$ = new NConstantExpression(*$1, *$2); delete $1; delete $2; }
               | T_INTEGER { $$ = new NConstantExpression(*$1); delete $1; }
               | T_INTEGER unit { $$ = new NConstantExpression(*$1, *$2); delete $1; delete $2; }
               ;

p_expr         : /* */ { $$ = new NExpression(); }
               | p_const_expr { $$ = $1; }
               | variable { $$ = $1; }
               | variable T_ASSIGN p_expr { $$ = new NAExpression(T_ASSIGN, $1, $3); }
               | p_expr T_PLUS p_expr { $$ = new NAExpression(T_PLUS, $1, $3); }
               | p_expr T_MINUS p_expr { $$ = new NAExpression(T_MINUS, $1, $3); }
               | p_expr T_STAR p_expr { $$ = new NAExpression(T_STAR, $1, $3); }
               | p_expr T_DIV p_expr { $$ = new NAExpression(T_DIV, $1, $3); }
               | p_expr T_LT p_expr { $$ = new NAExpression(T_LT, $1, $3); }
               | p_expr T_GT p_expr { $$ = new NAExpression(T_GT, $1, $3); }
               | p_expr T_LE p_expr { $$ = new NAExpression(T_LE, $1, $3); }
               | p_expr T_GE p_expr { $$ = new NAExpression(T_GE, $1, $3); }
               | p_expr T_EQ p_expr { $$ = new NAExpression(T_EQ, $1, $3); }
               | p_expr T_NE p_expr { $$ = new NAExpression(T_NE, $1, $3); }
               | p_expr T_OR p_expr { $$ = new NAExpression(T_OR, $1, $3); }
               | p_expr T_AND p_expr { $$ = new NAExpression(T_AND, $1, $3); }
               | p_expr T_LPAREN T_RPAREN { $$ = new NFunctionCallExpression($1); }
               | p_expr T_LPAREN p_expr_list T_RPAREN { $$ = new NFunctionCallExpression($1, $3); }
               | T_LPAREN p_expr T_RPAREN { $$ = $2; }
               ;

variable       : ident { $$ = new NVariableExpression(*$1); delete $1; }
               | ident T_LBRACKET p_range_expr T_RBRACKET { $$ = new NVariableExpression(*$1, $3); delete $1; }
               ;

p_range_expr   : p_expr { $$ = new NRangeExpression($1, nullptr, nullptr); }
               | p_expr T_COLON p_expr { $$ = new NRangeExpression($1, $3, nullptr); }
               | p_expr T_COLON p_expr T_COLON p_expr { $$ = new NRangeExpression($1, $3, $5); }
               ;

p_expr_list    : p_expr { $$ = new ExpressionList(); $$->emplace_back($1); }
               | p_expr_list T_COMMA p_expr { $1->emplace_back($3); }
               ;

/*
declaration_list : declaration
                 | declaration_list declaration
                 ;
*/

decl_stmt      : declaration_specifiers init_declarator_list { NDeclaraionStatement::UpdateType($2, $1); $$ = $2; }
               ;

declaration_specifiers : type_specifier
                       ;

type_specifier : T_FLOAT
               | T_INT
               | T_ARRAY
               | T_DICT
               ;

init_declarator_list : init_declarator { $$ = NNodeUtil::InitStatementList($1); }
                     | init_declarator_list T_COMMA init_declarator { $1->emplace_back($3); }
                     ;

init_declarator : declarator { $$ = new NDeclaraionStatement(); $$->SetIdentifier(*$1); delete $1;}
                | declarator T_ASSIGN initializer { $$ = new NDeclaraionStatement(); $$->SetIdentifier(*$1); $$->SetInitializer($3); delete $1; }
                ;

declarator     : ident
               ;

initializer    : p_expr { $$ = new NInitializerExpression($1); }
               | T_LBRACKET initializer_list T_RBRACKET { $$ = $2; /* $$ = new NInitializerExpression($2); */ }
               ;
 
initializer_list : initializer { $$ = $1; }
                 | initializer_list T_COMMA initializer { $1->Append($3); }
                 ;

unit           : T_IDENTIFIER
               ;
               
%%
