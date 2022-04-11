%{
#include <iostream>
#include <string>
#include <memory>

#include "node.h"
#include "lpp.y.hpp"

#define SAVE_TOKEN  yylval.String = new std::string(yytext, yyleng)
#define TOKEN(t)    (yylval.Token = t)

/* C_COMMENT https://gist.github.com/zneak/2ae80c1053dac67a9e487cb1d1d3c302 */
%}

%option noyywrap
%option yylineno
%option never-interactive
%x C_COMMENT

digit                         [0-9]

%%


[ \t\n]                       { } /* White space */
\"([^\\\"]|\\.)*\"            { yylval.String = new std::string(yytext + 1, yyleng - 2);  // remove first and last \", TODO: remove escape character
                                return T_STRING_LITERAL; } /* String literal */
"//"[^\n]*                    { } /* Single line comment */
"/*"                          { BEGIN(C_COMMENT); }
<C_COMMENT>"*/"               { BEGIN(INITIAL); }
<C_COMMENT>.                  { }
<C_COMMENT>\n                 { }
"for"                         { return TOKEN(T_FOR); }
"while"                       { return TOKEN(T_WHILE); }
"if"                          { return TOKEN(T_IF); }
"else"                        { return TOKEN(T_ELSE); }
"int"                         { SAVE_TOKEN; return T_INT; }
"float"                       { SAVE_TOKEN; return T_FLOAT; }
"array"                       { SAVE_TOKEN; return T_ARRAY; }
"dict"                        { SAVE_TOKEN; return T_DICT; }
"and"                         { return TOKEN(T_AND); }
"or"                          { return TOKEN(T_L_OR); }
"not"                         { return TOKEN(T_NOT); }

"protein"                     { return TOKEN(T_PROTEIN); }
"protein_complex"             { return TOKEN(T_PROTEIN_COMPLEX); }
"pathway"                     { return TOKEN(T_PATHWAY); }
"process"                     { return TOKEN(T_PROCESS); }
"experiment"                  { return TOKEN(T_EXPERIMENT); }
"organism"                    { return TOKEN(T_ORGANISM); }
"ecoli"                       { return TOKEN(T_ECOLI); }
"description"                 { return TOKEN(T_DESCRIPTION); }
"reaction_id"                 { return TOKEN(T_REACTION_ID); }
"reaction"                    { return TOKEN(T_REACTION); }
"property"                    { return TOKEN(T_PROPERTY); }
"using"                       { return TOKEN(T_USING); }
"module"                      { return TOKEN(T_MODULE); }
"cofactor"                    { return TOKEN(T_COFACTOR); }
"domain"                      { return TOKEN(T_DOMAIN); }
"PDB"                         { return TOKEN(T_PDB); }
"step"                        { return TOKEN(T_STEP); }
"sequence"                    { return TOKEN(T_SEQUENCE); } 
"threshold"                   { return TOKEN(T_THRESHOLD); }
"replication_origin"          { return TOKEN(T_REPLICATION_ORIGIN); }
"replication_terminus"        { return TOKEN(T_REPLICATION_TERMINUS); }
"ribosome_binding_site"       { return TOKEN(T_RIBOSOME_BINDING_SITE); }
"translation_terminator"      { return TOKEN(T_TRANSLATION_TERMINATOR); }
"polymerase"                  { return TOKEN(T_POLYMERASE); }
"ribosome"                    { return TOKEN(T_POLYMERASE); }
"initiation"                  { return TOKEN(T_INITIATION); }
"elongation"                  { return TOKEN(T_ELONGATION); }
"termination"                 { return TOKEN(T_TERMINATION); }
"container"                   { return TOKEN(T_CONTAINER); }
"petridish"                   { return TOKEN(T_PETRIDISH); }
"dna_{n}"                     { SAVE_TOKEN; return T_MOLECULE; }
"dna_{n+1}"                   { SAVE_TOKEN; return T_MOLECULE; }
"dnt"|"DNT"|"dNTP"            { SAVE_TOKEN; return T_MOLECULE; }
"rna_{n}"                     { SAVE_TOKEN; return T_MOLECULE; }
"rna_{n+1}"                   { SAVE_TOKEN; return T_MOLECULE; }
"nt"|"NT"|"NTP"               { SAVE_TOKEN; return T_MOLECULE; }
"ppi"                         { SAVE_TOKEN; return T_MOLECULE; }
"peptide_{n}"                 { SAVE_TOKEN; return T_MOLECULE; }
"peptide_{n+1}"               { SAVE_TOKEN; return T_MOLECULE; }
"aa"|"AA"                     { SAVE_TOKEN; return T_MOLECULE; }
"H_{2}O"                      { SAVE_TOKEN; return T_MOLECULE; }
"("                           { return TOKEN(T_LPAREN); }
")"                           { return TOKEN(T_RPAREN); }
"+"                           { return TOKEN(T_PLUS); }
"++"                          { return TOKEN(T_INC); }
"-"                           { return TOKEN(T_MINUS); }
"--"                          { return TOKEN(T_DEC); }
"-->"                         { return TOKEN(T_ARROW); }
"--|"                         { return TOKEN(T_INARROW); }
","                           { return TOKEN(T_COMMA); }
"."                           { return TOKEN(T_DOT); }
"<-->"                        { return TOKEN(T_BIARROW); }
"="                           { return TOKEN(T_ASSIGN); }
"|"                           { return TOKEN(T_OR); }
"&"                           { return TOKEN(T_AND); }
";"                           { return TOKEN(T_SEMIC); }
":"                           { return TOKEN(T_COLON); }
"{"                           { return TOKEN(T_LBRACE); }
"}"                           { return TOKEN(T_RBRACE); }
">"                           { return TOKEN(T_GT); }
"<"                           { return TOKEN(T_LT); }
">="                          { return TOKEN(T_GE); }
"<="                          { return TOKEN(T_LE); }
"=="                          { return TOKEN(T_EQ); }
"!="                          { return TOKEN(T_NE); }
"*"                           { return TOKEN(T_STAR); }
"/"                           { return TOKEN(T_DIV); }
"["                           { return TOKEN(T_LBRACKET); }
"]"                           { return TOKEN(T_RBRACKET); }
[a-zA-Z_][a-zA-Z0-9_\-.\{\}\^]*    { SAVE_TOKEN; return T_IDENTIFIER; }
{digit}+                      { SAVE_TOKEN; return T_INTEGER; }
[+-]?{digit}+("."{digit}+)?([Ee][+-]?{digit}+)?    { SAVE_TOKEN; return T_NUMBER; }

(([\(\[\{}]*\w+['\+\-\>']*[\)\]\}]*|[\(\[\{}]+['\+\-\>']+[\)\]\}])+?[\,\-]?\s?)+?	{ SAVE_TOKEN; return T_MOLECULE; } 	/*

	add_brackets_beginning = "[\(\[\{}]"
	add_words = "[A-Za-z\d]"
	add_characters = "[\'\+\-\>']"
	add_brackets_end = "[\)\]\}]"
	add_comma_or_dash = "[\,\-]"
	add_space = "\s"
	pattern_words = add_brackets_beginning+"*"+add_words+"+"+add_characters+"*"+add_brackets_end+"*"
	pattern_characters = add_brackets_beginning+"+"+add_characters+"+"+add_brackets_end+"+"
	compile_pattern = "("+"("+pattern_words+"|"+pattern_characters+")+?"+add_comma_or_dash+"?"+add_space+"?"+")+?"		*/

.                             { printf("Unknown token!\n"); yyterminate(); }

%%
