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
%x C_COMMENT

%%

[ \t\n]                       { } /* White space */
\"([^\\\"]|\\.)*\"            { yylval.String = new std::string(yytext + 1, yyleng - 2);  // remove first and last \", TODO: remove escape character
                                return T_STRING_LITERAL; } /* String literal */
"//".*$                       { } /* Single line comment */
"/*"                          { BEGIN(C_COMMENT); }
<C_COMMENT>"*/"               { BEGIN(INITIAL); }
<C_COMMENT>.                  { }
<C_COMMENT>\n                 { }
"protein"                     { return TOKEN(T_PROTEIN); }
"protein_complex"             { return TOKEN(T_PROTEIN_COMPLEX); }
"pathway"                     { return TOKEN(T_PATHWAY); }
"experiment"                  { return TOKEN(T_EXPERIMENT); }
"organism"                    { return TOKEN(T_ORGANISM); }
"description"                 { return TOKEN(T_DESCRIPTION); }
"reaction_id"                 { return TOKEN(T_REACTION_ID); }
"reaction"                    { return TOKEN(T_REACTION); }
"property"                    { return TOKEN(T_PROPERTY); }
"using"                       { return TOKEN(T_USING); }
"module"                      { return TOKEN(T_MODULE); }
"cofactor"                    { return TOKEN(T_COFACTOR); }
"domain"                      { return TOKEN(T_DOMAIN); }
"step"                        { return TOKEN(T_STEP); }
"sequence"                    { return TOKEN(T_SEQUENCE); } 

(([\(\[\{}]*\w+['\+\-\>']*[\)\]\}]*|[\(\[\{}]+['\+\-\>']+[\)\]\}])+?[\,\-]?\s?)+?	{ return TOKEN(T_MOLECULE); } 	/* 			

	add_brackets_beginning = "[\(\[\{}]" 
	add_words = "\w+"
	add_characters = "[\'\+\-\>']"
	add_brackets_end = "[\)\]\}]"
	add_comma_or_dash = "[\,\-]?"
	add_space = "\s?"
	pattern_words = add_brackets_beginning+"*"+add_words+add_characters+"*"+add_brackets_end+"*"
	pattern_characters = add_brackets_beginning+"+"+add_characters+"+"+add_brackets_end 
	compile_pattern = "("+"("+pattern_words+"|"+pattern_characters+")+?"+add_comma_or_dash+add_space+")+?"		*/ 

[a-zA-Z_][a-zA-Z0-9_\-+]*     { SAVE_TOKEN; return T_IDENTIFIER; }
[0-9]+                        { SAVE_TOKEN; return T_INTEGER; }
("-"|"+")?[0-9]+("."[0-9]+)?  { SAVE_TOKEN; return T_NUMBER; }
"("                           { return TOKEN(T_LPAREN); }
")"                           { return TOKEN(T_RPAREN); }
"+"                           { return TOKEN(T_PLUS); }
"-->"                         { return TOKEN(T_ARROW); }
","                           { return TOKEN(T_COMMA); }
"."                           { return TOKEN(T_DOT); }
"<-->"                        { return TOKEN(T_BIARROW); }
"="                           { return TOKEN(T_EQUAL); }
"|"                           { return TOKEN(T_OR); }
";"                           { return TOKEN(T_SEMIC); }
"{"                           { return TOKEN(T_LBRACE); }
"}"                           { return TOKEN(T_RBRACE); }
.                             { printf("Unknown token!\n"); yyterminate(); }

%%
