%{
#include <iostream>
#include <string>

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
                                return (T_STRING_LITERAL); } /* String literal */
"//".*$                       { } /* Single line comment */
"/*"                          { BEGIN(C_COMMENT); }
<C_COMMENT>"*/"               { BEGIN(INITIAL); }
<C_COMMENT>.                  { }
<C_COMMENT>\n                 { }
"protein"                     { return TOKEN(T_PROTEIN); }
"pathway"                     { return TOKEN(T_PATHWAY); }
"experiment"                  { return TOKEN(T_EXPERIMENT); }
"organism"                    { return TOKEN(T_ORGANISM); }
"atp"                         { SAVE_TOKEN; return (T_ATP); }
"co2"                         { SAVE_TOKEN; return (T_CO2); }
"CoA"                         { SAVE_TOKEN; return (T_COA); }
"oxaloacetate"                { SAVE_TOKEN; return (T_OXALOACETATE); }
"acetyl-CoA"                  { SAVE_TOKEN; return (T_ACETYL_COA); }
"citrate"                     { SAVE_TOKEN; return (T_CITRATE); }
"isocitrate"                  { SAVE_TOKEN; return (T_ISOCITRATE); }
"keto-glutarate"              { SAVE_TOKEN; return (T_KETO_GLUTARATE); }
"succinate"                   { SAVE_TOKEN; return (T_SUCCINATE); }
"fumarate"                    { SAVE_TOKEN; return (T_FUMARATE); }
"malate"                      { SAVE_TOKEN; return (T_MALATE); }
[a-zA-Z_][a-zA-Z0-9_-]*       { SAVE_TOKEN; return (T_IDENTIFIER); }
[0-9]+                        { SAVE_TOKEN; return (T_INTEGER); }
"("                           { return TOKEN(T_LPAREN); }
")"                           { return TOKEN(T_RPAREN); }
"+"                           { return TOKEN(T_PLUS); }
"-->"                         { return TOKEN(T_ARROW); }
","                           { return TOKEN(T_COMMA); }
"<-->"                        { return TOKEN(T_BIARROW); }
"="                           { return TOKEN(T_EQUAL); }
"|"                           { return TOKEN(T_OR); }
";"                           { return TOKEN(T_SEMIC); }
"{"                           { return TOKEN(T_LBRACE); }
"}"                           { return TOKEN(T_RBRACE); }
.                             { printf("Unknown token!\n"); yyterminate(); }

%%
