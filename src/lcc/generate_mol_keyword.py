#!/usr/bin/env python3
import os, sys
from argparse import ArgumentParser, FileType

vimtag = '"<!-- LPP -->'
lextag = '    /* <!-- LPP --> */'
yacctag_token = '/* <!-- LPP-TOKEN --> */'
yacctag_rule = '/* <!-- LPP-RULE --> */'


def generate_keyword(moltsv_fname, output_fname = None):
    if not output_fname or output_fname == "-":
        ofp = sys.stdout
    else:
        ofp = open(output_fname, "w")

    fp = open(moltsv_fname)
    # skip header
    fp.readline()

    mol_names = list()

    for line in fp:
        line = line.strip()
        fields = line.split()

        mol_name = fields[0].strip('"')
        if mol_name.startswith("CPD"):
            continue

        if mol_name.startswith("-"):
            continue

        mol_names.append(mol_name)

    mol_names.sort()
    for mol in mol_names:
        print(mol, file=ofp)


    return


def load_mol_keyword(keyword_fname):

    mol_names = list()

    fp = open(keyword_fname)

    for line in fp:
        line = line.strip()
        mol_names.append(line)

    return mol_names



def generate_vim(keyword_fname, vim_fname, output_fname = None):
    mol_names = load_mol_keyword(keyword_fname)

    if not output_fname or output_fname == "-":
        ofp = sys.stdout
    else:
        ofp = open(output_fname, "w")


    fp = open(vim_fname)

    for line in fp:
        line = line.strip()

        if line.startswith(vimtag):
            for mol in mol_names:
                print("syn match lppMolecule display \"{}\"".format(mol), file=ofp)

        else:
            print(line, file=ofp)


    return


def generate_lex(keyword_fname, lex_fname, output_fname = None):
    mol_names = load_mol_keyword(keyword_fname)

    if not output_fname or output_fname == "-":
        ofp = sys.stdout
    else:
        ofp = open(output_fname, "w")


    fp = open(lex_fname)

    for line in fp:
        line = line.rstrip()

        if line.startswith(lextag):
            for mol in mol_names:
                token_name = "T_M_" + mol.replace('-', '_').replace('+', '_')
                print("\"{}\"              {{ SAVE_TOKEN; return {}; }}".format(mol, token_name), file=ofp)

        else:
            print(line, file=ofp)

    return




def generate_yacc(keyword_fname, yacc_fname, output_fname = None):
    mol_names = load_mol_keyword(keyword_fname)

    if not output_fname or output_fname == "-":
        ofp = sys.stdout
    else:
        ofp = open(output_fname, "w")


    fp = open(yacc_fname)

    for line in fp:
        line = line.rstrip()

        if line.startswith(yacctag_token):
            for mol in mol_names:
                token_name = "T_M_" + mol.replace('-', '_').replace('+', '_')
                print("%token <String> {}".format(token_name), file=ofp)
        elif line.startswith(yacctag_rule):
            for mol in mol_names:
                token_name = "T_M_" + mol.replace('-', '_').replace('+', '_')
                print("               | {}         {{ $$ = new NMoleculeIdentifier({}, *$1); delete $1; }}".format(token_name, token_name), file=ofp)
        else:
            print(line, file=ofp)

    return



if __name__ == "__main__":
    parser = ArgumentParser(description='Generate LCC keyword files')

    parser.add_argument('-o', dest='output_fname', type=str, help='Output filename');
    parser.add_argument('-k', dest='keyword_fname', type=str, help='Keyword list file');

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-m', dest='moltsv_fname', type=str, help='Molecule list tsv file');
    group.add_argument('-v', dest='vim_fname', type=str, help='VIM Syntax template file');
    group.add_argument('-l', dest='lex_fname', type=str, help='Lex template file');
    group.add_argument('-y', dest='yacc_fname', type=str, help='Yacc template file');

    args = parser.parse_args()

    if (args.vim_fname or args.lex_fname or args.yacc_fname) and (not args.keyword_fname):
        parser.print_help()
        exit(1)


    if args.moltsv_fname:
        generate_keyword(args.moltsv_fname, args.output_fname)

    elif args.vim_fname:
        generate_vim(args.keyword_fname, args.vim_fname, args.output_fname)

    elif args.lex_fname:
        generate_lex(args.keyword_fname, args.lex_fname, args.output_fname)

    elif args.yacc_fname:
        generate_yacc(args.keyword_fname, args.yacc_fname, args.output_fname)


