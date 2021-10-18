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

	std::vector<std::shared_ptr<NMoleculeIdentifier>> *MolVec;
	std::vector<std::shared_ptr<NIdentifier>> *IdentVec;
	std::string *String;
	int Token;
}

%token <Token> T_PROTEIN T_PROTEIN_COMPLEX T_PATHWAY T_EXPERIMENT T_ORGANISM
%token <Token> T_DESCRIPTION T_REACTION T_REACTION_ID
%token <Token> T_PROPERTY T_USING T_MODULE
%token <Token> T_COFACTOR T_DOMAIN T_STEP T_SEQUENCE

%token <String> T_STRING_LITERAL

%token <String> T_ATP T_CO2 T_OXALOACETATE T_ACETYL_COA T_COA
%token <String> T_CITRATE T_ISOCITRATE T_KETO_GLUTARATE
%token <String> T_SUCCINATE T_FUMARATE T_MALATE
%token <String> T_M_1_AMINO_PROPAN_2_OL
%token <String> T_M_1_AMINO_PROPAN_2_ONE_3_PHOSPHATE
%token <String> T_M_1_CHLORO_24_DINITROBENZENE
%token <String> T_M_1_DEOXYXYLONOJIRIMYCIN
%token <String> T_M_1_ETHYLADENINE
%token <String> T_M_1_KETO_2_METHYLVALERATE
%token <String> T_M_1_PALMITOYLGLYCEROL_3_PHOSPHATE
%token <String> T_M_10_FORMYL_THF
%token <String> T_M_11_DEOXY_CORTISOL
%token <String> T_M_11_DEOXYCORTICOSTERONE
%token <String> T_M_15_DIDEOXY_15_IMINO_D_GALACTITOL
%token <String> T_M_17_DIAMINOHEPTANE
%token <String> T_M_2_3_DIHYDROXYBENZOATE
%token <String> T_M_2_3_DIHYDROXYPHENYL_PROPIONATE
%token <String> T_M_2_5_PHOSPHORIBOSYL_3_DEPHOSPHO_COA
%token <String> T_M_2_5_TRIPHOSPHORIBOSYL_3_DEPHOSPHO_
%token <String> T_M_2_ACETO_2_HYDROXY_BUTYRATE
%token <String> T_M_2_ACETO_LACTATE
%token <String> T_M_2_ALPHA_HYDROXYETHYL_THPP
%token <String> T_M_2_AMINO_3_OXO_4_PHOSPHONOOXYBUTYRATE
%token <String> T_M_2_AMINOACRYLATE
%token <String> T_M_2_AMINOMALONATE_SEMIALDEHYDE
%token <String> T_M_2_C_METHYL_D_ERYTHRITOL_4_PHOSPHATE
%token <String> T_M_2_CARBOXYMUCONATE
%token <String> T_M_2_D_THREO_HYDROXY_3_CARBOXY_ISOCAPROATE
%token <String> T_M_2_DEHYDRO_3_DEOXY_D_GALACTONATE
%token <String> T_M_2_DEHYDRO_3_DEOXY_D_GLUCONATE
%token <String> T_M_2_DEHYDROPANTOATE
%token <String> T_M_2_DEHYDROPANTOYL_LACTONE
%token <String> T_M_2_DEOXY_D_GLUCOSE
%token <String> T_M_2_DEOXY_D_GLUCOSE_6_PHOSPHATE
%token <String> T_M_2_DEOXYGLUCOSE_6_PHOSPHATE
%token <String> T_M_2_DEOXYRIBOSE
%token <String> T_M_2_DH_3_DO_D_ARABINONATE
%token <String> T_M_2_HYDROXY_2_METHYLPROPANENITRILE
%token <String> T_M_2_KETO_3_DEOXY_6_P_GLUCONATE
%token <String> T_M_2_KETO_3_DEOXY_D_GLUCARATE
%token <String> T_M_2_KETO_3_METHYL_VALERATE
%token <String> T_M_2_KETO_ISOVALERATE
%token <String> T_M_2_KETOGLUTARATE
%token <String> T_M_2_MERCAPTOETHANOL
%token <String> T_M_2_METHYLMALEATE
%token <String> T_M_2_O_ALPHA_MANNOSYL_D_GLYCERATE
%token <String> T_M_2_OCTAPRENYL_6_HYDROXYPHENOL
%token <String> T_M_2_OCTAPRENYL_6_METHOXYPHENOL
%token <String> T_M_2_OCTAPRENYLPHENOL
%token <String> T_M_2_OXOBUTANOATE
%token <String> T_M_2_PG
%token <String> T_M_2_PHOSPHO_4_CYTIDINE_5_DIPHOSPHO_2_C_MET
%token <String> T_M_2_PROTOCATECHUOYLPHLOROGLUCINOLCARBOXYLA
%token <String> T_M_2_THIOURIDINE
%token <String> T_M_23_DIMERCAPTOPROPAN_1_OL
%token <String> T_M_23_DIPHOSPHOGLYCERATE
%token <String> T_M_23_PENTANEDIONE
%token <String> T_M_246_TRINITROBENZENE
%token <String> T_M_25_DIDEHYDRO_D_GLUCONATE
%token <String> T_M_2C_METH_D_ERYTHRITOL_CYCLODIPHOSPHATE
%token <String> T_M_2K_4CH3_PENTANOATE
%token <String> T_M_3_4_DIHYDROXYBENZOATE
%token <String> T_M_3_5_ADP
%token <String> T_M_3_AMINO_12_PROPANEDIOL
%token <String> T_M_3_AMINOPROPYL_METHYL_PHOSPHINATE
%token <String> T_M_3_BETA_D_GLUCOSYLGLUCOSE
%token <String> T_M_3_BROMOPYRUVATE
%token <String> T_M_3_CARBOXY_3_HYDROXY_ISOCAPROATE
%token <String> T_M_3_CHLORO_D_ALANINE
%token <String> T_M_3_DEHYDRO_SHIKIMATE
%token <String> T_M_3_DEOXY_D_ARABINO_HEPTULOSONATE_7_P
%token <String> T_M_3_ENOLPYRUVYL_SHIKIMATE_5P
%token <String> T_M_3_ETHYLCATECHOL
%token <String> T_M_3_HYDROXY_PROPIONATE
%token <String> T_M_3_HYDROXYADIPYL_COA
%token <String> T_M_3_HYDROXYBENZOATE
%token <String> T_M_3_HYDROXYPHENYL_PROPIONATE
%token <String> T_M_3_HYDROXYPHENYLACETATE
%token <String> T_M_3_HYDROXYPHENYLACETYL_COA
%token <String> T_M_3_KETO_ADIPYL_COA
%token <String> T_M_3_KETO_L_GULONATE
%token <String> T_M_3_KETOBUTYRATE
%token <String> T_M_3_MERCAPTO_PYRUVATE
%token <String> T_M_3_METHYL_CROTONYL_COA
%token <String> T_M_3_Methyl_Adenines
%token <String> T_M_3_NITROBENZALDEHYDE
%token <String> T_M_3_OCTAPRENYL_4_HYDROXYBENZOATE
%token <String> T_M_3_OXOPALMITOYL_COA
%token <String> T_M_3_P_HYDROXYPYRUVATE
%token <String> T_M_3_P_SERINE
%token <String> T_M_3_PHENYLPROPIONATE
%token <String> T_M_3_SULFINOALANINE
%token <String> T_M_34_DIHYDROXYPHENYLACETALDEHYDE
%token <String> T_M_34_DIHYDROXYPHENYLACETYL_COA
%token <String> T_M_3FE_4S
%token <String> T_M_3OH_4P_OH_ALPHA_KETOBUTYRATE
%token <String> T_M_4_AMINO_4_DEOXY_L_ARABINOSE
%token <String> T_M_4_AMINO_4_DEOXYCHORISMATE
%token <String> T_M_4_AMINO_BUTYRALDEHYDE
%token <String> T_M_4_AMINO_BUTYRATE
%token <String> T_M_4_CARBOXYBENZALDEHYDE
%token <String> T_M_4_CYTIDINE_5_DIPHOSPHO_2_C
%token <String> T_M_4_HYDROXY_2_KETOVALERATE
%token <String> T_M_4_HYDROXY_BENZYL_ALCOHOL
%token <String> T_M_4_HYDROXY_BUTYRATE
%token <String> T_M_4_HYDROXY_L_PROLINE
%token <String> T_M_4_HYDROXYBENZALDEHYDE
%token <String> T_M_4_HYDROXYPHENYLACETATE
%token <String> T_M_4_HYDROXYPHENYLACETYL_COA
%token <String> T_M_4_OXALOMESACONATE
%token <String> T_M_4_P_PANTOTHENATE
%token <String> T_M_4_PHOSPHONOOXY_THREONINE
%token <String> T_M_4_hydroxybenzoate
%token <String> T_M_4FE_4S
%token <String> T_M_5_10_METHENYL_THF
%token <String> T_M_5_AMINO_5_DEOXY_D_GALACTOPYRANOSIDE
%token <String> T_M_5_AMINO_LEVULINATE
%token <String> T_M_5_AMINOPENTANOATE
%token <String> T_M_5_BETA_L_THREO_PENTAPYRANOSYL_4_ULOSE_
%token <String> T_M_5_CHLOROFORMYCIN
%token <String> T_M_5_DEHYDROGLUCONATE
%token <String> T_M_5_ETHYLTHIOADENOSINE
%token <String> T_M_5_FORMYL_THF
%token <String> T_M_5_HYDROXY_CTP
%token <String> T_M_5_HYDROXY_PIPECOLATE
%token <String> T_M_5_IODOPENTAPHOSPHONATE
%token <String> T_M_5_KETO_4_DEOXY_D_GLUCARATE
%token <String> T_M_5_METHYL_THF
%token <String> T_M_5_METHYLAMINOMETHYL_2_SELENOURIDINE
%token <String> T_M_5_METHYLAMINOMETHYL_2_THIOURIDINE
%token <String> T_M_5_METHYLTHIOADENOSINE
%token <String> T_M_5_METHYLTHIOFORMYCIN
%token <String> T_M_5_METHYLTHIOTUBERCIDIN
%token <String> T_M_5_N_PROPYLTHIOADENOSINE
%token <String> T_M_5_OXOPROLINE
%token <String> T_M_5_P_BETA_D_RIBOSYL_AMINE
%token <String> T_M_5_P_NITROPHENYLTHIOADENOSINE
%token <String> T_M_5_P_RIBOSYL_N_FORMYLGLYCINEAMIDE
%token <String> T_M_5_PHOSPHO_RIBOSYL_GLYCINEAMIDE
%token <String> T_M_5_PHOSPHORIBOSYL_5_AMINOIMIDAZOLE
%token <String> T_M_5_PHOSPHORIBOSYL_N_FORMYLGLYCINEAMIDINE
%token <String> T_M_6_AMINOPENICILLANATE
%token <String> T_M_6_DEOXY_D_GLUCOSE
%token <String> T_M_6_PYRUVOYL_5678_TETRAHYDROPTERIN
%token <String> T_M_6_SELENO_OCTANOATE
%token <String> T_M_6_THIO_OCTANOATE
%token <String> T_M_6R_6_FLUORO_EPSP
%token <String> T_M_6S_6_FLUORO_EPSP
%token <String> T_M_7_8_DIHYDROPTEROATE
%token <String> T_M_7_AMINOMETHYL_7_DEAZAGUANINE
%token <String> T_M_7_CYANO_7_DEAZAGUANINE
%token <String> T_M_7_METHYLGUANOSINE_5_PHOSPHATE
%token <String> T_M_8_AMINO_7_OXONONANOATE
%token <String> T_M_8_HYDROXYDEOXYGUANOSINE_5_TRIPHOSPHAT
%token <String> T_M_8_HYDROXYQUINOLINE
%token <String> T_M_8_SELENO_OCTANOATE
%token <String> T_M_8_THIO_OCTANOATE
%token <String> T_M_ACET
%token <String> T_M_ACETALD
%token <String> T_M_ACETAMIDE
%token <String> T_M_ACETOACETYL_COA
%token <String> T_M_ACETOL
%token <String> T_M_ACETYL_COA
%token <String> T_M_ACETYL_D_GLUCOSAMINYLDIPHOSPHO_UNDECAPRE
%token <String> T_M_ACETYL_GLU
%token <String> T_M_ACETYL_P
%token <String> T_M_ACETYLMALTOSE
%token <String> T_M_ACETYLSERINE
%token <String> T_M_ACRYLYL_COA
%token <String> T_M_ADENINE
%token <String> T_M_ADENINE_RING
%token <String> T_M_ADENOSINE
%token <String> T_M_ADENOSINE5TRIPHOSPHO5ADENOSINE
%token <String> T_M_ADENOSINE_DIPHOSPHATE_RIBOSE
%token <String> T_M_ADENOSYL_HOMO_CYS
%token <String> T_M_ADENOSYL_P4
%token <String> T_M_ADENOSYLCOBALAMIN
%token <String> T_M_ADENOSYLCOBALAMIN_5_P
%token <String> T_M_ADENOSYLCOBINAMIDE
%token <String> T_M_ADENOSYLCOBINAMIDE_GDP
%token <String> T_M_ADENOSYLCOBINAMIDE_P
%token <String> T_M_ADENYLOSUCC
%token <String> T_M_ADIPATE
%token <String> T_M_ADP
%token <String> T_M_ADP_D_GLUCOSE
%token <String> T_M_ADP_D_GLYCERO_D_MANNO_HEPTOSE
%token <String> T_M_ADP_GROUP
%token <String> T_M_ADP_L_GLYCERO_D_MANNO_HEPTOSE
%token <String> T_M_ADP_MANNOSE
%token <String> T_M_AG_
%token <String> T_M_AGMATHINE
%token <String> T_M_AICAR
%token <String> T_M_ALA_GLY
%token <String> T_M_ALL_TRANS_HEPTAPRENYL_DIPHOSPHATE
%token <String> T_M_ALL_TRANS_HEXAPRENYL_DIPHOSPHATE
%token <String> T_M_ALL_TRANS_PENTAPRENYL_DIPHOSPHATE
%token <String> T_M_ALLANTOATE
%token <String> T_M_ALLOLACTOSE
%token <String> T_M_ALPHA_D_GALACTOSE
%token <String> T_M_ALPHA_D_MANNOSYL_3_PHOSPHOGLYCERATE
%token <String> T_M_ALPHA_GLC_6_P
%token <String> T_M_ALPHA_GLUCOSE
%token <String> T_M_ALPHA_GLUCOSE_16_BISPHOSPHATE
%token <String> T_M_ALPHA_METHYLMETHIONINE
%token <String> T_M_ALPHA_NAPHTHYL_ACETATE
%token <String> T_M_ALPHA_RIBAZOLE
%token <String> T_M_ALPHA_RIBAZOLE_5_P
%token <String> T_M_ALTROSE
%token <String> T_M_AMINO_ACETONE
%token <String> T_M_AMINO_GROUP
%token <String> T_M_AMINO_HYDROXYMETHYL_METHYL_PYR_P
%token <String> T_M_AMINO_HYDROXYMETHYL_METHYLPYRIMIDINE_PP
%token <String> T_M_AMINO_OH_HYDROXYMETHYL_DIHYDROPTERIDINE
%token <String> T_M_AMINO_OXOBUT
%token <String> T_M_AMINO_RIBOSYLAMINO_1H_3H_PYR_DIONE
%token <String> T_M_AMINOMALONATE
%token <String> T_M_AMMONIA
%token <String> T_M_AMMONIUM
%token <String> T_M_AMP
%token <String> T_M_AMP_GROUP
%token <String> T_M_AMP_LYSINE
%token <String> T_M_ANTHRANILATE
%token <String> T_M_ANTIMONITE
%token <String> T_M_APS
%token <String> T_M_ARABINOSE
%token <String> T_M_ARABINOSE_5P
%token <String> T_M_ARG
%token <String> T_M_ARSENATE
%token <String> T_M_AS_5
%token <String> T_M_ASCORBATE
%token <String> T_M_ASN
%token <String> T_M_ATABRINE
%token <String> T_M_ATP
%token <String> T_M_Alpha_lactose
%token <String> T_M_B_ALANINE
%token <String> T_M_BA_2
%token <String> T_M_BACIMETHRIN
%token <String> T_M_BENZALDEHYDE
%token <String> T_M_BENZENE
%token <String> T_M_BENZOATE
%token <String> T_M_BENZOYLCOA
%token <String> T_M_BENZYL_ALCOHOL
%token <String> T_M_BETA_D_FRUCTOSE
%token <String> T_M_BETA_D_XYLOSE
%token <String> T_M_BETAINE
%token <String> T_M_BETAINE_ALDEHYDE_HYDRATE
%token <String> T_M_BETAINE_ALDEHYDE
%token <String> T_M_BIO_5_AMP
%token <String> T_M_BIOTIN
%token <String> T_M_BIS_P_NITROPHENOLPHOSPHATE
%token <String> T_M_BISGLYCEROPHOSPHOGLYCEROL
%token <String> T_M_BISOHMYR_GLC
%token <String> T_M_BISOHMYR_GLUCOSAMINYL_1P
%token <String> T_M_BORATE
%token <String> T_M_BOROHYDRIDE
%token <String> T_M_BR_
%token <String> T_M_BROMOACETATE
%token <String> T_M_BUTANAL
%token <String> T_M_BUTANEDIOL
%token <String> T_M_BUTANOL
%token <String> T_M_BUTHIONINE_SULFOXIMINE
%token <String> T_M_BUTYL_HYDROPEROXIDE
%token <String> T_M_BUTYLAMINE
%token <String> T_M_BUTYRIC_ACID
%token <String> T_M_BUTYRYL_COA
%token <String> T_M_C_DI_GMP
%token <String> T_M_C1
%token <String> T_M_C3
%token <String> T_M_C4
%token <String> T_M_C5
%token <String> T_M_C55_PP_GLCNAC_MANNACA
%token <String> T_M_C55_PP_GLCNAC_MANNACA_FUC4NAC
%token <String> T_M_C6
%token <String> T_M_CA_2
%token <String> T_M_CADAVERINE
%token <String> T_M_CAFFEOYLQUINATE
%token <String> T_M_CAMP
%token <String> T_M_CANAVANINE
%token <String> T_M_CAPSAICIN
%token <String> T_M_CARBAMATE
%token <String> T_M_CARBAMOYL_P
%token <String> T_M_CARBAMYUL_L_ASPARTATE
%token <String> T_M_CARBON_DIOXIDE
%token <String> T_M_CARBON_MONOXIDE
%token <String> T_M_CARBOXYETHYL_3_5_CYCLOHEXADIENE_1_2_DIOL
%token <String> T_M_CARBOXYL_GROUP
%token <String> T_M_CARBOXYMETHOXYLAMINE
%token <String> T_M_CARBOXYPHENYLAMINO_DEOXYRIBULOSE_P
%token <String> T_M_CARNITINE
%token <String> T_M_CARNOSINE
%token <String> T_M_CATECHOL
%token <String> T_M_CD_2
%token <String> T_M_CDP
%token <String> T_M_CDP_D_GLUCOSE
%token <String> T_M_CDP_GROUP
%token <String> T_M_CELLOBIOSE
%token <String> T_M_CEPHALOSPORIN_C
%token <String> T_M_CGMP
%token <String> T_M_CH33ADO
%token <String> T_M_CH4
%token <String> T_M_CHLORALAN_CPD
%token <String> T_M_CHLORATE
%token <String> T_M_CHLOROACETALDEHYDE
%token <String> T_M_CHOLANATE2
%token <String> T_M_CHOLATE
%token <String> T_M_CHOLINE
%token <String> T_M_CHORISMATE
%token <String> T_M_CIS_ACONITATE
%token <String> T_M_CIT
%token <String> T_M_CL_
%token <String> T_M_CMP
%token <String> T_M_CMP_GROUP
%token <String> T_M_CMP_KDO
%token <String> T_M_CO_2
%token <String> T_M_CO_A
%token <String> T_M_COA_GROUP
%token <String> T_M_COB_I_ALAMIN
%token <String> T_M_COBINAMIDE
%token <String> T_M_COFORMYCIN
%token <String> T_M_COPROPORPHYRINOGEN_III
%token <String> T_M_CR_3
%token <String> T_M_CR_6
%token <String> T_M_CREATINE_P
%token <String> T_M_CROTONO_BETAINE
%token <String> T_M_CROTONOBETAINYL_COA
%token <String> T_M_CROTONYL_COA
%token <String> T_M_CTP
%token <String> T_M_CU_
%token <String> T_M_CU_2
%token <String> T_M_CYS
%token <String> T_M_CYS_GLY
%token <String> T_M_CYSTINE
%token <String> T_M_CYTIDINE
%token <String> T_M_CYTOSINE
%token <String> T_M_D_4_HYDROXY_2_KETO_GLUTARATE
%token <String> T_M_D_6_P_GLUCONO_DELTA_LACTONE
%token <String> T_M_D_ALA_D_ALA
%token <String> T_M_D_ALANINE
%token <String> T_M_D_ALLOSE_6_PHOSPHATE
%token <String> T_M_D_ALLULOSE_6_PHOSPHATE
%token <String> T_M_D_ALTRONATE
%token <String> T_M_D_BETA_D_HEPTOSE_1_P
%token <String> T_M_D_BETA_D_HEPTOSE_17_DIPHOSPHATE
%token <String> T_M_D_CARNITINE
%token <String> T_M_D_CYSTEINE
%token <String> T_M_D_ERYTHRO_IMIDAZOLE_GLYCEROL_P
%token <String> T_M_D_GALACTARATE
%token <String> T_M_D_GALACTONATE
%token <String> T_M_D_GALACTONO_1_4_LACTONE
%token <String> T_M_D_GALACTOSYLAMINE
%token <String> T_M_D_GLT
%token <String> T_M_D_GLUCARATE
%token <String> T_M_D_LACTATE
%token <String> T_M_D_MANNONATE
%token <String> T_M_D_METHYL_MALONYL_COA
%token <String> T_M_D_MYO_INOSITOL_1_MONOPHOSPHATE
%token <String> T_M_D_RIBULOSE
%token <String> T_M_D_RIBULOSE_1_P
%token <String> T_M_D_SEDOHEPTULOSE_1_7_P2
%token <String> T_M_D_SEDOHEPTULOSE_7_P
%token <String> T_M_D_SERINE
%token <String> T_M_D_SORBITOL_6_P
%token <String> T_M_D_TAGATURONATE
%token <String> T_M_D_TARTRATE
%token <String> T_M_D_THREONINE
%token <String> T_M_D_TRYPTOPHAN
%token <String> T_M_D_TYROSINE
%token <String> T_M_D_XYLONATE
%token <String> T_M_D_XYLULOSE
%token <String> T_M_D_aspartate_peptidoglycan
%token <String> T_M_DADP
%token <String> T_M_DAMP
%token <String> T_M_DATP
%token <String> T_M_DCDP
%token <String> T_M_DCMP
%token <String> T_M_DCPIP
%token <String> T_M_DCTP
%token <String> T_M_DE_O_GLUCONATE
%token <String> T_M_DE_O_K_GLUCONATE
%token <String> T_M_DEAMIDO_NAD
%token <String> T_M_DECOYININE
%token <String> T_M_DEHYDRO_3_DEOXY_L_RHAMNONATE
%token <String> T_M_DEHYDRO_DEOXY_GALACTONATE_PHOSPHATE
%token <String> T_M_DEHYDROQUINATE
%token <String> T_M_DELTA1_PIPERIDEINE_2_6_DICARBOXYLATE
%token <String> T_M_DELTA3_ISOPENTENYL_PP
%token <String> T_M_DEMETHYLMENAQUINONE
%token <String> T_M_DEOXY_D_RIBOSE_1_PHOSPHATE
%token <String> T_M_DEOXYADENOSINE
%token <String> T_M_DEOXYCHOLATE
%token <String> T_M_DEOXYCYTIDINE
%token <String> T_M_DEOXYGUANOSINE
%token <String> T_M_DEOXYINOSINE
%token <String> T_M_DEOXYURIDINE
%token <String> T_M_DEOXYXYLULOSE_5P
%token <String> T_M_DEPHOSPHO_COA
%token <String> T_M_DETHIOBIOTIN
%token <String> T_M_DGDP
%token <String> T_M_DGMP
%token <String> T_M_DGTP
%token <String> T_M_DI_H_OROTATE
%token <String> T_M_DI_H_URACIL
%token <String> T_M_DIACETYL
%token <String> T_M_DIACETYLCHITOBIOSE_6_PHOSPHATE
%token <String> T_M_DIAMINO_OH_PHOSPHORIBOSYLAMINO_PYR
%token <String> T_M_DIAMINONONANOATE
%token <String> T_M_DICOUMAROL
%token <String> T_M_DIETHYLDITHIOCARBAMATE
%token <String> T_M_DIETHYLPYROCARBONATE
%token <String> T_M_DIGLUCO_DOCOSANOATE
%token <String> T_M_DIGLUCOACETYL_DOCOSANOATE
%token <String> T_M_DIGLUCODIACETYL_DOCOSANOATE
%token <String> T_M_DIHYDRO_DIOH_BENZOATE
%token <String> T_M_DIHYDRO_NEO_PTERIN
%token <String> T_M_DIHYDRO_THYMINE
%token <String> T_M_DIHYDROCOUMARIN
%token <String> T_M_DIHYDROFOLATE
%token <String> T_M_DIHYDROLIPOAMIDE
%token <String> T_M_DIHYDROMONAPTERIN_TRIPHOSPHATE
%token <String> T_M_DIHYDRONEOPTERIN_P
%token <String> T_M_DIHYDRONEOPTERIN_P3
%token <String> T_M_DIHYDROPTERIN_CH2OH_PP
%token <String> T_M_DIHYDROSIROHYDROCHLORIN
%token <String> T_M_DIHYDROXY_ACETONE_PHOSPHATE
%token <String> T_M_DIHYDROXY_BUTANONE_P
%token <String> T_M_DIHYDROXYACETONE
%token <String> T_M_DIHYDROXYNAPHTHOATE
%token <String> T_M_DIHYDROXYPENTANEDIONE
%token <String> T_M_DIISOPROPYL_FLUOROPHOSPHATE
%token <String> T_M_DIISOPROPYL_PHOSPHATE
%token <String> T_M_DIMETHYL_D_RIBITYL_LUMAZINE
%token <String> T_M_DIMETHYL_GLYCINE
%token <String> T_M_DIMETHYLARSINATE
%token <String> T_M_DIMETHYLBENZIMIDAZOLE
%token <String> T_M_DIMETHYLSUBERIMIDATE
%token <String> T_M_DIMETHYLSULFONIOACETATE
%token <String> T_M_DIMP
%token <String> T_M_DIPYRROMETHANE
%token <String> T_M_DITHIO_NITROBENZOATE
%token <String> T_M_DITHIOTHREITOL
%token <String> T_M_DITP
%token <String> T_M_DL_5_FLUOROTRYPTOPHAN
%token <String> T_M_DMSO
%token <String> T_M_DODECANOATE
%token <String> T_M_DOLICHOL_GROUP
%token <String> T_M_DOPAMINE
%token <String> T_M_DPG
%token <String> T_M_DTDP_D_GLUCOSE
%token <String> T_M_DTDP_DEOH_DEOXY_GLUCOSE
%token <String> T_M_DTDP_DEOH_DEOXY_MANNOSE
%token <String> T_M_DTDP_GROUP
%token <String> T_M_DTDP_RHAMNOSE
%token <String> T_M_DUDP
%token <String> T_M_DUMP
%token <String> T_M_DUTP
%token <String> T_M_ECTOINE
%token <String> T_M_EDTA
%token <String> T_M_EGTA
%token <String> T_M_ENOL_OXALOACETATE
%token <String> T_M_ENOL_PHENYLPYRUVATE
%token <String> T_M_ENTEROBACTIN
%token <String> T_M_ERYTHRITOL
%token <String> T_M_ERYTHRONATE_4P
%token <String> T_M_ERYTHROSE
%token <String> T_M_ERYTHROSE_4P
%token <String> T_M_ETHANAMINE
%token <String> T_M_ETHANOL_AMINE
%token <String> T_M_ETHIONINE
%token <String> T_M_ETHYL_2_METHYLACETOACETATE
%token <String> T_M_ETHYL_2R_METHYL_3S_HYDROXYBUTANOATE
%token <String> T_M_ETHYL_ACETOACETATE
%token <String> T_M_ETHYLACETATE
%token <String> T_M_ETHYLACETIMIDATE
%token <String> T_M_ETHYLENE_CMPD
%token <String> T_M_ETOH
%token <String> T_M_F_
%token <String> T_M_FAD
%token <String> T_M_FAD_STEM_GROUP
%token <String> T_M_FADH2
%token <String> T_M_FARNESYL_PP
%token <String> T_M_FARNESYLFARNESYLGERANIOL
%token <String> T_M_FE_2
%token <String> T_M_FE_3
%token <String> T_M_FERRIC_CITRATE_COMPLEX
%token <String> T_M_FERRIC_ENTEROBACTIN_COMPLEX
%token <String> T_M_FMN
%token <String> T_M_FMNH2
%token <String> T_M_FMPP
%token <String> T_M_FORMALDEHYDE
%token <String> T_M_FORMAMIDE
%token <String> T_M_FORMATE
%token <String> T_M_FORMYCIN_A
%token <String> T_M_FORMYL_COA
%token <String> T_M_FRU1P
%token <String> T_M_FRUCTOSE_16_DIPHOSPHATE
%token <String> T_M_FRUCTOSE_6P
%token <String> T_M_FRUCTOSEGLYCINE
%token <String> T_M_FRUCTOSELYSINE
%token <String> T_M_FRUCTOSELYSINE_6_PHOSPHATE
%token <String> T_M_FUCULOSE_1P
%token <String> T_M_FUM
%token <String> T_M_G3P
%token <String> T_M_GALACTITOL
%token <String> T_M_GALACTITOL_1_PHOSPHATE
%token <String> T_M_GALACTOSE
%token <String> T_M_GALACTOSE_1P
%token <String> T_M_GAMMA_BUTYROBETAINE
%token <String> T_M_GAMMA_BUTYROBETAINYL_COA
%token <String> T_M_GAMMA_GLUTAMYL_GAMMA_AMINOBUTYRALDEH
%token <String> T_M_GAMMA_GLUTAMYL_PUTRESCINE
%token <String> T_M_GAP
%token <String> T_M_GDP
%token <String> T_M_GDP_4_DEHYDRO_6_DEOXY_D_MANNOSE
%token <String> T_M_GDP_4_DEHYDRO_6_L_DEOXYGALACTOSE
%token <String> T_M_GDP_D_GLUCOSE
%token <String> T_M_GDP_GROUP
%token <String> T_M_GDP_MANNOSE
%token <String> T_M_GDP_TP
%token <String> T_M_GERANIAL
%token <String> T_M_GERANIOL
%token <String> T_M_GERANYL_PP
%token <String> T_M_GERANYLGERANYL_PP
%token <String> T_M_GLC
%token <String> T_M_GLC_1_P
%token <String> T_M_GLC_6_P
%token <String> T_M_GLC_D_LACTONE
%token <String> T_M_GLN
%token <String> T_M_GLOBOMYCIN
%token <String> T_M_GLT
%token <String> T_M_GLUCONATE
%token <String> T_M_GLUCOSAMINATE
%token <String> T_M_GLUCOSAMINE_1P
%token <String> T_M_GLUTACONATE
%token <String> T_M_GLUTAMATE_1_SEMIALDEHYDE
%token <String> T_M_GLUTAMIDE
%token <String> T_M_GLUTARATE
%token <String> T_M_GLUTATHIONE
%token <String> T_M_GLUTATHIONYLSPERMIDINE
%token <String> T_M_GLY
%token <String> T_M_GLYCERALD
%token <String> T_M_GLYCERATE
%token <String> T_M_GLYCEROL
%token <String> T_M_GLYCEROL_3P
%token <String> T_M_GLYCEROPHOSPHOGLYCEROL
%token <String> T_M_GLYCOL
%token <String> T_M_GLYCOLALDEHYDE
%token <String> T_M_GLYCOLLATE
%token <String> T_M_GLYOX
%token <String> T_M_GMP
%token <String> T_M_GMP_LYSINE_PHOSPHORAMIDATE
%token <String> T_M_GTP
%token <String> T_M_GUANINE
%token <String> T_M_GUANOSINE
%token <String> T_M_GUANOSINE_5DP_3DP
%token <String> T_M_GUANOSINE_TETRAPHOSPHATE
%token <String> T_M_H2CO3
%token <String> T_M_HCL
%token <String> T_M_HCN
%token <String> T_M_HCO3
%token <String> T_M_HEPTOSYL_KDO2_LIPID_IVA
%token <String> T_M_HEXANAL
%token <String> T_M_HEXANOATE
%token <String> T_M_HEXANOYL_COA
%token <String> T_M_HG_2
%token <String> T_M_HIS
%token <String> T_M_HISTAMINE
%token <String> T_M_HISTIDINAL
%token <String> T_M_HISTIDINOL
%token <String> T_M_HMP
%token <String> T_M_HOMO_CYS
%token <String> T_M_HOMO_SER
%token <String> T_M_HOMOARGININE
%token <String> T_M_HS
%token <String> T_M_HSCN
%token <String> T_M_HSO3
%token <String> T_M_HYDRAZINE
%token <String> T_M_HYDROGEN_MOLECULE
%token <String> T_M_HYDROGEN_PEROXIDE
%token <String> T_M_HYDROQUINONE
%token <String> T_M_HYDROQUINONE_O_BETA_D_GLUCOPYRANOSIDE
%token <String> T_M_HYDROXY_METHYL_BUTENYL_DIP
%token <String> T_M_HYDROXYL_GROUP
%token <String> T_M_HYDROXYLAMINE
%token <String> T_M_HYDROXYMALONATE
%token <String> T_M_HYDROXYMETHYLBILANE
%token <String> T_M_HYDROXYNAPHTHOQUINONE
%token <String> T_M_HYDROXYPROPANAL
%token <String> T_M_HYDROXYPYRIDINE_N_OXIDE
%token <String> T_M_HYDRPHENYLAC_CPD
%token <String> T_M_HYPOTAURINE
%token <String> T_M_HYPOXANTHINE
%token <String> T_M_IDP
%token <String> T_M_ILE
%token <String> T_M_IMIDAZOLE_ACETOL_P
%token <String> T_M_IMIDAZOLE_LACTATE
%token <String> T_M_IMIDAZOLE_PYRUVATE
%token <String> T_M_IMIDAZOLE_RING
%token <String> T_M_IMINO_GROUP
%token <String> T_M_IMINOASPARTATE
%token <String> T_M_IMP
%token <String> T_M_INDOLE
%token <String> T_M_INDOLE_3_GLYCEROL_P
%token <String> T_M_INDOLE_PYRUVATE
%token <String> T_M_INOSINE
%token <String> T_M_IODOACETAMIDE
%token <String> T_M_IODOACETATE
%token <String> T_M_IPRONIAZID
%token <String> T_M_ISOBUTANOL
%token <String> T_M_ISOCHORISMATE
%token <String> T_M_ISONIAZIDE
%token <String> T_M_ISOVALERATE
%token <String> T_M_ISOVALERYL_COA
%token <String> T_M_ITACONATE
%token <String> T_M_ITP
%token <String> T_M_K_
%token <String> T_M_K_HEXANOYL_COA
%token <String> T_M_KDO
%token <String> T_M_KDO_8P
%token <String> T_M_KDO_LIPID_IVA
%token <String> T_M_KDO2_LAUROYL_LIPID_IVA
%token <String> T_M_KDO2_LIPID_A
%token <String> T_M_KDO2_LIPID_IVA
%token <String> T_M_KDO2_LIPID_IVA_COLD
%token <String> T_M_KDO2_PALMITOLEOYL_LIPID_IVA
%token <String> T_M_KOJIC_ACID
%token <String> T_M_KYNURENATE
%token <String> T_M_L_1_GLYCERO_PHOSPHORYLCHOLINE
%token <String> T_M_L_1_GLYCEROPHOSPHORYLETHANOL_AMINE
%token <String> T_M_L_1_LECITHIN
%token <String> T_M_L_1_LYSOPHOSPHATIDATE
%token <String> T_M_L_2_AMINOHEXANOATE
%token <String> T_M_L_2_AMINOPENTANOIC_ACID
%token <String> T_M_L_ALA_GAMMA_D_GLU_DAP
%token <String> T_M_L_ALLO_THREONINE
%token <String> T_M_L_ALPHA_ALANINE
%token <String> T_M_L_ALPHA_AMINO_EPSILON_KETO_PIMELATE
%token <String> T_M_L_ARA4N_MODIFIED_KDO2_LIPID_A
%token <String> T_M_L_ARABITOL
%token <String> T_M_L_ARGININE_P
%token <String> T_M_L_ARGININO_SUCCINATE
%token <String> T_M_L_ASCORBATE_6_PHOSPHATE
%token <String> T_M_L_ASPARTATE
%token <String> T_M_L_ASPARTATE_SEMIALDEHYDE
%token <String> T_M_L_AZASERINE
%token <String> T_M_L_BETA_ASPARTYL_P
%token <String> T_M_L_CARNITINYL_COA
%token <String> T_M_L_CITRULLINE
%token <String> T_M_L_CYSTATHIONINE
%token <String> T_M_L_CYSTEATE
%token <String> T_M_L_DEHYDRO_ASCORBATE
%token <String> T_M_L_DELTA1_PYRROLINE_5_CARBOXYLATE
%token <String> T_M_L_DI_GMP
%token <String> T_M_L_FUCULOSE
%token <String> T_M_L_GALACTOSE
%token <String> T_M_L_GAMMA_GLUTAMYLCYSTEINE
%token <String> T_M_L_GLUTAMATE_5_P
%token <String> T_M_L_GLUTAMATE_GAMMA_SEMIALDEHYDE
%token <String> T_M_L_GLYCERALDEHYDE
%token <String> T_M_L_GLYCERALDEHYDE_3_PHOSPHATE
%token <String> T_M_L_GULONATE
%token <String> T_M_L_HISTIDINOL_P
%token <String> T_M_L_HOMOCYSTEATE
%token <String> T_M_L_IDONATE
%token <String> T_M_L_LACTATE
%token <String> T_M_L_LYXOSE
%token <String> T_M_L_ORNITHINE
%token <String> T_M_L_PANTOATE
%token <String> T_M_L_PENICILLAMINE
%token <String> T_M_L_PIPECOLATE
%token <String> T_M_L_RHAMNONATE
%token <String> T_M_L_RIBULOSE
%token <String> T_M_L_RIBULOSE_5_P
%token <String> T_M_L_SELENOCYSTEINE
%token <String> T_M_L_THREONINE_O_3_PHOSPHATE
%token <String> T_M_L_XYLULOSE
%token <String> T_M_L_XYLULOSE_5_P
%token <String> T_M_LACTALD
%token <String> T_M_LAUROYLCOA_CPD
%token <String> T_M_LEU
%token <String> T_M_LI_
%token <String> T_M_LINAMARIN
%token <String> T_M_LIPID_IV_A
%token <String> T_M_LIPOAMIDE
%token <String> T_M_LIPOIC_ACID
%token <String> T_M_LIPOL_AMP
%token <String> T_M_LIPOYL_AMP
%token <String> T_M_LL_DIAMINOPIMELATE
%token <String> T_M_LYS
%token <String> T_M_MAL
%token <String> T_M_MALEATE
%token <String> T_M_MALONATE
%token <String> T_M_MALONATE_S_ALD
%token <String> T_M_MALONYL_COA
%token <String> T_M_MALTOHEXAOSE
%token <String> T_M_MALTOPENTAOSE
%token <String> T_M_MALTOTETRAOSE
%token <String> T_M_MALTOTRIOSE
%token <String> T_M_MANNITOL
%token <String> T_M_MANNITOL_1P
%token <String> T_M_MANNOSE_1P
%token <String> T_M_MANNOSE_6P
%token <String> T_M_MELIBIOSE
%token <String> T_M_MELILOTATE
%token <String> T_M_MENADIOL
%token <String> T_M_MESO_DIAMINOPIMELATE
%token <String> T_M_MESO_TARTRATE
%token <String> T_M_MET
%token <String> T_M_METHYL_ACETYLPHOSPHONATE
%token <String> T_M_METHYL_BETA_D_GALACTOSIDE
%token <String> T_M_METHYL_GLYOXAL
%token <String> T_M_METHYL_GROUP
%token <String> T_M_METHYL_MALONYL_COA
%token <String> T_M_METHYLAMINE
%token <String> T_M_METHYLENE_THF
%token <String> T_M_METHYLMETHANETHIOSULFONATE
%token <String> T_M_METHYLNICOTINATE
%token <String> T_M_METOH
%token <String> T_M_MG_2
%token <String> T_M_MI_HEXAKISPHOSPHATE
%token <String> T_M_MI_PENTAKISPHOSPHATE
%token <String> T_M_MN_2
%token <String> T_M_MO_2
%token <String> T_M_MOENOMYCIN
%token <String> T_M_MONOMETHYL_ESTER_OF_TRANS_ACONITATE
%token <String> T_M_MONOTRANS_CIS_DECAPRENYL_GROUP
%token <String> T_M_MYO_INOSITOL
%token <String> T_M_N_23_DIHYDROXYBENZOYL_L_SERINE
%token <String> T_M_N_5_PHOSPHORIBOSYL_ANTHRANILATE
%token <String> T_M_N_ACETYL_D_GLUCOSAMINE
%token <String> T_M_N_ACETYL_D_GLUCOSAMINE_1_P
%token <String> T_M_N_ACETYL_D_MANNOSAMINE
%token <String> T_M_N_ACETYL_D_MANNOSAMINE_6P
%token <String> T_M_N_ACETYL_GLUTAMYL_P
%token <String> T_M_N_ACETYLANTHRANILATE
%token <String> T_M_N_ALPHA_ACETYLORNITHINE
%token <String> T_M_N_BROMOSUCCINIMIDE
%token <String> T_M_N_DELTA_PHOSPHONOACETYL_L_ORNITHINE
%token <String> T_M_N_ETHYLMALEIMIDE
%token <String> T_M_N_FORMYLMETHIONINE
%token <String> T_M_N_METHYLTRYPTOPHAN
%token <String> T_M_N_QUINOLIN_8_YLMETHANESULFONAMI
%token <String> T_M_N_SUCCINYL_2_AMINO_6_KETOPIMELATE
%token <String> T_M_N_SUCCINYLLL_2_6_DIAMINOPIMELATE
%token <String> T_M_N1_ACETYLSPERMINE
%token <String> T_M_N1_METHYLADENINE
%token <String> T_M_N2_SUCCINYLGLUTAMATE
%token <String> T_M_N2_SUCCINYLORNITHINE
%token <String> T_M_N3_METHYLCYTOSINE
%token <String> T_M_N5_METHYLGLUTAMINE
%token <String> T_M_NA_
%token <String> T_M_NACMUR
%token <String> T_M_NAD
%token <String> T_M_NAD_STEM_GROUP
%token <String> T_M_NADH
%token <String> T_M_NADP
%token <String> T_M_NADPH
%token <String> T_M_NAIR
%token <String> T_M_NI_2
%token <String> T_M_NIACINAMIDE
%token <String> T_M_NIACINE
%token <String> T_M_NICOTINAMIDE_NUCLEOTIDE
%token <String> T_M_NICOTINAMIDE_RIBOSE
%token <String> T_M_NICOTINATE_NUCLEOTIDE
%token <String> T_M_NITRATE
%token <String> T_M_NITRIC_OXIDE
%token <String> T_M_NITRITE
%token <String> T_M_NITROUS_OXIDE
%token <String> T_M_NMNH
%token <String> T_M_NOCARDICIN_A
%token <String> T_M_NUKED_PENTOSE_RING
%token <String> T_M_O_PHENANTHROLINE
%token <String> T_M_O_PHOSPHO_L_HOMOSERINE
%token <String> T_M_O_SUCCINYL_L_HOMOSERINE
%token <String> T_M_O_SUCCINYLBENZOATE
%token <String> T_M_OCTANOL
%token <String> T_M_OCTAPRENYL_DIPHOSPHATE
%token <String> T_M_OCTAPRENYL_METHOXY_BENZOQUINONE
%token <String> T_M_OCTAPRENYL_METHYL_METHOXY_BENZQ
%token <String> T_M_OCTAPRENYL_METHYL_OH_METHOXY_BENZQ
%token <String> T_M_OH
%token <String> T_M_OH_HEXANOYL_COA
%token <String> T_M_OH_MYRISTOYL
%token <String> T_M_OH_PYR
%token <String> T_M_OLEATE_CPD
%token <String> T_M_OLEOYL_COA
%token <String> T_M_OROTATE
%token <String> T_M_OROTIDINE_5_PHOSPHATE
%token <String> T_M_OXALACETIC_ACID
%token <String> T_M_OXALATE
%token <String> T_M_OXALO_SUCCINATE
%token <String> T_M_OXALYL_COA
%token <String> T_M_OXAMATE
%token <String> T_M_OXIDIZED_GLUTATHIONE
%token <String> T_M_OXOPENTENOATE
%token <String> T_M_OXYGEN_MOLECULE
%token <String> T_M_P_AMINO_BENZOATE
%token <String> T_M_P_BENZOQUINONE
%token <String> T_M_P_HYDROXY_PHENYLPYRUVATE
%token <String> T_M_P_NITROPHENOL
%token <String> T_M_P_NITROPHENYL_ACETATE
%token <String> T_M_P_RIBOSYL_4_SUCCCARB_AMINOIMIDAZOLE
%token <String> T_M_P3I
%token <String> T_M_P4I
%token <String> T_M_PALMITATE
%token <String> T_M_PALMITYL_COA
%token <String> T_M_PANTETHEINE_P
%token <String> T_M_PANTOTHENATE
%token <String> T_M_PANTOYL_LACTONE
%token <String> T_M_PAPS
%token <String> T_M_PB_2
%token <String> T_M_PENICILLIN_G
%token <String> T_M_PENTACHLOROPHENOL
%token <String> T_M_PENTOSE_RING
%token <String> T_M_PHE
%token <String> T_M_PHENANTHRENE_RING
%token <String> T_M_PHENYL_PYRUVATE
%token <String> T_M_PHENYLACETALDEHYDE
%token <String> T_M_PHENYLACETATE
%token <String> T_M_PHENYLETHYLAMINE
%token <String> T_M_PHENYLGLYOXAL
%token <String> T_M_PHENYLHYDANTOIN
%token <String> T_M_PHENYLHYDRAZINE
%token <String> T_M_PHMB
%token <String> T_M_PHOSPHATE_GROUP
%token <String> T_M_PHOSPHATIDYLETHANOLAMINE_KDO2
%token <String> T_M_PHOSPHO_ENOL_PYRUVATE
%token <String> T_M_PHOSPHONATE
%token <String> T_M_PHOSPHOPANTOTHEINE_GROUP
%token <String> T_M_PHOSPHORIBOSYL_AMP
%token <String> T_M_PHOSPHORIBOSYL_ATP
%token <String> T_M_PHOSPHORIBOSYL_CARBOXY_AMINOIMIDAZOLE
%token <String> T_M_PHOSPHORIBOSYL_FORMAMIDO_CARBOXAMIDE
%token <String> T_M_PHOSPHORIBOSYL_FORMIMINO_AICAR_P
%token <String> T_M_PHOSPHORIBULOSYL_FORMIMINO_AICAR_P
%token <String> T_M_PHOSPHORYL_CHOLINE
%token <String> T_M_PHOSPHORYL_ETHANOLAMINE
%token <String> T_M_PHTHALATE
%token <String> T_M_PI
%token <String> T_M_PICOLINATE
%token <String> T_M_PNP_PROPIONATE
%token <String> T_M_PORPHOBILINOGEN
%token <String> T_M_PORPHYRIN_RING
%token <String> T_M_PPI
%token <String> T_M_PQQ
%token <String> T_M_PRECURSOR_Z
%token <String> T_M_PREPHENATE
%token <String> T_M_PRO
%token <String> T_M_PROPANE_1_2_DIOL
%token <String> T_M_PROPANOL
%token <String> T_M_PROPIONATE
%token <String> T_M_PROPIONYL_COA
%token <String> T_M_PROPIONYL_P
%token <String> T_M_PROTOHEME
%token <String> T_M_PROTON
%token <String> T_M_PROTOPORPHYRINOGEN
%token <String> T_M_PROTOPORPHYRIN_IX
%token <String> T_M_PRPP
%token <String> T_M_PSEUDOURIDINE_5_P
%token <String> T_M_PSICOSE
%token <String> T_M_PSICOSELYSINE
%token <String> T_M_PTERIDINE_RING
%token <String> T_M_PURINE_RING
%token <String> T_M_PUTRESCINE
%token <String> T_M_PYRAZINAMIDE
%token <String> T_M_PYRAZINOIC_ACID
%token <String> T_M_PYRAZOLE
%token <String> T_M_PYRIDINE_RING
%token <String> T_M_PYRIDOXAL
%token <String> T_M_PYRIDOXAL_PHOSPHATE
%token <String> T_M_PYRIDOXAMINE
%token <String> T_M_PYRIDOXAMINE_5P
%token <String> T_M_PYRIDOXINE
%token <String> T_M_PYRIDOXINE_5P
%token <String> T_M_PYRIMIDINE_RING
%token <String> T_M_PYROPHOSPHATE_GROUP
%token <String> T_M_PYRUVATE
%token <String> T_M_QUINATE
%token <String> T_M_QUINOLINATE
%token <String> T_M_R__ALLANTOIN
%token <String> T_M_R_4_PHOSPHOPANTOTHENOYL_L_CYSTEINE
%token <String> T_M_RB_
%token <String> T_M_REDUCED_MENAQUINONE
%token <String> T_M_RHAMNOSE
%token <String> T_M_RHAMNULOSE_1P
%token <String> T_M_RIBITOL
%token <String> T_M_RIBOFLAVIN
%token <String> T_M_RIBOSE_15_BISPHOSPHATE
%token <String> T_M_RIBOSE_1P
%token <String> T_M_RIBOSE_TRIPHOSPHATE
%token <String> T_M_RIBULOSE_5P
%token <String> T_M_S_2_AMINOETHYL_L_CYSTEINE
%token <String> T_M_S_24_DINITROPHENYLGLUTATHIONE
%token <String> T_M_S_3_HYDROXYBUTANOYL_COA
%token <String> T_M_S_ACETYLDIHYDROLIPOAMIDE
%token <String> T_M_S_ACETYLPHOSPHOPANTETHEINE
%token <String> T_M_S_ADENOSYL_4_METHYLTHIO_2_OXOBUTANOATE
%token <String> T_M_S_ADENOSYLMETHIONINAMINE
%token <String> T_M_S_ADENOSYLMETHIONINE
%token <String> T_M_S_ALLANTOIN
%token <String> T_M_S_CARBOXYMETHYL_D_CYSTEINE
%token <String> T_M_S_CITRAMALATE
%token <String> T_M_S_FORMYCINYLHOMOCYSTEINE
%token <String> T_M_S_HYDROXYMETHYLGLUTATHIONE
%token <String> T_M_S_LACTOYL_GLUTATHIONE
%token <String> T_M_S_METHYL_L_CYSTEINE
%token <String> T_M_S_NITROSOGLUTATHIONE
%token <String> T_M_S_TUBERCIDINYLHOMOCYSTEINE
%token <String> T_M_S2O3
%token <String> T_M_S2O4
%token <String> T_M_SARCOSINE
%token <String> T_M_SE_2
%token <String> T_M_SELENALYSINE
%token <String> T_M_SELENATE
%token <String> T_M_SELENITE
%token <String> T_M_SELENOHOMOCYSTEINE
%token <String> T_M_SELENOLIPOATE
%token <String> T_M_SELENOMETHIONINE
%token <String> T_M_SEMICARBAZIDE
%token <String> T_M_SEPO3
%token <String> T_M_SER
%token <String> T_M_SERYL_AMP
%token <String> T_M_SHIKIMATE
%token <String> T_M_SHIKIMATE_5P
%token <String> T_M_SINEFUNGIN
%token <String> T_M_SIROHEME
%token <String> T_M_SIROHYDROCHLORIN
%token <String> T_M_SN_GLYCEROL_1_PHOSPHATE
%token <String> T_M_SO3
%token <String> T_M_SOLANESYL_PYROPHOSPHATE
%token <String> T_M_SORBITOL
%token <String> T_M_SPERMIDINE
%token <String> T_M_SPERMINE
%token <String> T_M_SS_DIMETHYL_BETA_PROPIOTHETIN
%token <String> T_M_STEARIC_ACID
%token <String> T_M_STEAROYL_COA
%token <String> T_M_STERONE_RING
%token <String> T_M_SUC
%token <String> T_M_SUC_COA
%token <String> T_M_SUCC_S_ALD
%token <String> T_M_SUCCINYL_OH_CYCLOHEXADIENE_COOH
%token <String> T_M_SUCROSE
%token <String> T_M_SUCROSE_6P
%token <String> T_M_SULFATE
%token <String> T_M_SULFUR_DIOXIDE
%token <String> T_M_SUPER_OXIDE
%token <String> T_M_T2_C4_DECADIENYL_COA
%token <String> T_M_T2_DECENOYL_COA
%token <String> T_M_TAGATOSE
%token <String> T_M_TAGATOSE_1_6_DIPHOSPHATE
%token <String> T_M_TAGATOSE_6_PHOSPHATE
%token <String> T_M_TARTRATE
%token <String> T_M_TARTRONATE_S_ALD
%token <String> T_M_TAURINE
%token <String> T_M_TDP
%token <String> T_M_TDP_D_FUCOSAMINE
%token <String> T_M_TDP_FUC4NAC
%token <String> T_M_TEREPHTHALATE
%token <String> T_M_TETRADECANOYL_COA
%token <String> T_M_TETRAHYDROTHIOPHENE
%token <String> T_M_TETRAHYDROTHIOPHENE_1_OXIDE
%token <String> T_M_TETRANITROMETHANE
%token <String> T_M_THF
%token <String> T_M_THIABENDAZOLE
%token <String> T_M_THIALYSINE
%token <String> T_M_THIAMINE
%token <String> T_M_THIAMINE_P
%token <String> T_M_THIAMINE_PYROPHOSPHATE
%token <String> T_M_THIOGLYCOLATE
%token <String> T_M_THIOL_GROUP
%token <String> T_M_THR
%token <String> T_M_THREO_DS_ISO_CITRATE
%token <String> T_M_THYMIDINE
%token <String> T_M_THYMINE
%token <String> T_M_THZ
%token <String> T_M_THZ_P
%token <String> T_M_TMP
%token <String> T_M_TRANS_2_HEXENAL
%token <String> T_M_TRANS_CIS_UNDECAPRENYL_GROUP
%token <String> T_M_TREHALOSE
%token <String> T_M_TREHALOSE_6P
%token <String> T_M_TRIMENTHLAMINE_N_O
%token <String> T_M_TRIMETHYLAMINE
%token <String> T_M_TRIMETHYLSULFONIUM
%token <String> T_M_TRIS
%token <String> T_M_TRP
%token <String> T_M_TTP
%token <String> T_M_TUNGSTATE
%token <String> T_M_TYR
%token <String> T_M_TYRAMINE
%token <String> T_M_UBIQUINONE_10
%token <String> T_M_UBIQUINONE_2
%token <String> T_M_UBIQUINONE_6
%token <String> T_M_UBIQUINONE_8
%token <String> T_M_UDP
%token <String> T_M_UDP_4_AMINO_4_DEOXY_L_ARABINOSE
%token <String> T_M_UDP_AA_GLUTAMATE
%token <String> T_M_UDP_AAGM_DIAMINOHEPTANEDIOATE
%token <String> T_M_UDP_ACETYL_CARBOXYVINYL_GLUCOSAMINE
%token <String> T_M_UDP_D_GALACTO_14_FURANOSE
%token <String> T_M_UDP_D_GLUCOSAMINE
%token <String> T_M_UDP_GLUCURONATE
%token <String> T_M_UDP_GROUP
%token <String> T_M_UDP_L_ARA4_FORMYL_N
%token <String> T_M_UDP_MANNAC
%token <String> T_M_UDP_MANNACA
%token <String> T_M_UDP_MURNAC_TETRAPEPTIDE
%token <String> T_M_UDP_N_ACETYL_D_GLUCOSAMINE
%token <String> T_M_UDP_N_ACETYLMURAMATE
%token <String> T_M_UDP_OHMYR_ACETYLGLUCOSAMINE
%token <String> T_M_UDP_OHMYR_GLUCOSAMINE
%token <String> T_M_UMP
%token <String> T_M_UNDECAPRENYL_DIPHOSPHATE
%token <String> T_M_UNDECAPRENYL_P
%token <String> T_M_URACIL
%token <String> T_M_URATE
%token <String> T_M_UREA
%token <String> T_M_URIDINE
%token <String> T_M_UROPORPHYRINOGEN_III
%token <String> T_M_UROPORPHYRIN_III
%token <String> T_M_UTP
%token <String> T_M_V_5
%token <String> T_M_VAL
%token <String> T_M_VANILLATE
%token <String> T_M_X_PHOSPHATE_GROUP
%token <String> T_M_X_PYROPHOSPHATE_GROUP
%token <String> T_M_XANTHINE
%token <String> T_M_XANTHOSINE
%token <String> T_M_XANTHOSINE_5_PHOSPHATE
%token <String> T_M_XDP
%token <String> T_M_XTP
%token <String> T_M_XYLITOL
%token <String> T_M_XYLOSE
%token <String> T_M_XYLULOSE_5_PHOSPHATE
%token <String> T_M_ZN_2
%token <String> T_M_glycogen_monomer
%token <String> T_M_undecaprenyl

%token <String> T_NUMBER T_INTEGER T_DOUBLE T_IDENTIFIER

%token <Token> T_LPAREN T_RPAREN T_LBRACE T_RBRACE T_COMMA T_DOT
%token <Token> T_PLUS T_MINUS T_ARROW T_BIARROW
%token <Token> T_EQUAL T_OR T_SEMIC

%type <Ident> ident gen_ident
%type <Block> program stmts block
%type <Block> pathway_block pathway_stmts
%type <Block> protein_block protein_stmts
%type <Stmt> stmt protein_decl protein_complex_decl pathway_decl
%type <MolVec> mol_expr
%type <MolIdent> mol_ident
%type <Reaction> protein_decl_args
%type <Reaction> gen_reaction_decl_args protein_step_decl_args
%type <ChainReaction> pathway_expr pathway_decl_args
%type <Stmt> organism_decl experiment_decl pathway_stmt protein_stmt
%type <Stmt> pathway_description_stmt pathway_reaction_id_stmt pathway_reaction_stmt
%type <Stmt> protein_cofactor_stmt protein_domain_stmt protein_step_stmt protein_sequence_stmt
%type <IdentVec> ident_list protein_cofactor_decl_args protein_domain_decl_args gen_expr protein_sequence_decl_args
%type <IdentVec> protein_complex_decl_args
%type <Stmt> using_stmt
%type <ChainReaction> chain_reaction_decl_args
%type <ChainReactionExpr> chain_expr

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
               | protein_complex_decl T_SEMIC
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

protein_complex_decl : T_PROTEIN_COMPLEX ident T_EQUAL protein_complex_decl_args { $$ = new NProteinComplexDeclaration(*$2, *$4); delete $2; delete $4; }
                     ;

protein_complex_decl_args : gen_expr { $$ = $1; }
                          ;

pathway_decl   : T_PATHWAY ident T_EQUAL pathway_decl_args { $$ = new NPathwayDeclaration(*$2, $4); delete $2; }
               | T_PATHWAY ident pathway_block { $$ = new NPathwayDeclaration(*$2, $3); delete $2; }
               ;


pathway_decl_args : chain_reaction_decl_args { $$ = $1; }
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

chain_reaction_decl_args : /* blank */ { $$ = new NChainReaction(); }
                         | chain_expr { $$ = new NChainReaction(); $$->Add($1); }
                         | chain_reaction_decl_args T_ARROW chain_expr { $1->Add($3, 1); }
                         | chain_reaction_decl_args T_BIARROW chain_expr { $1->Add($3, 2); }
                         ;

chain_expr : gen_ident { $$ = new NChainReactionExpression(); $$->Add(*$1); delete $1; }
           | chain_expr T_PLUS gen_ident { $1->Add(*$3, 1); delete $3; }
           | chain_expr T_OR gen_ident { $1->Add(*$3, 2); delete $3; }
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
