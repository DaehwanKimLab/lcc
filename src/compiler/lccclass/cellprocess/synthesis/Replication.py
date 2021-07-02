# dna replication: DNA Polymerase III holozyme,  DNA Polymerase II (back-up)
'''
https://en.wikipedia.org/wiki/DNA_polymerase_III_holoenzyme
The replisome is composed of the following:

2 DNA Pol III enzymes, each comprising α, ε and θ subunits. (It has been proven that there is a third copy of Pol III at the replisome.[1])
the α subunit (encoded by the dnaE gene) has the polymerase activity.
the ε subunit (dnaQ) has 3'→5' exonuclease activity.
the θ subunit (holE) stimulates the ε subunit's proofreading.
2 β units (dnaN) which act as sliding DNA clamps, they keep the polymerase bound to the DNA.
2 τ units (dnaX) which act to dimerize two of the core enzymes (α, ε, and θ subunits).
1 γ unit (also dnaX) which acts as a clamp loader for the lagging strand Okazaki fragments, helping the two β subunits to form a unit and bind to DNA. The γ unit is made up of 5 γ subunits which include 3 γ subunits, 1 δ subunit (holA), and 1 δ' subunit (holB). The δ is involved in copying of the lagging strand.
Χ (holC) and Ψ (holD) which form a 1:1 complex and bind to γ or τ. X can also mediate the switch from RNA primer to DNA.[2]

1000 nucleotides per second
'''

'''
When Chr bp + 1,
+ 2 PPi
- 1 (use ACGT (frequency table for stoichiometry))

rate: + 1000 bp / 1 second
'''

def SetUpReactions(ProGen):
    Reactions = list()

    Reaction = dict()
    # Reaction No.1

    RXNType = 'Biological Event'
    RXNEquation = '2 dNTP -> 1 ChromosomeSize + 2 PPi'
    RXNRate = '600~1000 events per second'
    RXNTrigger = 'DNA polymerase III, core enzyme >= 4'

    # # The final product would be the following:
    Type = 'Biological Event'

    Stoich_MolIDs = ['dNTP', 'Chromosome1_Rep2', 'PPI[c]']
    Stoich_Coeffs = [-2, 1, 2]

    Rate_Min = 600
    Rate_Max = 1000

    Trigger_MolIDs = ['CPLX0-2361'] # DNA polymerase III, core enzyme
    Trigger_Thresholds = ['4']
    Trigger_Conditions = ['>='] # Greater than or equal to

    # Generate a reaction dictionary with above inputs
    Reaction['Type'] = Type
    Reaction['Stoichiometry'] = [Stoich_MolIDs, Stoich_Coeffs]
    Reaction['Rate'] = [Rate_Min, Rate_Max]
    Reaction['Trigger'] = [Trigger_MolIDs, Trigger_Thresholds, Trigger_Conditions]

    Reaction_SetUp = ProGen.SetUpReaction(Reaction)
    # MolIdxs have been added to Stoichiometry and Trigger variables
    Reactions.append(Reaction_SetUp)

    return Reactions

def Write_CellProcess(Writer, ProGen, ProcessID, Reactions):
    ProGen.GenerateCellProcess(Writer, ProcessID, Reactions)
