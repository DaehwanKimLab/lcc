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
+ 2 Pi
- 1 (use ACGT (frequency table for stoichiometry))

rate: + 1000 bp / 1 second
'''

def Write_Transcription_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Transcription_Init():"):
        Writer.Statement("")

        # Get the master index for the count of DNA Pol III holozymes
        ID_DNAPIIIHolozyme = CompilerData.Comp_Complex.ComplexName2ID['DNA polymerase III, holoenzyme']
        Idx_DNAPIIIHolozyme = CompilerData.Comp_Master.ID2Idx_Master[ID_DNAPIIIHolozyme]
        Writer.Variable_('Cel.Idx_DNAPIIIHolozyme', Idx_DNAPIIIHolozyme)

        # Get the master index for the count of DNA Pol III holozymes


def Write_Transcription_Loop(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Transcription_Init():"):
        Writer.Statement("")


