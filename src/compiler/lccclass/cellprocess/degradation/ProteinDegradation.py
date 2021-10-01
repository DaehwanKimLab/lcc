# protein damage
'''
N-degron and C-degron pathways of protein degradation
Alexander Varshavsky
PNAS January 8, 2019 116 (2) 358-366; first published January 8, 2019; https://doi.org/10.1073/pnas.1816596116
https://www.pnas.org/content/pnas/116/2/358.full.pdf

Bacterial fMet/N-degron pathway:
In 2015, it was found that Nt-fMet residues of Degraded bacterial proteins
can act as bacterial N-degrons, termed fMet/N-degrons (Fig. 1C) (38)

Bacterial Leu/N-Degron Pathway:
This pathway comprises the following components: (i) ClpAP, a proteasome-like, ATP-dependent protease;
(ii) ClpS, the 12-kDa Leu/N-recognin that binds to Nt-Leu, -Phe, -Trp, or -Tyr and delivers bound substrates
to the ClpAP protease; (iii) Aat, an L/F-transferase that employs Leu-tRNA or Phe-tRNA as a cosubstrate
to conjugate largely Leu (and occasionally Phe) to the N-termini of proteins bearing Nt-Lys or Nt-Arg (Fig. 1D);
and (iv) Bpt, an L-transferase that employs Leu-tRNA to conjugate Leu to Nt-Asp, -Glu,
and (possibly) oxidized -Cys (Fig. 1D).
'''

'''
Regulated Proteolysis in Bacteria
Samar A. Mahmoud and Peter Chien
Annual Review of Biochemistry Vol. 87:677-696 (Volume publication date June 2018)
https://doi.org/10.1146/annurev-biochem-062917-012848

In bacteria, regulated proteolysis is carried out by energydependent AAA+ (ATPases associated with cellular activities) proteases that use the power of ATP
hydrolysis to recognize, unfold, translocate, and degrade substrates. Several energy-dependent
proteases exist in bacteria: Lon, ClpXP, ClpAP, ClpCP, ClpEP, HslUV, and FtsH (2–4).
'''

'''
EC Number: 3.4.21.92
a protein + ATP + 2 H2O → 2 a peptide + ADP + phosphate + 3 H+

Component enzyme of ClpAP : ClpP serine protease
Gene:	clpP	Accessions: EG10158 (EcoCyc), b0437, ECK0431

ssrA-mediated protein degradation by ClpXP and ClpAP proteases
'''

'''
https://biocyc.org/gene?orgid=ECOLI&id=G6463-MONOMER
ClpS binds to the ClpA amino-terminus and affects the specificity of protein degradation
by the ClpAP chaperone-protease complex, possibly by interfering with interactions between substrate
and ClpA [Dougan02]. ClpS stimulates ClpAP recognition and degradation of aggregated protein substrates
while it inhibits degradation of non-aggregated substrates, including ClpA [Dougan02].
ClpS also links ClpAP to the N-end rule degradation pathway,
binding to the amino-terminal destabilizing residues in N-end substrates and targeting those substrates
for degradation by ClpAP [Erbse06].
'''

'''
proteins.tsv search for potentially relevant proteins

['protease'containing protein names]
GeneID:  EG10255 	ProtID:  EG10255-MONOMER 	Name:  ecotin monomer; serine protease inhibitor
GeneID:  EG10397 	ProtID:  EG10397-MONOMER 	Name:  intramembrane serine protease GlpG
GeneID:  EG10435 	ProtID:  EG10435-MONOMER 	Name:  regulator of FtsH protease
GeneID:  EG10436 	ProtID:  EG10436-MONOMER 	Name:  regulator of FtsH protease
GeneID:  EG10673 	ProtID:  EG10673-MONOMER 	Name:  outer membrane protease VII (outer membrane protein 3b)
GeneID:  EG10741 	ProtID:  EG10741-MONOMER 	Name:  protease involved in Microcin B17 maturation and in sensitivity to the DNA gyrase inhibitor LetD
GeneID:  EG10760 	ProtID:  EG10760-MONOMER 	Name:  tail-specific protease
GeneID:  EG10786 	ProtID:  EG10786-MONOMER 	Name:  protease III
GeneID:  EG10823 	ProtID:  EG10823-MONOMER 	Name:  DNA strand exchange and recombination protein with protease and nuclease activity
GeneID:  EG10968 	ProtID:  EG10968-MONOMER 	Name:  protease IV
GeneID:  EG11260 	ProtID:  EG11260-MONOMER 	Name:  predicted ATP-dependent protease
GeneID:  EG11506 	ProtID:  EG11506-MONOMER 	Name:  ATP-dependent zinc metalloprotease FtsH
GeneID:  EG11652 	ProtID:  EG11652-MONOMER 	Name:  DegS serine endoprotease
GeneID:  EG11676 	ProtID:  EG11676-MONOMER 	Name:  peptidase component of the HslVU protease
GeneID:  EG11881 	ProtID:  EG11881-MONOMER 	Name:  ATPase component of the HslVU protease
GeneID:  EG12436 	ProtID:  EG12436-MONOMER 	Name:  RseP zinc protease, signal peptide peptidase
GeneID:  G6265 	ProtID:  G6265-MONOMER 	Name:  predicted protease, membrane anchored
GeneID:  G6463 	ProtID:  G6463-MONOMER 	Name:  specificity factor for ClpA-ClpP chaperone-protease complex
GeneID:  G6493 	ProtID:  G6493-MONOMER 	Name:  putative ATP-dependent protease
GeneID:  G7414 	ProtID:  G7414-MONOMER 	Name:  hydrogenase 3 maturation protease
GeneID:  G7653 	ProtID:  G7653-MONOMER 	Name:  predicted protease
GeneID:  G7682 	ProtID:  G7682-MONOMER 	Name:  serine endoprotease, periplasmic
GeneID:  G7689 	ProtID:  G7689-MONOMER 	Name:  protease involved in Microcin B17 maturation and in sensitivity to the DNA gyrase inhibitor LetD

['peptidase' containing protein names]
GeneID:  EG10013 	ProtID:  EG10013-MONOMER 	Name:  murein DD-endopeptidase YebA
GeneID:  EG10201 	ProtID:  EG10201-MONOMER 	Name:  D-alanyl-D-alanine carboxypeptidase, fraction A; penicillin-binding protein 5
GeneID:  EG10202 	ProtID:  EG10202-MONOMER 	Name:  D-alanyl-D-alanine endopeptidase
GeneID:  EG10212 	ProtID:  EG10212-MONOMER 	Name:  dipeptidyl carboxypeptidase II
GeneID:  EG10374 	ProtID:  EG10374-MONOMER 	Name:  &gamma;-glutamyltranspeptidase
GeneID:  EG10530 	ProtID:  EG10530-MONOMER 	Name:  leader peptidase (signal peptidase I)
GeneID:  EG10535 	ProtID:  EG10535-MONOMER 	Name:  Lit, cell death peptidase; phage exclusion; e14 prophage
GeneID:  EG10548 	ProtID:  EG10548-MONOMER 	Name:  prolipoprotein signal peptidase II
GeneID:  EG10570 	ProtID:  EG10570-MONOMER 	Name:  methionine aminopeptidase
GeneID:  EG10694 	ProtID:  EG10694-MONOMER 	Name:  aminopeptidase A/I
GeneID:  EG10695 	ProtID:  EG10695-MONOMER 	Name:  peptidase D
GeneID:  EG10696 	ProtID:  EG10696-MONOMER 	Name:  aminopeptidase N
GeneID:  EG10698 	ProtID:  EG10698-MONOMER 	Name:  proline dipeptidase
GeneID:  EG10956 	ProtID:  EG10956-MONOMER 	Name:  predicted inner membrane peptidase
GeneID:  EG11004 	ProtID:  EG11004-MONOMER 	Name:  oligopeptidase B
GeneID:  EG11253 	ProtID:  EG11253-MONOMER 	Name:  L,D-transpeptidase YcbB
GeneID:  EG11291 	ProtID:  EG11291-MONOMER 	Name:  outer membrane metallopeptidase
GeneID:  EG11359 	ProtID:  EG11359-MONOMER 	Name:  leader peptidase, integral membrane protein
GeneID:  EG11441 	ProtID:  EG11441-MONOMER 	Name:  oligopeptidase A
GeneID:  EG11676 	ProtID:  EG11676-MONOMER 	Name:  peptidase component of the HslVU protease
GeneID:  EG11744 	ProtID:  EG11744-MONOMER 	Name:  putative zinc peptidase
GeneID:  EG11802 	ProtID:  EG11802-MONOMER 	Name:  predicted maturation peptidase for hydrogenase 2
GeneID:  EG11882 	ProtID:  EG11882-MONOMER 	Name:  conserved hypothetical protein of the NlpC/P60 peptidase superfamily
GeneID:  EG11893 	ProtID:  RPOA-MONOMER 	Name:  DD-carboxypeptidase, penicillin-binding protein 6b
GeneID:  EG11920 	ProtID:  EG11920-MONOMER 	Name:  peptidase E, a dipeptidase where amino-terminal residue is aspartate
GeneID:  EG12107 	ProtID:  EG12107-MONOMER 	Name:  prepilin peptidase dependent protein
GeneID:  EG12254 	ProtID:  EG12254-MONOMER 	Name:  predicted zinc-dependent peptidase
GeneID:  EG12310 	ProtID:  EG12310-MONOMER 	Name:  aminopeptidase B
GeneID:  EG12436 	ProtID:  EG12436-MONOMER 	Name:  RseP zinc protease, signal peptide peptidase
GeneID:  EG12557 	ProtID:  SGCB-MONOMER 	Name:  KpLE2 phage-like element; predicted endoglucanase with Zn-dependent exopeptidase domain
GeneID:  EG12867 	ProtID:  EG12867-MONOMER 	Name:  DD-endopeptidase / DD-carboxypeptidase
GeneID:  G6101 	ProtID:  G6101-MONOMER 	Name:  predicted aminopeptidase
GeneID:  G6111 	ProtID:  G6111-MONOMER 	Name:  predicted lipoprotein and C40 family peptidase
GeneID:  G6311 	ProtID:  G6311-MONOMER 	Name:  DLP12 prophage; predicted murein endopeptidase
GeneID:  G6422 	ProtID:  G6422-MONOMER 	Name:  L,D-transpeptidase YbiS
GeneID:  G6571 	ProtID:  G6571-MONOMER 	Name:  L,D-transpeptidase YcfS
GeneID:  G6621 	ProtID:  G6621-MONOMER 	Name:  L,D-carboxypeptidase A
GeneID:  G6686 	ProtID:  G6686-MONOMER 	Name:  Rac prophage; predicted defective peptidase
GeneID:  G6746 	ProtID:  G6746-MONOMER 	Name:  predicted peptidase
GeneID:  G6782 	ProtID:  G6782-MONOMER 	Name:  D-Ala-D-Ala dipeptidase
GeneID:  G6856 	ProtID:  G6856-MONOMER 	Name:  predicted peptidase
GeneID:  G6892 	ProtID:  G6892-MONOMER 	Name:  murein DD-endopeptidase YdhO
GeneID:  G6904 	ProtID:  G6904-MONOMER 	Name:  L,D-transpeptidase YnhG
GeneID:  G7073 	ProtID:  G7073-MONOMER 	Name:  L,D-transpeptidase ErfK
GeneID:  G7118 	ProtID:  G7118-MONOMER 	Name:  predicted peptidase
GeneID:  G7147 	ProtID:  G7147-MONOMER 	Name:  murein LD-carboxypeptidase / murein DD-endopeptidase
GeneID:  G7178 	ProtID:  G7178-MONOMER 	Name:  predicted peptidase
GeneID:  G7247 	ProtID:  G7247-MONOMER 	Name:  broad-specificity examinopeptidase
GeneID:  G7248 	ProtID:  G7248-MONOMER 	Name:  aminopeptidase
GeneID:  G7328 	ProtID:  G7328-MONOMER 	Name:  predicted peptidase
GeneID:  G7491 	ProtID:  G7491-MONOMER 	Name:  predicted peptidase
GeneID:  G7539 	ProtID:  G7539-MONOMER 	Name:  prepilin peptidase
GeneID:  G7652 	ProtID:  G7652-MONOMER 	Name:  predicted peptidase (collagenase-like)

[Clp containing protein names]
GeneID:  EG10156 	ProtID:  EG10156-MONOMER 	Name:  ClpA
GeneID:  EG10157 	ProtID:  EG10157-MONOMER 	Name:  ClpB chaperone
GeneID:  EG10158 	ProtID:  EG10158-MONOMER 	Name:  ClpP
GeneID:  EG10159 	ProtID:  EG10159-MONOMER 	Name:  ClpX
GeneID:  G6463 	ProtID:  G6463-MONOMER 	Name:  specificity factor for ClpA-ClpP chaperone-protease complex
'''

# TODO: incorporate half life flat data for choosing proteins to degrade

def Write_CellProcess(Writer, Comp, ProGen, ProcessID):

    # Molecule indices for Molecular IDs

    # Protease subunit
    Idx_ClpP_SerineProtease = Comp.Master.ID2Idx_Master['EG10158-MONOMER']

    # Substrate specific adaptors
    Idx_ClpA = Comp.Master.ID2Idx_Master['EG10156-MONOMER']   # Hexamer
    Idx_ClpX = Comp.Master.ID2Idx_Master['EG10159-MONOMER']   # Hexamer

    # Idx_ClpB = Comp.Master.ID2Idx_Master['EG10157-MONOMER']
    # Idx_ClpS = Comp.Master.ID2Idx_Master['G6463-MONOMER']


    Idx_H2O = Comp.Master.ID2Idx_Master['WATER[c]']
    Idx_ATP = Comp.Master.ID2Idx_Master['ATP[c]']

    Idx_ADP = Comp.Master.ID2Idx_Master['ADP[c]']
    Idx_Pi = Comp.Master.ID2Idx_Master['PI[c]']
    Idx_Proton = Comp.Master.ID2Idx_Master['PROTON[c]']

    # Set up Idx to Idx system in cell.py

    # Temporary parameters

    # half life ~20 h, which is 2.5 % per hour = 0.025 / 3600
    # half life when starvation induced = 2.5 to 6 % per hour
    Rate_ProteinDegradation = 0.025 / 3600

    # Temporary references
    NUniq_Proteins = Comp.Protein.NUniq_Proteins

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        with Writer.Statement("def Init_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Idx_Proteins", 0)
            Writer.Variable_("self.Idx_RndProteinsDegraded", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Rate_ProteinDegradation", 0)
            Writer.BlankLine()
            Writer.Variable_("self.Count_ProteinsToBeDegraded", 0)
            Writer.Variable_("self.Count_AAsToBeReleased", 0)
            Writer.Variable_("self.Count_AAsToBeReleasedTotal", 0)
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):

            # Local indices
            Writer.VarRange_("self.Idx_Proteins", 0, NUniq_Proteins)

            # Master indices
            Writer.Variable_("self.Cel.Idx_H2O", Idx_H2O)
            Writer.Variable_("self.Cel.Idx_ATP", Idx_ATP)
            Writer.Variable_("self.Cel.Idx_ADP", Idx_ADP)
            Writer.Variable_("self.Cel.Idx_Pi", Idx_Pi)
            Writer.Variable_("self.Cel.Idx_Proton", Idx_Proton)
            Writer.BlankLine()

            Writer.Variable_("self.Rate_ProteinDegradation", Rate_ProteinDegradation)
            # Writer.Variable_("self.Cel.Rate_ProteinDegradation_Matrix", Rate_ProteinDegradation, Shape=[NUniq_mRNAs, NMax_ProteasesPermRNA])
            Writer.BlankLine()


            # Writer.Random()
            # Writer.Statement("random_id = Random()");

            # Writer.InitZeros("self.Count_AAsToBeReplenished", NUniq_AAs, 'int32')

            # Writer.VarFill__("self.Cel.Len_ProteinsDegraded", [NUniq_Proteins, NMax_ProteasesPermRNA], -1)
            # Writer.Overwrite("self.Cel.Len_ProteinsDegradedMax", "self.Cel.Len_Proteins")
            # Writer.BlankLine()
            #
            # Writer.InitZeros("self.Count_ProteinsDegradedNew", NUniq_Proteins, 'int32')
            # Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def ExecuteProcess(self):"):
            Writer.Statement("self.ResetVariables()")
            Writer.Statement("self.Initiation()")
            Writer.Statement("self.Degradation()")
            # Writer.Statement("self.Termination()   # Not Implemented")
            Writer.BlankLine()

        with Writer.Statement("def ResetVariables(self):"):
            Writer.InitZeros("self.Count_ProteinsDegraded", NUniq_Proteins, 'int32')
            Writer.BlankLine()

        # TODO: Model Protease Specificity Interaction
        with Writer.Statement("def ProteaseSubstrateBinding(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def Denaturation(self):"):

            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def Translocation(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def SelectProteinsToDegrade(self):"):
            # Determine the number of proteins to degrade
            Writer.Gather___("Count_Proteins", "self.Cel.Counts", "self.Cel.Idx_Master_Proteins")
            Writer.ReduceSum("Count_ProteinsTotal", "Count_Proteins", 0)
            Writer.Statement("N_ProteinsToBeDegraded = self.DetermineAmountFromRate(Count_ProteinsTotal, self.Rate_ProteinDegradation)")
            Writer.Statement("self.Idx_RndProteinsDegraded = self.PickRandomIndexFromPool_Weighted_Local(N_ProteinsToBeDegraded, self.Idx_Proteins, Count_Proteins)")
            Writer.BlankLine()

        with Writer.Statement("def GetCountOfDegradedProteins(self):"):
            Writer.OnesLike_("OnesForRndProteins", "self.Idx_RndProteinsDegraded", 'int32')
            Writer.ScatNdUpd("self.Count_ProteinsDegraded", "self.Count_ProteinsDegraded", "self.Idx_RndProteinsDegraded", "OnesForRndProteins")
            Writer.BlankLine()

        with Writer.Statement("def Initiation(self):"):
            # Denaturation and translocation steps are the rate-limiting steps in protein degradation
            # These two steps are currently expressed as an average half life value in selecting proteins to degrade
            Writer.Statement("self.ProteaseSubstrateBinding()   # Not implemented yet")
            Writer.Statement("self.Denaturation()   # Not implemented yet")
            Writer.Statement("self.Translocation()   # Not implemented yet")
            Writer.BlankLine()
            Writer.Statement("self.SelectProteinsToDegrade()")
            Writer.Statement("self.GetCountOfDegradedProteins()")
            Writer.BlankLine()

        with Writer.Statement("def ReduceCountOfDegradedProteins(self):"):
            Writer.Statement("self.Count_ProteinsDegraded = self.CorrectCountGettingBelowZeroAfterRemoval(self.Cel.Idx_Master_Proteins, self.Count_ProteinsDegraded)")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Master_Proteins, -self.Count_ProteinsDegraded)")
            Writer.BlankLine()

        with Writer.Statement("def DetermineAAsToBeReleased(self):"):
            # AA Count per protein * count of proteins to be degraded
            Writer.Reshape__("Count_ProteinsDegraded", "self.Count_ProteinsDegraded", [-1, 1])
            Writer.MatrixMul("self.Count_AAsToBeReleased", "self.Cel.Count_AAsInProteins", "Count_ProteinsDegraded")
            Writer.BlankLine()

        with Writer.Statement("def ReplenishAAs(self):"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_AAs, self.Count_AAsToBeReleased)")
            Writer.BlankLine()

        with Writer.Statement("def UpdateByproducts(self):"):
            Writer.ReduceSum("Count_AAsToBeReleasedTotal", "self.Count_AAsToBeReleased", 0)
            Writer.Overwrite("self.Count_AAsToBeReleasedTotal", "Count_AAsToBeReleasedTotal")
            Writer.ReduceSum("Count_RndProteinsDegradedTotal", "self.Count_ProteinsDegraded", 0)
            Writer.BlankLine()

            Writer.Subtract_("N_ProteolysisEvents", "Count_AAsToBeReleasedTotal", "Count_RndProteinsDegradedTotal")
            Writer.Multiply_("N_ProteolysisEventsDouble", "N_ProteolysisEvents", 2)
            Writer.Multiply_("N_ProteolysisEventsTriple", "N_ProteolysisEvents", 3)
            Writer.BlankLine()

            # a protein + ATP + 2 H2O → 2 peptide + ADP + phosphate + 3 Proton
            # TODO: This part may be expressed in a process including a common ATP hydrolysis routine in a fewer lines

            # Consume ATP per proteolysis event
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_ATP, -N_ProteolysisEvents)")

            # Consume two H2O's per proteolysis event
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_H2O, -N_ProteolysisEventsDouble)")

            # Release ADP per proteolysis event
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_ADP, N_ProteolysisEvents)")

            # Release Pi per proteolysis event
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Pi, N_ProteolysisEvents)")

            # Release three protons per proteolysis event
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_Proton, N_ProteolysisEventsTriple)")

            Writer.BlankLine()

        with Writer.Statement("def Degradation(self):"):
            Writer.Statement("self.ReduceCountOfDegradedProteins()")
            Writer.Statement("self.DetermineAAsToBeReleased()")
            Writer.Statement("self.ReplenishAAs()")
            Writer.Statement("self.UpdateByproducts()")
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.PrintStrg("===== Protein Degradation ===== ")
            Writer.PrintStVa("# of Proteins Degraded",
                             "tf.math.reduce_sum(self.Count_ProteinsDegraded)")
            Writer.PrintStVa("# of AA release",
                             "self.Count_AAsToBeReleasedTotal")
            Writer.BlankLine()
