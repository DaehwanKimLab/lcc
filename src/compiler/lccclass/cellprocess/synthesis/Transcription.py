# transcription
# use sigma factor


def Write_CellProcess(Writer, Comp, ProGen, ProcessID):
    # ProGen.GenerateCellProcess(Writer, ProcessID)

    # Debugging switches
    Debugging_TranscriptionInitiation = False
    Debugging_TranscriptElongation = False
    Debugging_TranscriptionTermination = False

    Debugging_TranscriptionInitiation = True
    Debugging_TranscriptElongation = True
    Debugging_TranscriptionTermination = True

    Debugging_Transcription = bool(Debugging_TranscriptionInitiation or Debugging_TranscriptElongation or Debugging_TranscriptionTermination)

    # Molecule indices for Molecular IDs
    Idx_RNAP = Comp.Master.ID2Idx_Master[Comp.Complex.Name2ID_Complexes['RNA polymerase, core enzyme']]
    Count_RNAPFree = int(Comp.Master.Count_Master[Idx_RNAP])

    Idx_NTPs = ProGen.BuildingBlockIdxs('NTP')
    Idx_PPi = Comp.Master.ID2Idx_Master['PPI[c]']

    # Set up Idx to Idx system in cell.py

    # Temporary parameters
    Rate_RNAElongation = 60
    NMax_RNAPsPerGene = 10

    # Temporary references
    NUniq_Genes = Comp.Gene.NUniq_Genes
    NUniq_RNAs = Comp.RNA.NUniq_RNAs

    with Writer.Statement("class F%s(FCellProcess):" % ProcessID):
        ProGen.Init_Common(Writer)

        # Override the abstract method
        with Writer.Statement("def SetUp_ProcessSpecificVariables(self):"):
            Writer.Variable_("self.Cel.Idx_NTPs", Idx_NTPs)
            Writer.Variable_("self.Cel.Idx_PPi", Idx_PPi)
            Writer.BlankLine()

            Writer.Variable_("self.Cel.Rate_RNAElongation", Rate_RNAElongation)
            Writer.Variable_("self.Cel.Rate_RNAElongation_Scalar", Rate_RNAElongation)
            Writer.Variable_("self.Cel.Rate_RNAElongation_Vector", Rate_RNAElongation, NUniq_Genes)
            Writer.VarRepeat("self.Cel.Rate_RNAElongation_Matrix", "self.Cel.Rate_RNAElongation_Vector", [NMax_RNAPsPerGene])
            Writer.Reshape__("self.Cel.Rate_RNAElongation_Matrix", "self.Cel.Rate_RNAElongation_Matrix", [NUniq_Genes, NMax_RNAPsPerGene])
            Writer.BlankLine()

            Writer.Variable_("self.Cel.Count_RNAPFree", Count_RNAPFree)
            Writer.InitZeros("self.Cel.Count_RNAPsBoundPerRNA", NUniq_RNAs)
            Writer.Variable_("self.Cel.Count_ElongationCompletedPerRNA", NUniq_RNAs)
            Writer.BlankLine()

            Writer.InitZeros("self.Cel.Len_RNAsNascentElongated", [NUniq_RNAs, NMax_RNAPsPerGene], 'int32')
            Writer.InitZeros("self.Cel.Len_RNAsNascent", [NUniq_RNAs, NMax_RNAPsPerGene])
            Writer.VarRepeat("self.Cel.Len_RNAsNascentMax", "self.Cel.Len_RNAs", NMax_RNAPsPerGene)
            Writer.Reshape__("self.Cel.Len_RNAsNascentMax", "self.Cel.Len_RNAsNascentMax", [NUniq_RNAs, NMax_RNAPsPerGene])
            Writer.BlankLine()

        # Override the abstract method
        with Writer.Statement("def ExecuteProcess(self):"):
            if Debugging_Transcription:
                Writer.Statement("self.Debugging_Transcription()")
            Writer.Comment__("The following processes run independently, based on the numbers from the previous simulation step")
            Writer.Statement("self.TranscriptionInitiation()")
            Writer.Statement("self.TranscriptElongation()")
            Writer.Statement("self.TranscriptionTermination()")
            Writer.BlankLine()

        with Writer.Statement("def SigmaFactorBinding(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def SelectRNAsToTranscribe(self):"):
            # Replace with weighted random later
            Writer.RndIdxUni("self.Cel.Idx_RndRNAsNascent", "self.Cel.Count_RNAPFree", "self.Cel.Idx_Master_RNAsNascent")

            # Weighted random distribution
            # Writer.OperGathr("self.Cel.Count_GeneCopies", "self.Cel.Counts", "self.Cel.Idx_Master_Genes")
            # Writer.RndIdxWgh("self.Cel.Idx_RndRNAsNascent", "self.Cel.Count_RNAPFree", "self.Cel.Count_GeneCopies")
            Writer.BlankLine()

        with Writer.Statement("def PlaceRNAPsToSelectedNascentRNAs(self):"):
            Writer.InitOnes_("OnesForRndRNAs", "self.Cel.Count_RNAPFree")
            # Distribute RNAP to selected RNA Nascent
            Writer.OperScAdd("self.Cel.DeltaCounts", "self.Cel.Idx_RndRNAsNascent", "OnesForRndRNAs")
            Writer.BlankLine()

        with Writer.Statement("def ClearFreeRNAPs(self):"):
            Writer.ClearVal_("self.Cel.Count_RNAPFree", 'int32')
            Writer.BlankLine()

        with Writer.Statement("def DistributeRNAPToGenes(self):"):
            Writer.Statement("self.SelectRNAsToTranscribe()")
            Writer.Statement("self.PlaceRNAPsToSelectedNascentRNAs()")
            Writer.Statement("self.ClearFreeRNAPs()")
            Writer.BlankLine()

        with Writer.Statement("def TranscriptionInitiation(self):"):
            Writer.Statement("self.SigmaFactorBinding()  # Not implemented")
            Writer.Statement("self.DistributeRNAPToGenes()")
            Writer.BlankLine()




        with Writer.Statement("def GetNascentRNACounts(self):"):
            Writer.OperGathr("self.Cel.Count_RNAsNascent", "self.Cel.Counts", "self.Cel.Idx_Master_RNAsNascent")
            # Temporary code
            Writer.InitOnes_("self.Cel.Count_RNAsNascent", NUniq_RNAs, 'int32')

            Writer.BlankLine()

        with Writer.Statement("def DetermineNascentRNAIndices(self):"):
            # Get indices of Nascent RNA in length table by using tf.range? or repeat (Ones x repeat by Count then turn into boolean)
            # TO DO: use a single API for the routine to fill 1s from the left then zeros

            Writer.Variable_("RNAPArrayForAllTranscripts", 0)
            with Writer.Statement("for i, Count in enumerate(self.Cel.Count_RNAsNascent):"):
                Writer.InitOnes_("RNAP_Present", "Count", 'int32')
                Writer.InitZeros("RNAP_Abscent", "%s - Count" % NMax_RNAPsPerGene, 'int32')
                Writer.OperCncat("RNAPArrayForCurrentTranscript", "RNAP_Present", "RNAP_Abscent")
                with Writer.Statement("if i == 0:"):
                    Writer.Statement("RNAPArrayForAllTranscripts = RNAPArrayForCurrentTranscript")
                    Writer.Statement("continue")
                Writer.OperCncat("RNAPArrayForAllTranscripts", "RNAPArrayForAllTranscripts", "RNAPArrayForCurrentTranscript")
            Writer.Reshape__("RNAPArrayForAllTranscripts", "RNAPArrayForAllTranscripts", [-1, NMax_RNAPsPerGene])

            Writer.OperElGr_("self.Cel.Bool_RNAsNascentElongating", "RNAPArrayForAllTranscripts", 0)
            Writer.Cast_____("self.Cel.Bin_RNAsNascentElongating", "self.Cel.Bool_RNAsNascentElongating", 'int32')
            Writer.OperRdSum("self.Cel.Count_RNAsNascentElongating", "self.Cel.Bin_RNAsNascentElongating")  # For summary, the number of elongating RNAs with RNAP
            Writer.BlankLine()

        with Writer.Statement("def GetNascentRNAsIndicesInLengthMatrix(self):"):
            Writer.Statement("self.GetNascentRNACounts()")
            Writer.Statement("self.DetermineNascentRNAIndices()")
            Writer.BlankLine()

        with Writer.Statement("def ElongateNascentRNAsInLengthMatrix(self):"):
            # Rate Matrix * bin(boolean) + length table (overelongated not counted yet)
            Writer.OperElMul("self.Cel.Rate_RNAElongation_Matrix", "self.Cel.Bin_RNAsNascentElongating", "self.Cel.Rate_RNAElongation")
            Writer.OperElAdd("self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentElongated",
                             "self.Cel.Rate_RNAElongation_Matrix")
            Writer.BlankLine()

        with Writer.Statement("def DetermineOverElongatedNascentRNAs(self):"):
            # compare with len MAX to get boolean matrix for overelongated (>= MAX)
            Writer.OperElGrE("self.Cel.Bool_RNAsNascentOverElongated", "self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentMax")
            Writer.Cast_____("self.Cel.Bin_RNAsNascentOverElongated", "self.Cel.Bool_RNAsNascentOverElongated", 'int32')
            Writer.OperRdSum("self.Cel.Count_RNAsNascentOverElongated", "self.Cel.Bin_RNAsNascentOverElongated")  # For summary, the number of completed/overelongated RNAs
            Writer.BlankLine()

        with Writer.Statement("def DetermineAmountOfOverElongatedNTPs(self):"):
            Writer.OperElMul("MaxLengthForOverElongated", "self.Cel.Bin_RNAsNascentOverElongated", "self.Cel.Len_RNAsNascentMax")
            Writer.OperElMul("LengthForOverElongatedOnly", "self.Cel.Bin_RNAsNascentOverElongated", "self.Cel.Len_RNAsNascentElongated")
            Writer.OperElSub("self.Cel.Len_RNAsNascentOverElongated", "LengthForOverElongatedOnly", "MaxLengthForOverElongated")
            Writer.OperRdSum("self.Cel.Count_NTPsOverElongated", "self.Cel.Len_RNAsNascentOverElongated")  # For summary, over consumped NTPs to be corrected
            Writer.BlankLine()

        with Writer.Statement("def CorrectOverElongationInLengthMatrix(self):"):
            # Update nascent RNA length by removing over elongated length
            Writer.OperElSub("self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentElongated", "self.Cel.Len_RNAsNascentOverElongated")
            # subtract overelongated from the len matrix (now all completed elongation has its max value)
            Writer.BlankLine()

        with Writer.Statement("def CorrectElongatedLengthForTranscript(self):"):
            # Initial elongation rate
            Writer.OperElMul("self.Cel.Rate_RNAElongation_Matrix_Corrected", "self.Cel.Rate_RNAElongation_Matrix", "self.Cel.Bin_RNAsNascentElongating")
            # Adjusted elongation rate (initial - elongated)
            Writer.OperElSub("self.Cel.Rate_RNAElongation_Matrix_Corrected", "self.Cel.Rate_RNAElongation_Matrix_Corrected", "self.Cel.Len_RNAsNascentOverElongated")
            # Calculate the length of elongation per gene
            Writer.OperRdSum("self.Cel.Count_RNAElongationLengthPerGene", "self.Cel.Rate_RNAElongation_Matrix_Corrected", 1)  # For summary, total elongation length per Gene
            # Calculate
            Writer.OperRdSum("self.Cel.Count_RNAElongationLengthTotal", "self.Cel.Rate_RNAElongation_Matrix_Corrected")  # For summary, total elongation length in this sim step
            Writer.BlankLine()

        with Writer.Statement("def ElongateTranscripts(self):"):
            # Nascent RNAs will be elongated.
            Writer.Statement("self.GetNascentRNAsIndicesInLengthMatrix()")
            Writer.Statement("self.ElongateNascentRNAsInLengthMatrix()")
            Writer.Statement("self.DetermineOverElongatedNascentRNAs()")
            Writer.Statement("self.DetermineAmountOfOverElongatedNTPs()")
            Writer.Statement("self.CorrectOverElongationInLengthMatrix()")
            Writer.Statement("self.CorrectElongatedLengthForTranscript()")
            Writer.BlankLine()

        with Writer.Statement("def GetRawNTPConsumption(self):"):
            # Prepare NT elongation length for each gene matrix for multiplication
            Writer.Reshape__("self.Cel.Count_RNAElongationLengthPerGene", "self.Cel.Count_RNAElongationLengthPerGene", [-1, 1])
            Writer.Cast_____("self.Cel.Count_RNAElongationLengthPerGene", "self.Cel.Count_RNAElongationLengthPerGene", 'float32')
            # NTP Frequency per gene * elongating transcript (in NT count) per gene
            Writer.OperMXMul("NTPConsumption_Raw", "self.Cel.Freq_NTsInRNAs", "self.Cel.Count_RNAElongationLengthPerGene")
            Writer.ReturnVar("NTPConsumption_Raw")
            Writer.BlankLine()

        with Writer.Statement("def GetRoundedNTPConsumption(self, NTPConsumption_Raw):"):
            Writer.RoundInt_("NTPConsumption_Rounded", "NTPConsumption_Raw")
            Writer.ReturnVar("NTPConsumption_Rounded")
            Writer.BlankLine()

        with Writer.Statement("def CalculateDiscrepancy(self, NTPConsumption_Rounded):"):
            # Determine the difference between rounded vs corrected elongation
            Writer.OperRdSum("Sum_Rounded", "NTPConsumption_Rounded")
            Writer.Cast_____("Sum_Rounded", "Sum_Rounded", 'int32')
            Writer.OperRdSum("Sum_CorrectedElongationRate", "self.Cel.Rate_RNAElongation_Matrix_Corrected")
            Writer.OperElSub("DeltaSum", "Sum_CorrectedElongationRate", "Sum_Rounded")
            Writer.ReturnVar("DeltaSum")
            Writer.BlankLine()

        with Writer.Statement("def AdjustNTPConsumption(self, NTPConsumption_Rounded, Discrepancy):"):
            # Get equal amount of NTPs if discrepancy is greater than or equal to 4 or less than or eqqal to -4.
            Writer.OperElQuo("N_NTPSets", "Discrepancy", "self.Cel.NUniq_NTPs")
            Writer.VarRepeat("N_NTPSets", "N_NTPSets", "self.Cel.NUniq_NTPs")
            Writer.OperElAdd("NTPConsumption_MissingSet", "NTPConsumption_Rounded", "N_NTPSets")
            Writer.BlankLine()

            # Get random NTP for the remainder (replace with weighted random NTP based on ACGU in all genes)
            Writer.OperElRem("N_NTPRemainder", "Discrepancy", "self.Cel.NUniq_NTPs")
            Writer.Reshape__("N_NTPRemainder", "N_NTPRemainder", -1)
            Writer.InitZeros("NTPConsumption_MissingRemainder", "self.Cel.NUniq_NTPs", 'int32')
            Writer.RndNumUni("Idx_Remainder", "N_NTPRemainder", "0", "self.Cel.NUniq_NTPs")
            Writer.InitOnes_("OnesForRemainder", "N_NTPRemainder", 'int32')
            Writer.OperScAdd("NTPConsumption_MissingRemainder", "Idx_Remainder", "OnesForRemainder")
            Writer.BlankLine()

            # Return the adjusted NTP Consumption
            Writer.OperElAdd("NTPConsumption_Adjusted", "NTPConsumption_MissingSet", "NTPConsumption_MissingRemainder")
            Writer.OperRdSum("TotalNTPConsumption", "NTPConsumption_Adjusted")
            Writer.Statement("tf.debugging.assert_equal(TotalNTPConsumption, self.Cel.Count_RNAElongationLengthTotal)")
            Writer.ReturnVar("NTPConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def DetermineNTPConsumption(self):"):
            Writer.Statement("NTPConsumption_Raw = self.GetRawNTPConsumption()")
            Writer.Statement("NTPConsumption_Rounded = self.GetRoundedNTPConsumption(NTPConsumption_Raw)")
            Writer.Statement("Discrepancy = self.CalculateDiscrepancy(NTPConsumption_Rounded)")
            Writer.Statement("NTPConsumption_Adjusted = self.AdjustNTPConsumption(NTPConsumption_Rounded, Discrepancy)")
            Writer.ReturnVar("NTPConsumption_Adjusted")
            Writer.BlankLine()

        with Writer.Statement("def ConsumeNTPs(self):"):
            Writer.Statement("NTPConsumption_Final = self.DetermineNTPConsumption()")
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_NTPs, NTPConsumption_Final)")
            Writer.BlankLine()

        with Writer.Statement("def ReleasePPi(self):"):
            Writer.Statement("self.AddToDeltaCounts(self.Cel.Idx_PPi, self.Cel.Count_RNAElongationLengthTotal)")
            Writer.BlankLine()

        with Writer.Statement("def TranscriptElongation(self):"):
            Writer.Statement("self.ElongateTranscripts()")
            Writer.Statement("self.ConsumeNTPs()")
            Writer.Statement("self.ReleasePPi()")
            Writer.BlankLine()

        with Writer.Statement("def TranscriptionTermination(self):"):
            # identify all NMAX values by comparison (boolean)

            # add up Trues for RNAP release (RNAP FREE)

            # Reset all NMAX values to zero

            # boolean for reducing RNA nascent (scatter for deltacounts)

            # boolean for increasing RNA (scatter for deltacounts)

            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def Debugging_Transcription(self):"):
            if Debugging_TranscriptionInitiation:
                Writer.Statement("self.Debugging_TranscriptionInitiation()")
            if Debugging_TranscriptElongation:
                Writer.Statement("self.Debugging_TranscriptElongation()")
            if Debugging_TranscriptionTermination:
                Writer.Statement("self.Debugging_TranscriptionTermination()")
            Writer.BlankLine()

        with Writer.Statement("def Debugging_TranscriptionInitiation(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def Debugging_TranscriptElongation(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def Debugging_TranscriptionTermination(self):"):
            Writer.Pass_____()
            Writer.BlankLine()

        with Writer.Statement("def ViewProcessSummary(self):"):
            Writer.Pass_____()
            Writer.BlankLine()


# def SetUpReactions(ProGen):
#     Reactions = list()
#
#     Reaction = None
#
#     # # Reaction No.1
#     #
#     # RXNType = 'Polymerization'
#     # RXNEquation = '4 NTP -> 1 RNA size + 2 PPi'
#     # RXNRate = '60 +- 20 events per second'
#     # RXNTrigger = 'DNA polymerase III, core enzyme >= 4'
#     #
#     # # The final product would be the following:
#     #
#     # # Type = 'Polymerization' or 'Biochemical Reaction'
#     #
#     # # Transcriptional elongation
#     #
#     # RNAP = ProGen.Comp.Complex.Name2ID_Complexes['RNA polymerase, core enzyme']
#     #
#     # for i, ID_RNANascent in enumerate(ProGen.Comp.RNA.ID_RNAsNascent):
#     #
#     #     Type = 'Polymerization'
#     #
#     #     Stoich_MolIDs = ['NTP', ID_RNANascent, 'PPI[c]']
#     #     Stoich_Coeffs = [-1, 1, 2]
#     #
#     #     Rate_Mean = 60  # basepairs per second, accounting for both directions
#     #     Rate_SD = 20
#     #     Rate_UnitTime = 'Second'
#     #
#     #     Trigger_MolIDs = [RNAP]  # 'RNA polymerase, core enzyme'
#     #     Trigger_Thresholds = ['1']  # To be replaced
#     #     Trigger_Conditions = ['>=']  # Greater than or equal to
#     #
#     #     # Generate a reaction dictionary with above inputs
#     #     Reaction['Type'] = Type
#     #     Reaction['Stoichiometry'] = [Stoich_MolIDs, Stoich_Coeffs]
#     #     Reaction['Rate'] = [Rate_Mean, Rate_SD, Rate_UnitTime]
#     #     # Reaction['Trigger'] = [Trigger_MolIDs, Trigger_Thresholds, Trigger_Conditions]
#
#     Reaction_SetUp = Reaction
#     # MolIdxs have been added to Stoichiometry and Trigger variables
#     Reactions.append(Reaction_SetUp)
#
#     return Reactions





'''

                    # RNA selection to be transcribe by RNAP
# [DYNAMIC -> Loop]     # Draw indexes of RNAs to be transcribed

                    # REACTANTS: NTPs (A, C, G, U)
# [STATIC  -> Init]     # Get indexes for NTPs
# [DYNAMIC -> Loop]     # Consume NTPs

                    # REACTANT STOICHIOMETRY:
# [STATIC  -> Init]     # Get frequency of NTPs of the selected RNA

                    # PRODUCTS: Nascent RNA, mature RNA, Pi
# [STATIC  -> Init]     # Get index for Pi
# [DYNAMIC -> Loop]     # Increment 2 Pi's
                        # If new started, 
# [DYNAMIC -> Loop]         # increment Nascent RNA Count
                        # If completed,
# [DYNAMIC -> Loop]         # Increment RNA Count
# [DYNAMIC -> Loop]         # Consume nascent RNA Count

                    # Rate:
# [STATIC  -> Init]     # A fixed elongation Rate by RNAP, OR
# [STATIC  -> Init]     # Known transcription efficiency for each RNA applied

'''
'''
        # Select RNAs to be transcribe by RNAP
            # Which RNAs? --> List of Genes (indices to draft from)
            # Relative Abundance of Gene copy number --> DNA replication
            # How many RNAs to select? --> Number of RNAP

                                           inactive vs active (free vs engaged)
            Writer.SelectIdx('SelectedIdxRNA', 'N_RNAP', 'Idx_RNA', Weights='Counts_gene')
            # def SelectIdx(self, VariableName, N_MoleculesToDistribute, Indexes, Weights='None'):

        # REACTANTS:
            # NTPs (A, C, G, U)
            IDs_NTP = ['ATP', 'CTP', 'GTP', 'UTP'] # or = CompilerData.NTPs

            # Get Indexes for NTPs
            # def IDlst2Idx(self, MolIndexList, MolList, Mol2Index):
            Writer.IDlst2Idx('Idx_NTPs', IDs_NTP, 'ID2Idx_Master')


            # Consume NTPs
            # def OperGathr(self, VariableName, Source, Index):
            Writer.OperGathr("Counts_NTPs", "Counts_Master", "Idx_NTPs")


        def OperGathr(self, VariableName, Target, Index):
            self.PrepGathr(Target, Index)
            self.Statement("%s = tf.gather(%s, %s)" % (VariableName, Target, Index))

        def PrepGathr(self, Target, Index):
            self.Reshape__(Target, Target, -1)
            self.Reshape__(Index, Index, [-1, 1])

        def Reshape__(self, DestVar, SrcVar, Shape):
            Line = '%s = tf.reshape(%s, %s)' % (DestVar, SrcVar, str(Shape))
            self.Statement(Line)
            self.DebugPVar(DestVar)

        # REACTANT STOICHIOMETRY:
            # Frequency of A, C, G, U of the selected RNA



        # PRODUCTS:
            # Nascent RNA (RNAP processed length on RNA)
            # RNA Count if completed
            # 2 Pi's

        Writer.Increment(TargetMX, N_MoleculesToDistribute, Indexes, IncrementValue, WeightMX='None'):

        def Increment(self, TargetMX, N_MoleculesToDistribute, Indexes, IncrementValue, WeightMX='None'):
            self.GenValues('Values', N_MoleculesToDistribute, Indexes, Weights=WeightMX)
            self.DebugPStV("Before Increment:", TargetMX)
            with self.Statement("for i in range(%s):" % N_MoleculesToDistribute):
                self.Statement("j = Values[i]")
                self.Statement("SelectedFromIndex = %s[j]" % Indexes)
                self.Statement("SelectedFromIndex = tf.reshape(SelectedFromIndex, [-1, 1])")
                self.Statement("TargetMXDataType = tf.shape(%s).dtype" % TargetMX)
                self.Statement("One = tf.ones(1, TargetMXDataType)")
                self.Statement("%s = tf.tensor_scatter_nd_add(%s, SelectedFromIndex, One * %s)" % (
                    TargetMX, TargetMX, IncrementValue))
            self.DebugPStV("After Increment:", TargetMX)

        def GenValues(self, VariableName, N_MoleculesToDistribute, Indexes, Weights='None'):
            self.InitZeros(VariableName, 0)
            if Weights:
                self.OperGathr("Weights", Weights, Indexes)
                self.Statement("Weights = Weights / len(Weights)")
            else:
                self.InitZeros("Weights", len(Indexes), Type='int32')
            with self.Statement("for i in range(%s):" % N_MoleculesToDistribute):
                self.Statement("Values = tf.data.experimental.sample_from_datasets(%s, weights=Weights)" % Indexes)
                self.Statement("%s = tf.concat([%s, Value], 0)" % (VariableName, VariableName))
            self.DebugPStV("Values Generated:", VariableName)
            return VariableName

        # PRODUCT STOICHIOMETRY:
        # 1
        # 1 if completed
        # 2

        # Rate:
        # Arbitrarilly set or get a table for each RNA'''
