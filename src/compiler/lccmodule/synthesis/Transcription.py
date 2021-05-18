# transcription
# use sigma factor

def Write_Transcription_Init(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Transcription_Init():"):
        Writer.Statement("")

def Write_Transcription_Loop(Writer, CompilerData):
    Writer.BlankLine()
    with Writer.Statement("def Transcription_Init():"):
        Writer.Statement("")

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
            Writer.IDlst2Idx('Idx_NTPs', IDs_NTP, 'MasterID2Index')


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
            self.InitArrayWithZero(VariableName, 0)
            if Weights:
                self.OperGathr("Weights", Weights, Indexes)
                self.Statement("Weights = Weights / len(Weights)")
            else:
                self.InitArrayWithOne("Weights", len(Indexes), Type='int32')
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
