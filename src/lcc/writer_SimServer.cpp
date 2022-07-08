#include "writer.h"

using namespace std;

void FWriter::GenerateVisObjects(std::ofstream& ofs, int indents, std::string ObjectFamilyName, std::string N_Objects_Str) {
    std::string in = "    ";
    std::string IN = std::string(4 * indents, ' ');

    std::string Idx = ObjectFamilyName + "_i";

    ofs << IN+ "# For " << ObjectFamilyName << " VisObjects" << endl;
    ofs << IN+ "VisObjects_" << ObjectFamilyName << " = {}" << endl;
    ofs << IN+ "for " << Idx << " in range(" << N_Objects_Str << "):" << endl;
    ofs << IN+ in+ "VisObjID = " << Idx << " + 1   # ID 0 is used for static objects" << endl;
    ofs << IN+ in+ "VisObjects_" << ObjectFamilyName << "[VisObjID] = lccsimulation_pb2.MVisObjectData(" << endl;
    ofs << IN+ in+ in+ "ID=VisObjID," << endl;
    ofs << IN+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_" << Utils::UpperCaseStr(ObjectFamilyName) << ", " << endl;
    ofs << IN+ in+ in+ "Position=ZeroVec,   # TODO: Requires a check with visualization" << endl;
    ofs << IN+ in+ in+ "Rotation=ZeroVec,   # TODO: Requires a check with visualization" << endl;
    ofs << IN+ in+ in+ "# Scale =, # Scale doesn't change, leave it out" << endl;
    ofs << IN+ in+ in+ "# Color =, # Color doesn't change, leave it out" << endl;
    ofs << IN+ in+ ")" << endl;
    ofs << endl;
}


void FWriter::SimServer() {
    std::cout << "Generating simulation server..." << std::endl;

    // write SimServer.py
    std::ofstream ofs(Option.SimServerFile.c_str());
    std::string endl = "\n";

    ofs << "# Compile 'Models/CentralDogma/4_Replication_Transcription_Translation.lpp'" << endl;
    ofs << "# Use option '--maxgenes <number> to control the maximum number of genes imported from the genome'" << endl;
    ofs << endl;
    ofs << "## lcc" << endl;
    ofs << "import random" << endl;
    ofs << "import SimModule" << endl;
    ofs << "import SimFunctions as SimF" << endl;
    //ofs << "import SimVis2D as SimV" << endl;
    ofs << "import SimState as SimS" << endl;
    ofs << "import numpy as np" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import time" << endl;
    ofs << endl;
    ofs << "## gRPC" << endl;
    ofs << "import asyncio" << endl;
    ofs << "from concurrent import futures" << endl;
    ofs << "import logging" << endl;
    ofs << endl;
    ofs << "import grpc" << endl;
    ofs << endl;
    ofs << "import os, sys" << endl;
    ofs << "os.system('python -m grpc_tools.protoc -I. --python_out=. --grpc_python_out=. ./protos/lccsimulation.proto')" << endl;
    ofs << "sys.path.insert(0, './protos')" << endl;
    ofs << "from protos import lccsimulation_pb2" << endl;
    ofs << "from protos import lccsimulation_pb2_grpc" << endl;
    ofs << endl;

    ofs << "def RandomInt(maxInt):" << endl;
    ofs << in+ "return np.random.randint(0, high=maxInt)" << endl;
    ofs << endl;

    ofs << "def GenRandomColor():" << endl;
    ofs << in+ "return lccsimulation_pb2.MVector3(X=RandomInt(256),Y=RandomInt(256),Z=RandomInt(256))" << endl;
    ofs << endl;

    ofs << "def GenGreenToRedColor(ValueBTWZeroAndOne):" << endl;
    ofs << in+ "return lccsimulation_pb2.MVector3(X = int(255 * ValueBTWZeroAndOne), Y = int(255 * (1 - ValueBTWZeroAndOne)), Z = 0)" << endl;
    ofs << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto Organisms = Context.GetUniqueContainers_LocationList("Organism");

    auto Chromosomes = Context.GetSubList_MoleculeList("Chromosome");
    auto DNAPs = Context.GetSubList_MoleculeList("DNAP");
    auto RNAPs = Context.GetSubList_MoleculeList("RNAP");
    auto Ribosomes = Context.GetSubList_MoleculeList("Ribosome");

    ofs << endl;
    ofs << "class LCCSimulation(lccsimulation_pb2_grpc.LCCSimulationServicer):" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.IsRunning = False" << endl;
    ofs << endl;
    ofs << in+ in+ "self.State = None" << endl;
    ofs << in+ in+ "self.Data = None" << endl;
    ofs << in+ in+ "self.DataManager = None" << endl;
    ofs << in+ in+ "self.SimM = None" << endl;
    ofs << endl;
    for (auto organism : Organisms) {
        ofs << in+ in+ "self." << organism->Name << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ "def Initialize(self, request, context):" << endl;
    ofs << in+ in+ "print('initializing')" << endl;
    ofs << in+ in+ "# Load model" << endl;
    ofs << in+ in+ "self.State = SimS.FState()" << endl;
    ofs << in+ in+ "self.Data = SimS.FDataset()" << endl;
    ofs << in+ in+ "self.DataManager = SimS.FDataManager()" << endl;
    ofs << in+ in+ "self.SimM = SimModule.FSimulation(self.State, self.Data, self.DataManager)" << endl;
    ofs << endl;
    ofs << in+ in+ "# Initialize model" << endl;
    ofs << in+ in+ "self.SimM.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ "SimUnitTime = 0.2" << endl;

    if (!Context.ThresholdList.empty()) {
        if (Context.ThresholdList[0].first == "Am") {
            ofs << in+ in+ "self.SimM.Receptivity(20000)   # Pass time" << endl;
            ofs << endl;
        }
    }

    ofs << in+ in+ "# Setup Static Objects" << endl;
    
    ofs << in+ in+ "InitVisObjects = []" << endl;
    ofs << in+ in+ "ZeroVec = lccsimulation_pb2.MVector3(X=0, Y=0, Z=0)" << endl;
    ofs << in+ in+ "UnitVec = lccsimulation_pb2.MVector3(X=1, Y=1, Z=1)" << endl;
    ofs << in+ in+ "White = lccsimulation_pb2.MVector3(X=255, Y=255, Z=255)" << endl;
    ofs << in+ in+ "Green = lccsimulation_pb2.MVector3(X=0, Y=255, Z=50)" << endl;
    ofs << in+ in+ "Red = lccsimulation_pb2.MVector3(X=255, Y=0, Z=50)" << endl;
    ofs << in+ in+ "Blue = lccsimulation_pb2.MVector3(X=143, Y=186, Z=255)" << endl;
    ofs << in+ in+ "Yellow = lccsimulation_pb2.MVector3(X=255, Y=255, Z=102)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Static Petri Dish" << endl;
    ofs << in+ in+ "InitVisObjects.append(lccsimulation_pb2.MVisObjectData(" << endl;
    ofs << in+ in+ in+ in+ "ID=0, " << endl;
    ofs << in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_PETRI_DISH," << endl;
    ofs << in+ in+ in+ in+ "Position=lccsimulation_pb2.MVector3(X=100, Y=100, Z=0)," << endl;
    ofs << in+ in+ in+ in+ "Rotation=ZeroVec," << endl;
    ofs << in+ in+ in+ in+ "Scale=lccsimulation_pb2.MVector3(X=self.SimM.State.Dimension_X, Y=self.SimM.State.Dimension_Y, Z=1)," << endl;
    ofs << in+ in+ in+ in+ "Color=Yellow))" << endl;
    ofs << endl;

    // temporary static molecule 'distribution'
// TODO: take dynamic distribution


    if (!MolLoc.empty()) {
        ofs << in+ in+ "# Static Molecule Distribution Data" << endl;

        for (int i = 0; i < MolLoc.size(); i++) {
            ofs << in+ in+ "# " << MolLoc[i]->Name << endl;
            ofs << in+ in+ "for Pos in SimF.InitializeStaticParticles(" << MolLoc[i]->Coord[0] << ", " << MolLoc[i]->Coord[1] << ", Particle_N=5000, Particle_PerLayer=20, Particle_SpreadFactor=1.2):" << endl;
            ofs << in+ in+ in+ "NewObj = lccsimulation_pb2.MVisObjectData(" << endl;
            ofs << in+ in+ in+ in+ "ID=" << i + 1 << ", " << endl;
            ofs << in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_GLUCOSE,   # same as M_SPHERE?" << endl;
            ofs << in+ in+ in+ in+ "Position=lccsimulation_pb2.MVector3(X=Pos[0], Y=Pos[1], Z=0)," << endl;
            ofs << in+ in+ in+ in+ "Rotation=ZeroVec," << endl;
            ofs << in+ in+ in+ in+ "Scale=lccsimulation_pb2.MVector3(X=3,Y=3,Z=3)," << endl;
            ofs << in+ in+ in+ in+ "Color=Blue)" << endl;
            ofs << in+ in+ in+ "" << endl;
            ofs << in+ in+ in+ "InitVisObjects.append(NewObj)" << endl;
            ofs << in+ in+ "" << endl;
        }
    }

    ofs << in+ in+ "# Setup Dynamic Objects" << endl;

    if (!Organisms.empty()) {
        
        ofs << in+ in+ "# Dynamic Organisms" << endl;

        for (auto& organism : Organisms) {
            auto Organism = dynamic_cast<FOrganism*>(organism);
            ofs << in+ in+ "X, Y, Angle = self.SimM.GetPositionXYAngleByName('" << Organism->Name << "')" << endl;
             ofs << in+ in+ "for i in range(len(X)):" << endl;
            ofs << in+ in+ in+ "PosVec = lccsimulation_pb2.MVector3(X=X[i], Y=Y[i], Z=0)" << endl;
            ofs << in+ in+ in+ "RotVec = lccsimulation_pb2.MVector3(X=0, Y=Angle[i] * (180/np.pi), Z=0)" << endl;
            ofs << in+ in+ in+ "" << endl;
            ofs << in+ in+ in+ "InitVisObjects.append(lccsimulation_pb2.MVisObjectData(" << endl;
            ofs << in+ in+ in+ in+ "ID = i + 1," << endl;
            ofs << in+ in+ in+ in+ "ObjType = lccsimulation_pb2.VisObjectType.M_" << Utils::UpperCaseStr(Organism->Species) << "," << endl;
            ofs << in+ in+ in+ in+ "Position = PosVec," << endl;
            ofs << in+ in+ in+ in+ "Rotation = RotVec," << endl;
            ofs << in+ in+ in+ in+ "Scale = lccsimulation_pb2.MVector3(X=0.2,Y=0.2,Z=0.2)," << endl;
            ofs << in+ in+ in+ in+ "Color = Green," << endl;
            ofs << in+ in+ in+ "))" << endl;
            ofs << endl;
        }
    }

    ofs << in+ in+ "# DNA Init Data" << endl;
    ofs << in+ in+ "DNA_Init = []" << endl;
    ofs << endl;

    if (!Chromosomes.empty()) {

        ofs << in+ in+ "# DNA Annotation" << endl;
        ofs << endl;

        ofs << in+ in+ "# Temporary halving" << endl;
        ofs << in+ in+ "Sequence=self.SimM.State.OpenFASTADatabase(r'./Database/EscherichiaColi.fasta')[0]" << endl;
        ofs << in+ in+ "UniprotLinkTag='https://www.uniprot.org/uniprotkb/'" << endl;
        ofs << endl;

        ofs << in+ in+ "DNA_Annotations = lccsimulation_pb2.MDNA_AnnotationData(" << endl;
        //ofs << in+ in+ in+ "Sequence='ACGT'," << endl;
        //ofs << in+ in+ in+ "Sequence=self.SimM.State.OpenFASTADatabase(r'./Database/EscherichiaColi.fasta')," << endl;
        ofs << in+ in+ in+ "Sequence=Sequence," << endl;
        ofs << in+ in+ in+ "Gene_StartIndex_bp=self.State.Pos_Gene_Start_bp[0]," << endl;
        ofs << in+ in+ in+ "Gene_EndIndex_bp=self.State.Pos_Gene_End_bp[0]," << endl;
        ofs << in+ in+ in+ "Gene_Symbol=self.State.Name_Genes," << endl;
        ofs << in+ in+ in+ "Gene_PDBID=self.State.Name_PDBIDs," << endl;
        ofs << in+ in+ in+ "Gene_UniprotLink=[UniprotLinkTag + i for i in self.State.Name_UniprotIDs]," << endl;
        ofs << in+ in+ ")" << endl;
        ofs << endl;
        
        ofs << in+ in+ "# DNA Position" << endl;
        ofs << endl;

        // temporary code for visualization
        ofs << in+ in+ "DNAScale = 7" << endl; 
        ofs << endl;
        ofs << in+ in+ "# Center the DNA" << endl;
        ofs << in+ in+ "pos_avg = self.State.Pos_Ref[0]" << endl;
        ofs << in+ in+ "for pos in self.State.Pos_Ref:" << endl;
        ofs << in+ in+ in+ "pos_avg += pos" << endl;
        ofs << endl;
        ofs << in+ in+ "pos_avg /= len(self.State.Pos_Ref)" << endl;
        ofs << endl;
        // 
        
        ofs << in+ in+ "Positions = []" << endl;
        ofs << in+ in+ "for pos in self.State.Pos_Ref:" << endl;
        ofs << in+ in+ in+ "pos_temp = DNAScale * (pos - pos_avg)" << endl; // temporary code for visualization
        ofs << in+ in+ in+ "Positions.append(lccsimulation_pb2.MVector3(X=pos_temp[0], Y=pos_temp[1], Z=pos_temp[2]))" << endl;
        ofs << endl;
        ofs << in+ in+ "DNA_Positions = lccsimulation_pb2.MDNA_PositionData(" << endl;
        ofs << in+ in+ in+ "Points=Positions" << endl;
        ofs << in+ in+ ")" << endl;
        ofs << endl;
        
        ofs << in+ in+ "DNA_Init.append(lccsimulation_pb2.MDNA_InitData(" << endl;
        ofs << in+ in+ in+ "DNA_Annotations=DNA_Annotations," << endl;
        ofs << in+ in+ in+ "DNA_Positions=DNA_Positions," << endl;
        ofs << in+ in+ "))" << endl;
    }

    ofs << in+ in+ "# Name Init Data" << endl;
    ofs << in+ in+ "Name_Init = []" << endl;
    ofs << endl;

    ofs << in+ in+ "Name_Init.append(lccsimulation_pb2.MNames(" << endl;
    ofs << in+ in+ in+ "Name_Count_All=self.State.Mol_Names," << endl;
    ofs << in+ in+ "))" << endl;
    ofs << endl;

    ofs << in+ in+ "# Idx Init Data" << endl;
    ofs << in+ in+ "Idx_Init = []" << endl;
    ofs << endl;

    ofs << in + in + "# Organize Indices" << endl;
    ofs << endl;

    ofs << in+ in+ "Idx_Genes = self.State.Idx_Template_Transcription" << endl;
    ofs << in+ in+ "Idx_RNAs = self.State.Idx_Target_Transcription" << endl;
    ofs << in+ in+ "Idx_mRNAs = self.State.Idx_Template_Translation" << endl;
    ofs << in+ in+ "Idx_Proteins = self.State.Idx_Target_Translation" << endl;
    ofs << endl;

    ofs << in+ in+ "Dict_mRNA2Protein = {}" << endl;
    ofs << in+ in+ "if np.any(Idx_Proteins):" << endl;
    ofs << in+ in+ in+ "for i in range(Idx_Proteins.shape[0]):" << endl;
    ofs << in+ in+ in+ in+ "Dict_mRNA2Protein[Idx_mRNAs[i]] = Idx_Proteins[i]" << endl;
    ofs << endl;

    ofs << in+ in+ "Idx_Gene2RNA = {}" << endl;
    ofs << in+ in+ "Idx_RNA2Protein = {}" << endl;
    ofs << in+ in+ "Idx_Gene2Protein = {}" << endl;
    ofs << in+ in+ "Idx_ProteinOrmRNA2GeneOrRNA_Local = {}" << endl;
    ofs << endl;
    ofs << in+ in+ "if np.any(Idx_Genes):" << endl;
    ofs << in+ in+ in+ "nth_protein_or_mRNA = 0" << endl;
    ofs << in+ in+ in+ "for i in range(Idx_Genes.shape[0]):" << endl;
    ofs << in+ in+ in+ in+ "idx_gene = Idx_Genes[i]" << endl;
    ofs << in+ in+ in+ in + "idx_rna = Idx_RNAs[i]" << endl;
    ofs << in+ in+ in+ in + "Idx_Gene2RNA[idx_gene] = idx_rna" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "idx_protein = -1   # -1 for non-coding genes (only generates non-coding type RNA (not mRNA), hence no protein)" << endl;
    ofs << in+ in+ in+ in+ "if idx_rna in Idx_mRNAs:" << endl;
    ofs << in+ in+ in+ in+ in+ "idx_protein = Dict_mRNA2Protein[idx_rna]" << endl;
    ofs << in+ in+ in+ in+ in+ "Idx_ProteinOrmRNA2GeneOrRNA_Local[nth_protein_or_mRNA] = i" << endl;
    ofs << in+ in+ in+ in+ in+ "nth_protein_or_mRNA += 1" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "Idx_RNA2Protein[idx_rna] = idx_protein" << endl;
    ofs << in+ in+ in+ in+ "Idx_Gene2Protein[idx_gene] = idx_protein" << endl;
    ofs << endl;

    ofs << in+ in+ "Idx_Init.append(lccsimulation_pb2.MIdx(" << endl;
    ofs << in+ in+ in+ "Gene=Idx_Genes," << endl;
    ofs << in+ in+ in+ "RNA=Idx_RNAs," << endl;
    ofs << in+ in+ in+ "mRNA=Idx_mRNAs," << endl;
    ofs << in+ in+ in+ "Protein=Idx_Proteins," << endl;
    ofs << in+ in+ in+ "Gene2RNA=Idx_Gene2RNA," << endl;
    ofs << in+ in+ in+ "RNA2Protein=Idx_RNA2Protein," << endl;
    ofs << in+ in+ in+ "Gene2Protein=Idx_Gene2Protein," << endl;
    ofs << in+ in+ in+ "ProteinOrmRNA2GeneOrRNA_Local=Idx_ProteinOrmRNA2GeneOrRNA_Local," << endl;
    ofs << in+ in+ "))" << endl;
    ofs << endl;

    ofs << in+ in+ "return lccsimulation_pb2.MInitData(InitObjects=InitVisObjects, InitDNA=DNA_Init, InitName=Name_Init, InitIdx=Idx_Init)" << endl;
    ofs << endl;

    // TODO: stream run
    ofs << in+ "def Run(self, request, context):" << endl;
    ofs << in+ in+ "print('running')" << endl;
    ofs << in+ in+ "self.IsRunning = True" << endl;
    ofs << endl;

    ofs << in+ in+ "# Setup Static Objects" << endl;

    ofs << in+ in+ "InitVisObjects = []" << endl;
    ofs << in+ in+ "ZeroVec = lccsimulation_pb2.MVector3(X=0, Y=0, Z=0)" << endl;
    ofs << in+ in+ "UnitVec = lccsimulation_pb2.MVector3(X=1, Y=1, Z=1)" << endl;
    ofs << endl;



    //ofs << in+ in+ "# Setup Central Dogma Numpy Arrays" << endl;

    //ofs << in+ in+ "if np.any(self.State.Pos_Pol_Replication):" << endl;
    //ofs << in+ in+ in+ "Pos_Pol_Replication = self.State.Pos_Pol_Replication" << endl;
    //ofs << in+ in+ "if np.any(self.State.Pos_Pol_Transcription):" << endl;
    //ofs << in+ in+ in+ "Pos_Pol_Transcription = self.State.Pos_Pol_Transcription" << endl;
    //ofs << in+ in+ "if np.any(self.State.Pos_Pol_Translation):" << endl;
    //ofs << in+ in+ in+ "Pos_Pol_Translation = self.State.Pos_Pol_Translation" << endl;
    //ofs << in+ in+ "if np.any(self.State.Dir_Pol_Replication):" << endl;
    //ofs << in+ in+ in+ "Dir_Pol_Replication = self.State.Dir_Pol_Replication" << endl;
    //ofs << in+ in+ "if np.any(self.State.Dir_Pol_Transcription):" << endl;
    //ofs << in+ in+ in+ "Dir_Pol_Transcription = self.State.Dir_Pol_Transcription" << endl;
    //ofs << endl;


    ofs << in+ in+ "def response_messages():" << endl;
    ofs << endl;
    
    ofs << in+ in+ in+ "SimUnitTime = 0.2" << endl;
    ofs << in+ in+ in+ "MinElapsedTime = 0.2" << endl;
    ofs << in+ in+ in+ "ElapsedTime = 0" << endl;
    ofs << in+ in+ in+ "PrevTime = datetime.now()" << endl;
    ofs << in+ in+ in+ "SimStep = 0" << endl;
    ofs << endl;
    
    ofs << in+ in+ in+ "while(True):" << endl;
    ofs << in+ in+ in+ in+ "if not self.IsRunning:" << endl;
    ofs << in+ in+ in+ in+ in+ "break" << endl;
    ofs << in+ in+ in+ in+ "start = time.time()" << endl;
    ofs << endl;
    
    ofs << in+ in+ in+ in+ "CurrTime = datetime.now()" << endl;
    ofs << in+ in+ in+ in+ "ElapsedTime += (CurrTime - PrevTime).total_seconds()" << endl;
    ofs << in+ in+ in+ in+ "PrevTime = CurrTime" << endl;
    ofs << endl;
    
    ofs << in+ in+ in+ in+ "while ElapsedTime >= SimUnitTime:" << endl;
    ofs << in+ in+ in+ in+ in+ "self.SimM.Receptivity(20)" << endl;
    ofs << in+ in+ in+ in+ in+ "ElapsedTime -= SimUnitTime" << endl;
    ofs << endl;
    
    ofs << in+ in+ in+ in+ "CurrentSimStep = self.SimM.GetSimStep()" << endl;
    ofs << in+ in+ in+ in+ "if SimStep == CurrentSimStep:" << endl;
    ofs << in+ in+ in+ in+ in+ "continue" << endl;
    ofs << in+ in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ in+ "SimStep = self.SimM.GetSimStep()" << endl;
    ofs << endl;    
    
    ofs << in+ in+ in+ in+ "# Record data SimUnitState" << endl;
    
    ofs << in+ in+ in+ in+ "# Temporary: The following state is at the organism level (highest container level)" << endl;

    ofs << in+ in+ in+ in+ "VisObjects = {} # map from id --> VisObjectData" << endl;
    ofs << in+ in+ in+ in+ "Counts = {} # map from id --> VisObjectData" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Setup Default Numpy Arrays" << endl;

    ofs << in+ in+ in+ in+ "ZeroArray = np.zeros([self.State.Count_All.shape[0], 1])" << endl;
    ofs << in+ in+ in+ in+ "NegOneArray = np.full([self.State.Count_All.shape[0], 1], -1)" << endl;

    ofs << in+ in+ in+ in+ "Pos_Pol_Replication = NegOneArray" << endl;
    ofs << in+ in+ in+ in+ "Pos_Pol_Transcription = NegOneArray" << endl;
    ofs << in+ in+ in+ in+ "Pos_Pol_Translation = NegOneArray" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "Dir_Pol_Replication = ZeroArray" << endl;
    ofs << in+ in+ in+ in+ "Dir_Pol_Transcription = ZeroArray" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "Pos_Pol_Template_Transcription = NegOneArray" << endl;
    ofs << in+ in+ in+ in+ "Pos_Pol_Template_Translation = NegOneArray" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "Count_Nascent_Chromosome = ZeroArray" << endl;
    ofs << in+ in+ in+ in+ "Count_Nascent_Gene = ZeroArray" << endl;
    ofs << in+ in+ in+ in+ "Count_Nascent_RNA = ZeroArray" << endl;
    ofs << in+ in+ in+ in+ "Count_Nascent_Protein = ZeroArray" << endl;
    ofs << endl;


    std::vector<std::string> VisObjectFamilyListInOrganism, Processes;
    std::string VisObjectFamily, Process;
    std::string DNAP, RNAP, Ribosome;

    DNAP = "DNAP";
    RNAP = "RNAP";
    Ribosome = "Ribosome";

    VisObjectFamilyListInOrganism.push_back(DNAP);
    VisObjectFamilyListInOrganism.push_back(RNAP);
    VisObjectFamilyListInOrganism.push_back(Ribosome);

    Processes = { "Replication", "Transcription", "Translation" };

    for (auto& process : Processes) {
        ofs << in+ in+ in+ in+ process + "s = {} # map from id --> VisObjectData" << endl;
    }
    ofs << endl;

    if (!Organisms.empty()) {
        for (auto& organism : Organisms) {
            auto Organism = dynamic_cast<FOrganism*>(organism);

            // For an organism position

            ofs << in+ in+ in+ in+ "# " << Organism->Name << " State" << endl;
            ofs << in+ in+ in+ in+ "ReplicationCompletionRate = self.SimM.GetReplicationCompletionRateByCompartmentName('" << Organism->Name << "')" << endl;

            ofs << in+ in+ in+ in+ "X, Y, Angle = self.SimM.GetPositionXYAngleByName('" << organism->Name << "')" << endl;
            ofs << in+ in+ in+ in+ "for i in range(len(X)):   # i here is an organism ID" << endl;
            ofs << in+ in+ in+ in+ in+ "PosVec = lccsimulation_pb2.MVector3(X=X[i], Y=Y[i], Z=0)" << endl;
            ofs << in+ in+ in+ in+ in+ "RotVec = lccsimulation_pb2.MVector3(X=0, Y=Angle[i] * (180/np.pi), Z=0)" << endl;
            ofs << in+ in+ in+ in+ in+ "GrowthVec = lccsimulation_pb2.MVector3(X=0.2 * (1 + ReplicationCompletionRate[i]), Y=0.2, Z=0.2)" << endl;
            ofs << in+ in+ in+ in+ in+ "ColorVec = GenGreenToRedColor(ReplicationCompletionRate[i])" << endl;
            ofs << endl;
            ofs << in+ in+ in+ in+ in+ "ObjID = i + 1   # ID 0 is used for static objects" << endl;
            ofs << in+ in+ in+ in+ in+ "VisObjects[ObjID] = lccsimulation_pb2.MVisObjectData(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_" << Utils::UpperCaseStr(Organism->Species) << "," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Position=PosVec," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Rotation=RotVec," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Scale=GrowthVec," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Color=ColorVec," << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;
            ofs << endl;

            // Count_All

            ofs << in+ in+ in+ in+ in+ "Counts[ObjID] = lccsimulation_pb2.MState_Count(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Count_All=self.State.Count_All[i]," << endl;

            std::string Text_Count_Nascent_Chromosome = "self.State.Count_Nascent_Chromosome";
            std::string Text_Count_Nascent_Gene = "self.State.Count_Nascent_Gene";
            std::string Text_Count_Nascent_RNA = "self.State.Count_Nascent_RNA";
            std::string Text_Count_Nascent_Protein = "self.State.Count_Nascent_Protein";
            if (DNAPs.empty()) {
                Text_Count_Nascent_Chromosome = "Count_Nascent_Chromosome";
                Text_Count_Nascent_Gene = "Count_Nascent_Gene";
            }
            if (RNAPs.empty()) {
                Text_Count_Nascent_Gene = "Count_Nascent_Gene";
                Text_Count_Nascent_RNA = "Count_Nascent_RNA";
            }
            if (Ribosomes.empty()) {
                Text_Count_Nascent_Protein = "Count_Nascent_Protein";
            }

            ofs << in+ in+ in+ in+ in+ in+ "Count_Nascent_Chromosome=" << Text_Count_Nascent_Chromosome << "[i]," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Count_Nascent_Gene=" << Text_Count_Nascent_Gene << "[i], " << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Count_Nascent_RNA=" << Text_Count_Nascent_RNA << "[i]," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Count_Nascent_Protein=" << Text_Count_Nascent_Protein << "[i]," << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;
            ofs << endl;
            
            // Organize Central dogma info in the organism

            for (int i = 0; i < VisObjectFamilyListInOrganism.size(); i++) {
                std::string Text_Pos_Pol = "self.State.Pos_Pol_";
                if ((VisObjectFamilyListInOrganism[i] == "DNAP") & (DNAPs.empty())) {
                    Text_Pos_Pol = "Pos_Pol_";
                } 
                else if ((VisObjectFamilyListInOrganism[i] == "RNAP") & (RNAPs.empty())) {
                    Text_Pos_Pol = "Pos_Pol_";
                }
                else if ((VisObjectFamilyListInOrganism[i] == "Ribosome") & (Ribosomes.empty())) {
                    Text_Pos_Pol = "Pos_Pol_";
                }
                GenerateVisObjects(ofs, 5, VisObjectFamilyListInOrganism[i], Text_Pos_Pol + Processes[i] + ".shape[1]");
            }

            // Packaging DNA Replication message in the organism

            VisObjectFamily = DNAP;
            Process = Processes[0];


            ofs << in+ in+ in+ in+ in+ Process << "s[ObjID] = lccsimulation_pb2.MState_" << Process << "(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Objects_" << VisObjectFamily << " = VisObjects_" << VisObjectFamily << "," << endl;
            
            std::string Text_Pos_Pol = "self.State.Pos_Pol_";
            std::string Text_Dir_Pol = "self.State.Dir_Pol_";
            if (DNAPs.empty()) {
                Text_Pos_Pol = "Pos_Pol_";
                Text_Dir_Pol = "Dir_Pol_";
            }

            ofs << in+ in+ in+ in+ in+ in+ "Pos_" << VisObjectFamily << "_bp = " << Text_Pos_Pol << Process << "[i], " << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Dir_" << VisObjectFamily << " = " << Text_Dir_Pol << Process << "[i]," << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;

            // Packaging Transcription in the organism

            VisObjectFamily = RNAP;
            Process = Processes[1];


            ofs << in+ in+ in+ in+ in+ Process << "s[ObjID] = lccsimulation_pb2.MState_" << Process << "(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Objects_" << VisObjectFamily << " = VisObjects_" << VisObjectFamily << "," << endl;

            Text_Pos_Pol = "self.State.Pos_Pol_";
            Text_Dir_Pol = "self.State.Dir_Pol_";
            if (RNAPs.empty()) {
                Text_Pos_Pol = "Pos_Pol_";
                Text_Dir_Pol = "Dir_Pol_";
            }

            ofs << in+ in+ in+ in+ in+ in+ "Pos_" << VisObjectFamily << "_bp = " << Text_Pos_Pol << Process << "[i], " << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Dir_" << VisObjectFamily << " = " << Text_Dir_Pol << Process << "[i]," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Pos_" << VisObjectFamily << "_Gene = " << Text_Pos_Pol << "Template_" << Process << "[i], " << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;

            // Packaging Translation in the organism

            VisObjectFamily = Ribosome;
            Process = Processes[2];

            ofs << in+ in+ in+ in+ in+ Process << "s[ObjID] = lccsimulation_pb2.MState_" << Process << "(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Objects_" << VisObjectFamily << " = VisObjects_" << VisObjectFamily << "," << endl;

            Text_Pos_Pol = "self.State.Pos_Pol_";
            if (Ribosomes.empty()) {
                Text_Pos_Pol = "Pos_Pol_";
            }

            ofs << in+ in+ in+ in+ in+ in+ "Pos_" << VisObjectFamily << "_nt = " << Text_Pos_Pol << Process << "[i], " << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Pos_" << VisObjectFamily << "_mRNA = " << Text_Pos_Pol << "Template_" << Process << "[i], " << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;

        }
    }

    ofs << endl;

    // Packaging current state
    ofs << in+ in+ in+ in+ "CurState = lccsimulation_pb2.MSimUnitState(";
    // current state items
    ofs << "SimulationStep=self.SimM.GetSimStep(), ";
    ofs << "SimulatedTime=self.SimM.GetSimTime(), ";
    ofs << "Objects=VisObjects, ";
    ofs << "Counts=Counts, ";
    for (auto& process : Processes) {
        ofs << process + "=" + process + "s, ";
    } ofs << ")" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "response = lccsimulation_pb2.MRunData(State=CurState, Info='Any arbitrary message')" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "time.sleep(max(1./14 - (time.time() - start), 0))" << endl;
    ofs << in+ in+ in+ in+ "yield response" << endl;
    ofs << endl;
    ofs << in+ in+ "return response_messages()" << endl;
    ofs << endl;

    ofs << in+ "def Pause(self, request, context):" << endl;
    ofs << in+ in+ "print('pausing')" << endl;
    ofs << in+ in+ "self.IsRunning = False" << endl;
    ofs << in+ in+ "return lccsimulation_pb2.MControlSimulationResponse()" << endl;
    ofs << endl;

    ofs << in+ "def Stop(self, request, context):" << endl;
    ofs << in+ in+ "print('stopping')" << endl;
    ofs << in+ in+ "return lccsimulation_pb2.MControlSimulationResponse()" << endl;
    ofs << endl;

    ofs << "def main():   # add verbose" << endl;
    ofs << in+ "server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))" << endl;
    ofs << in+ "lccsimulation_pb2_grpc.add_LCCSimulationServicer_to_server(LCCSimulation(), server)" << endl;
    ofs << endl;
    ofs << in+ "server.add_insecure_port('[::]:50051')" << endl;
    ofs << in+ "server.start()" << endl;
    ofs << in+ "server.wait_for_termination()" << endl;
    ofs << endl;
    
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "main()" << endl;
    ofs << endl;


    std::cout << "  Simulation Server has been generated: ";

}
