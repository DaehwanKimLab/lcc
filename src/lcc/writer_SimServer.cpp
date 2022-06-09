#include "writer.h"

using namespace std;

void FWriter::SimServer() {
    std::cout << "Generating simulation server..." << std::endl;

    // write SimServer.py
    std::ofstream ofs(Option.SimServerFile.c_str());
    std::string endl = "\n";

    ofs << endl;
    ofs << "## lcc" << endl;
    ofs << "import random" << endl;
    ofs << "import SimModule" << endl;
    ofs << "import SimVis2D as SimV" << endl;
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

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto Organisms = Context.GetUniqueContainers_LocationList("Organism");

    auto Chromosomes = Context.GetSubList_MoleculeList("Chromosome");
    auto DNAPs = Context.GetSubList_MoleculeList("DNAP");

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
    ofs << in+ in+ "self.State = SimModule.FState()" << endl;
    ofs << in+ in+ "self.Data = SimModule.FDataset()" << endl;
    ofs << in+ in+ "self.DataManager = SimModule.FDataManager()" << endl;
    ofs << in+ in+ "self.SimM = SimModule.FSimulation(self.State, self.Data, self.DataManager)" << endl;
    ofs << endl;
    ofs << in+ in+ "# Initialize model" << endl;
    ofs << in+ in+ "self.SimM.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ "SimUnitTime = 0.2" << endl;

    if (!Context.ThresholdList.empty()) {
        if (Context.ThresholdList[0].first == "Am") {
            ofs << in+ in+ "SimM.Receptivity(20000)   # Pass time" << endl;
            ofs << endl;
        }
    }

    ofs << in+ in+ "# Setup Static Objects" << endl;
    
    ofs << in+ in+ "InitVisObjects = []" << endl;
    ofs << in+ in+ "ZeroVec = lccsimulation_pb2.MVector3(X=0, Y=0, Z=0)" << endl;
    ofs << in+ in+ "UnitVec = lccsimulation_pb2.MVector3(X=1, Y=1, Z=1)" << endl;
    ofs << in+ in+ "White = lccsimulation_pb2.MVector3(X=255, Y=255, Z=255)" << endl;
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
    ofs << in+ in+ in+ in+ "Color=Blue))" << endl;
    ofs << endl;

    // temporary static molecule 'distribution'
// TODO: take dynamic distribution

//for (int i = 0; i < MolLoc.size(); i++) {
//    // instantiate
//    ofs << in+ in+ MolLoc[i]->Name << " = SimV.FMolecule('" << MolLoc[i]->Name << "', '" << Identity[i] << "', ";

    if (!MolLoc.empty()) {
    
        ofs << in+ in+ "# Static Glucose" << endl;
        ofs << in+ in+ "for Pos in L.Particle_XY_Static:" << endl;
        ofs << in+ in+ in+ "NewObj = lccsimulation_pb2.MVisObjectData(" << endl;
        ofs << in+ in+ in+ in+ "ID=0, " << endl;
        ofs << in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_GLUCOSE," << endl;
        ofs << in+ in+ in+ in+ "Position=lccsimulation_pb2.MVector3(X=Pos[0], Y=Pos[1], Z=0)," << endl;
        ofs << in+ in+ in+ in+ "Rotation=ZeroVec," << endl;
        ofs << in+ in+ in+ in+ "Scale=lccsimulation_pb2.MVector3(X=10,Y=10,Z=10)," << endl;
        ofs << in+ in+ in+ in+ "Color=White)" << endl;
        ofs << in+ in+ in+ "" << endl;
        ofs << in+ in+ in+ "InitVisObjects.append(NewObj)" << endl;
        ofs << in+ in+ "" << endl;

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
            ofs << in+ in+ in+ in+ "Color = GenRandomColor()," << endl;
            ofs << in+ in+ in+ "))" << endl;
            ofs << endl;
        }
    }

    ofs << in+ in+ "# DNA Init Data" << endl;
    ofs << in+ in+ "DNA_Init = []" << endl;
    ofs << endl;

    if (!Chromosomes.empty()) {

        ofs << in+ in+ "# DNA Annotation" << endl;
        ofs << in+ in+ "DNA_Annotations = []" << endl;
        ofs << endl;

        ofs << in+ in+ "# Temporary halving" << endl;
        ofs << in+ in+ "Sequence=self.SimM.State.OpenFASTADatabase(r'./Database/EscherichiaColi.fasta')[0]" << endl;
        ofs << endl;

        ofs << in+ in+ "DNA_Annotations.append(lccsimulation_pb2.MDNA_AnnotationData(" << endl;
        //ofs << in+ in+ in+ "Sequence='ACGT'," << endl;
        //ofs << in+ in+ in+ "Sequence=self.SimM.State.OpenFASTADatabase(r'./Database/EscherichiaColi.fasta')," << endl;
        ofs << in+ in+ in+ "Sequence=Sequence," << endl;
        ofs << in+ in+ in+ "Gene_StartIndex_bp=self.State.Pos_Gene_Start_bp[0]," << endl;
        ofs << in+ in+ in+ "Gene_EndIndex_bp=self.State.Pos_Gene_End_bp[0]," << endl;
        ofs << in+ in+ in+ "Gene_Symbol=self.State.Name_Genes," << endl;
        ofs << in+ in+ "))" << endl;
        ofs << endl;
        
        ofs << in+ in+ "# DNA Position" << endl;
        ofs << in+ in+ "DNA_Positions = []" << endl;
        ofs << endl;
        ofs << in+ in+ "Positions = []" << endl;
        ofs << in+ in+ "for pos in self.State.Pos_Ref:" << endl;
        ofs << in+ in+ in+ "Positions.append(lccsimulation_pb2.MVector3(X=pos[0], Y=pos[1], Z=pos[2]))" << endl;
        ofs << endl;
        ofs << in+ in+ "DNA_Positions.append(lccsimulation_pb2.MDNA_PositionData(" << endl;
        ofs << in+ in+ in+ "Points=Positions" << endl;
        ofs << in+ in+ "))" << endl;
        ofs << endl;
        
        ofs << in+ in+ "DNA_Init = lccsimulation_pb2.MDNA_InitData(" << endl;
        ofs << in+ in+ in+ "DNA_Annotations=DNA_Annotations," << endl;
        ofs << in+ in+ in+ "DNA_Positions=DNA_Positions," << endl;
        ofs << in+ in+ ")" << endl;
    }

    ofs << in+ in+ "return lccsimulation_pb2.MInitData(InitObjects=InitVisObjects, InitDNA=DNA_Init)" << endl;
    ofs << endl;

    // TODO: stream run
    ofs << in+ "def Run(self, request, context):" << endl;
    ofs << in+ in+ "print('running')" << endl;
    ofs << in+ in+ "self.IsRunning = True" << endl;
    ofs << endl;
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
    ofs << in+ in+ in+ in+ "VisObjects = {} # map from id --> VisObjectData" << endl;
    ofs << in+ in+ in+ in+ "DNAReplications = {} # map from id --> VisObjectData" << endl;
    ofs << endl;

    if (!Organisms.empty()) {
        for (auto& organism : Organisms) {
            auto Organism = dynamic_cast<FOrganism*>(organism);
            ofs << in+ in+ in+ in+ "ReplicationCompletionRate = self.SimM.GetReplicationCompletionRate('" << Organism->Name << "')" << endl;

            ofs << in+ in+ in+ in+ "X, Y, Angle = self.SimM.GetPositionXYAngleByName('" << organism->Name << "')" << endl;
            ofs << in+ in+ in+ in+ "for i in range(len(X)):" << endl;
            ofs << in+ in+ in+ in+ in+ "PosVec = lccsimulation_pb2.MVector3(X=X[i], Y=Y[i], Z=0)" << endl;
            ofs << in+ in+ in+ in+ in+ "RotVec = lccsimulation_pb2.MVector3(X=0, Y=Angle[i] * (180/np.pi), Z=0)" << endl;
            ofs << endl;
            ofs << in+ in+ in+ in+ in+ "ObjID = i + 1   # ID 0 is used for static objects" << endl;
            ofs << in+ in+ in+ in+ in+ "VisObjects[ObjID] = lccsimulation_pb2.MVisObjectData(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_" << Utils::UpperCaseStr(Organism->Species) << "," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Position=PosVec," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "Rotation=RotVec," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "# Scale =, # Scale doesn't change, leave it out" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "# Color =, # Color doesn't change, leave it out" << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;
            ofs << endl;
            ofs << in+ in+ in+ in+ in+ "DNAReplications[ObjID] = lccsimulation_pb2.MState_DNAReplication(" << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
            ofs << in+ in+ in+ in+ in+ in+ "ReplicationCompletionRate = ReplicationCompletionRate[i][0] * 100," << endl;
            ofs << in+ in+ in+ in+ in+ ")" << endl;
        }
    }

    ofs << endl;
    ofs << in+ in+ in+ in+ "CurState = lccsimulation_pb2.MSimUnitState(SimulationStep=self.SimM.GetSimStep(), SimulatedTime=self.SimM.GetSimTime(), Objects=VisObjects, DNAReplication=DNAReplications)" << endl;
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
