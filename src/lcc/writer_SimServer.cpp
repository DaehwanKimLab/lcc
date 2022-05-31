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
    ofs << in+ "return lccsimulation_pb2.Vector3(X=RandomInt(256),Y=RandomInt(256),Z=RandomInt(256))" << endl;
    ofs << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto OrgNames = Context.GetUniqueNames_LocationList("Organism");

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
    for (int i = 0; i < OrgNames.size(); i++) {
        ofs << in+ in+ "self." << OrgNames[i] << " = None" << endl;
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
    ofs << in+ in+ "PetriDish = SimV.FEnvironment()" << endl;
    ofs << endl;

    if (!MolLoc.empty()) {
        // Distribution settings (hardcoding)
        //                                              L,              qL
        std::vector<std::string>    Identity            {"Glucose",     "Autoinducer",     "Mol3"     };
        std::vector<std::string>    Color               {"Blue",      "Black",     "Yellow"     };
        std::vector<int>            MaxBrightness       {200,           150,        50          };

        std::vector<std::string>    Pattern             {"heatmap",     "heatmap",  "particles" };
        std::vector<std::string>    NormalizationType   {"linear",      "linear",   "particles" };
        std::vector<int>            ReductionFactor     {5,             2,          1           };

        std::vector<std::string>    MaxStatic;
        std::vector<std::string>    ContourLine;

    //    if (Option.bDebug)
    //    {
    //                                Color =             {"Green",       "Blue",     "Green"     };
    //                                MaxStatic =         {"True",        "False",    "particles" };
    //                                ContourLine =       {"True",        "False",    "False" };
    //    }

        // Instantiate Molecules for Distribution
        // auto& MolLoc = Context.GetSubList_LocationList("Molecule");

        for (int i = 0; i < MolLoc.size(); i++) {
            // instantiate
            ofs << in+ in+ MolLoc[i]->Name << " = SimV.FMolecule('" << MolLoc[i]->Name << "', '" << Identity[i] << "', ";
            ofs << MolLoc[i]->Coord[0] << ", " << MolLoc[i]->Coord[1] << ")" << endl;

            // set color
            ofs << in+ in+ MolLoc[i]->Name << ".SetColor('" << Color[i] << "', " << MaxBrightness[i] << ")" << endl;

            float DrawingThreshold = 0;
            std::string bRaw = "False";
            if (MolLoc[i]->Coord[0] == -1) {
                DrawingThreshold = Context.GetInitialCountByName_CountList(MolLoc[i]->Name);
                bRaw = "True";
            }

            // set pattern
            ofs << in+ in+ MolLoc[i]->Name << ".SetPattern('" << Pattern[i] << "', '" << NormalizationType[i] << "', " << ReductionFactor[i] << ", " << Utils::SciFloat2Str(DrawingThreshold) << ", " << bRaw;
            if (!MaxStatic.empty())     { ofs << ", maxstatic=" << MaxStatic[i]; }
            if (!ContourLine.empty())   { ofs << ", contourline=" << ContourLine[i]; }
            ofs << ",)" << endl;

            ofs << endl;
        }
    }

    if (!OrgNames.empty()) {
            // Distribution settings (hardcoding)
        std::vector<std::string>    Radar_Switch;

    //    if (Option.bDebug)
    //    {
        Radar_Switch =      {"False",};
    //    }

        // Instantiate Organisms
        for (int i = 0; i < OrgNames.size(); i++) {
            ofs << in+ in+ "self." << OrgNames[i] << " = SimV.FOrganism('" << OrgNames[i] << "', " << "'Ecoli'" << ")" << endl; // TODO: Get Species later

            // set radar
            if (!Radar_Switch.empty())     { ofs << in+ in+ "self." << OrgNames[i] << ".SetRadar(switch=" << Radar_Switch[i] << ")" << endl; }
        }
        ofs << endl;

        // Initialize Organisms
        for (auto& OrganismName : OrgNames) {
            ofs << in+ in+ "self." << OrganismName << ".Initialize()" << endl;
    //        ofs << "# ";
            if (!Context.ThresholdList.empty()) {
                if (Context.ThresholdList[0].first == "Am") {
                    ofs << in+ in+ "self." << OrganismName << ".Receptivity(20000)   # Pass time" << endl;
                }
            }
        }
        ofs << endl;
    }

    ofs << in + in + "# Setup Static Objects" << endl;
    
    ofs << in + in + "InitVisObjects = []" << endl;
    ofs << in + in + "ZeroVec = lccsimulation_pb2.Vector3(X=0, Y=0, Z=0)" << endl;
    ofs << in + in + "UnitVec = lccsimulation_pb2.Vector3(X=1, Y=1, Z=1)" << endl;
    ofs << in + in + "White = lccsimulation_pb2.Vector3(X=255, Y=255, Z=255)" << endl;
    ofs << in + in + "Blue = lccsimulation_pb2.Vector3(X=143, Y=186, Z=255)" << endl;
    ofs << in + in + "Yellow = lccsimulation_pb2.Vector3(X=255, Y=255, Z=102)" << endl;
    ofs << endl;

    ofs << in + in + "# Static Petri Dish" << endl;
    ofs << in + in + "InitVisObjects.append(lccsimulation_pb2.VisObjectData(" << endl;
    ofs << in + in + in + in + "ID=0, " << endl;
    ofs << in + in + in + in + "ObjType=lccsimulation_pb2.VisObjectType.M_PETRI_DISH," << endl;
    ofs << in + in + in + in + "Position=lccsimulation_pb2.Vector3(X=100, Y=100, Z=0)," << endl;
    ofs << in + in + in + in + "Rotation=ZeroVec," << endl;
    ofs << in + in + in + in + "Scale=lccsimulation_pb2.Vector3(X=self.SimM.State.Dimension_X, Y=self.SimM.State.Dimension_Y, Z=1)," << endl;
    ofs << in + in + in + in + "Color=Blue))" << endl;
    ofs << endl;

    // temporary static molecule 'distribution'
// TODO: take dynamic distribution

//for (int i = 0; i < MolLoc.size(); i++) {
//    // instantiate
//    ofs << in + in + MolLoc[i]->Name << " = SimV.FMolecule('" << MolLoc[i]->Name << "', '" << Identity[i] << "', ";

    if (!MolLoc.empty()) {
    
        ofs << in + in + "# Static Glucose" << endl;
        ofs << in + in + "for Pos in L.Particle_XY_Static:" << endl;
        ofs << in + in + in + "NewObj = lccsimulation_pb2.VisObjectData(" << endl;
        ofs << in + in + in + in + "ID=0, " << endl;
        ofs << in + in + in + in + "ObjType=lccsimulation_pb2.VisObjectType.M_GLUCOSE," << endl;
        ofs << in + in + in + in + "Position=lccsimulation_pb2.Vector3(X=Pos[0], Y=Pos[1], Z=0)," << endl;
        ofs << in + in + in + in + "Rotation=ZeroVec," << endl;
        ofs << in + in + in + in + "Scale=lccsimulation_pb2.Vector3(X=10,Y=10,Z=10)," << endl;
        ofs << in + in + in + in + "Color=White)" << endl;
        ofs << in + in + in + "" << endl;
        ofs << in + in + in + "InitVisObjects.append(NewObj)" << endl;
        ofs << in + in + "" << endl;

    }

    ofs << in + in + "# Setup Dynamic Objects" << endl;

    if (!OrgNames.empty()) {
        
        ofs << in+ in+ "# Dynamic Ecoli" << endl;

        for (auto& orgName : OrgNames) {
            ofs << in+ in+ "for i in range(len(self." << orgName << ".X)):" << endl;
            ofs << in+ in+ in+ "PosVec = lccsimulation_pb2.Vector3(X=self." << orgName << ".X[i], Y=self." << orgName << ".Y[i], Z=0)" << endl;
            ofs << in+ in+ in+ "RotVec = lccsimulation_pb2.Vector3(X=0, Y=self." << orgName << ".Angle[i] * (180/np.pi), Z=0)" << endl;
            ofs << in+ in+ in+ "" << endl;
            ofs << in+ in+ in+ "InitVisObjects.append(lccsimulation_pb2.VisObjectData(" << endl;
            ofs << in+ in+ in+ in+ "ID=i + 1," << endl;
            ofs << in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_ECOLI," << endl;
            ofs << in+ in+ in+ in+ "Position=PosVec," << endl;
            ofs << in+ in+ in+ in+ "Rotation=RotVec," << endl;
            ofs << in+ in+ in+ in+ "Scale = lccsimulation_pb2.Vector3(X=0.2,Y=0.2,Z=0.2)," << endl;
            ofs << in+ in+ in+ in+ "Color = GenRandomColor()," << endl;
            ofs << in+ in+ in+ "))" << endl;
            ofs << in+ in+ "" << endl;
        }

    }

    ofs << in+ in+ "# DNA Init Data" << endl;
    ofs << in+ in+ "DNA_Init = []" << endl;
    ofs << endl;

    if (!Chromosomes.empty()) {

        ofs << in+ in+ "# DNA Annotation" << endl;
        ofs << in+ in+ "DNA_Annotations = []" << endl;
        ofs << endl;
        ofs << in+ in+ "DNA_Annotations.append(lccsimulation_pb2.DNA_AnnotationData(" << endl;
        ofs << in+ in+ in+ "Sequence='ACGT'," << endl;
        ofs << in+ in+ in+ "Gene_StartIndex_nm=self.State.Pos_Gene_Start_nm," << endl;
        ofs << in+ in+ in+ "Gene_EndIndex_nm=self.State.Pos_Gene_End_nm," << endl;
        ofs << in+ in+ in+ "Gene_Symbol=self.State.Name_Genes," << endl;
        ofs << in+ in+ "))" << endl;
        ofs << endl;
        
        ofs << in+ in+ "# DNA Position" << endl;
        ofs << in+ in+ "DNA_Positions = []" << endl;
        ofs << endl;
        ofs << in+ in+ "Positions = []" << endl;
        ofs << in+ in+ "for pos in self.State.Pos_Ref:" << endl;
        ofs << in+ in+ in+ "Positions.append(lccsimulation_pb2.Vector3(X=pos[0], Y=pos[1], Z=pos[2]))" << endl;
        ofs << endl;
        ofs << in+ in+ "DNA_Positions.append(lccsimulation_pb2.DNA_PositionData(" << endl;
        ofs << in+ in+ in+ "Points=Positions" << endl;
        ofs << in+ in+ "))" << endl;
        ofs << endl;
        
        ofs << in+ in+ "DNA_Init = lccsimulation_pb2.DNA_InitData(" << endl;
        ofs << in+ in+ in+ "DNA_Annotations=DNA_Annotations," << endl;
        ofs << in+ in+ in+ "DNA_Positions=DNA_Positions," << endl;
        ofs << in+ in+ ")" << endl;
    }

    ofs << in + in + "return lccsimulation_pb2.InitData(InitObjects=InitVisObjects, InitDNA=DNA_Init)" << endl;
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

    for (auto& OrganismName : OrgNames) {
        ofs << in+ in+ in+ in+ in+ "self." << OrganismName << ".Receptivity(20)" << endl;
        ofs << in+ in+ in+ in+ in+ "self." << OrganismName << ".SetPosition(self.SimM.GetPositionXYAngleByName('" << OrganismName << "'))" << endl;
        if (!DNAPs.empty()) {
            ofs << in+ in+ in+ in+ in+ "self." << OrganismName << ".SetReplicationCompletionRate(self.SimM.GetReplicationCompletionRate('" << OrganismName << "'))" << endl;
        }
        ofs << in+ in+ in+ in+ in+ "self." << OrganismName << ".IncrementSimCount()" << endl;
    }
    ofs << endl;
    ofs << in+ in+ in+ in+ in+ "ElapsedTime -= SimUnitTime" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Record data SimUnitState" << endl;
    ofs << in+ in+ in+ in+ "VisObjects = {} # map from id --> VisObjectData" << endl;
    ofs << endl;


    ofs << in+ in+ in+ in+ "for i in range(len(self.E.X)):" << endl;
    ofs << in+ in+ in+ in+ in+ "PosVec = lccsimulation_pb2.Vector3(X=self.E.X[i], Y=self.E.Y[i], Z=0)" << endl;
    ofs << in+ in+ in+ in+ in+ "RotVec = lccsimulation_pb2.Vector3(X=0, Y=self.E.Angle[i] * (180/np.pi), Z=0)" << endl;
    ofs << in+ in+ in+ in+ in+ "" << endl;
    ofs << in+ in+ in+ in+ in+ "ObjID = i + 1; # ID 0 is used for static objects" << endl;
    ofs << in+ in+ in+ in+ in+ "VisObjects[ObjID] = lccsimulation_pb2.VisObjectData(" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "ID=ObjID," << endl;
    ofs << in+ in+ in+ in+ in+ in+ "ObjType=lccsimulation_pb2.VisObjectType.M_ECOLI," << endl;
    ofs << in+ in+ in+ in+ in+ in+ "Position=PosVec," << endl;
    ofs << in+ in+ in+ in+ in+ in+ "Rotation=RotVec," << endl;
    ofs << in+ in+ in+ in+ in+ in+ "# Scale =, # Scale doesn't change, leave it out" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "# Color =, # Color doesn't change, leave it out" << endl;
    ofs << in+ in+ in+ in+ in+ ")" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "CurState = lccsimulation_pb2.SimUnitState(Time=self.E.SimCount, Objects=VisObjects)" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "response = lccsimulation_pb2.RunData(State=CurState, Info='Any arbitrary message')" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "time.sleep(max(1./14 - (time.time() - start), 0))" << endl;
    ofs << in+ in+ in+ in+ "yield response" << endl;
    ofs << endl;
    ofs << in+ in+ "return response_messages()" << endl;
    ofs << endl;

    ofs << in+ "def Pause(self, request, context):" << endl;
    ofs << in+ in+ "print('pausing')" << endl;
    ofs << in+ in+ "self.IsRunning = False" << endl;
    ofs << in+ in+ "return lccsimulation_pb2.ControlSimulationResponse()" << endl;
    ofs << endl;

    ofs << in+ "def Stop(self, request, context):" << endl;
    ofs << in+ in+ "print('stopping')" << endl;
    ofs << in+ in+ "return lccsimulation_pb2.ControlSimulationResponse()" << endl;
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
