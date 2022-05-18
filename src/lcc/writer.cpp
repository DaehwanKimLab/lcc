#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SetUpDefaultVariables(int InN_MoleculesAllowed, std::string InName_Pseudo, float InFloat_Default, int InInt_Default) {
    in = "    ";
    N_MoleculesAllowed = InN_MoleculesAllowed;
    Name_Pseudo = InName_Pseudo;
    Float_Init = InFloat_Default; // random initialized float
    Int_Init = InInt_Default; // random initialized int
}

std::string GetRegType(std::string Type)
{
    std::string RegType;

    if      (Type.find("Inhibition_Allosteric") != string::npos)   { RegType = "Regulatory_Inhibition_Allosteric"; }
    else if (Type.find("Inhibition_Competitive") != string::npos)  { RegType = "Regulatory_Inhibition_Competitive"; }
    else if (Type.find("Activation_Allosteric") != string::npos)   { RegType = "Regulatory_Activation_Allosteric"; }
    else if (Type.find("Unregulated") != string::npos)             { RegType = "Unregulated"; }

    return RegType;
}

void FWriter::Initialize_StandardReaction(ofstream& ofs, std::string Type, std::string NameSpace_Pathway)
{
    std::string TypeText = Type;
    if (NameSpace_Pathway != "") {
        TypeText += "_" + NameSpace_Pathway;
    }

    ofs << in+ in+ "# " << TypeText << endl;

    // Standard vs. MichaelisMenten
    ofs << in+ in+ "self.Const_k_Reactant_" << TypeText << " = None" << endl;
    ofs << in+ in+ "self.Const_k_Product_" << TypeText << " = None" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << TypeText << " = None" << endl;
    }
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Product_" << i << "_" << TypeText << " = None" << endl;
    }
    ofs << endl;

    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Const_n_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << TypeText << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << TypeText << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << TypeText << " = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_StandardReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList, std::string NameSpace_Pathway)
{
    std::string TypeText = Type;
    if (NameSpace_Pathway != "") {
        TypeText += "_" + NameSpace_Pathway;
    }

    if (NameSpace_Pathway != "") {
        ReactionSubList = Context.FilterByPathway_ReactionList(ReactionSubList, NameSpace_Pathway);
    }

    // Types allowed: "Standard_Unregulated", "Standard_Inhibition", "Standard_Activation"
    // TODO: Multiplexing for regulation
    std::string RegType = GetRegType(Type);

    // for standard reactions
    std::vector<float> k;
    std::vector<float> krev;

    // common variables
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList(Name_Pseudo);
    std::vector<std::vector<int>> Idx_Reactants(N_MoleculesAllowed);
    std::vector<std::vector<int>> Idx_Products(N_MoleculesAllowed);

    // for regulatory mechanisms
    std::vector<float> K;
    std::vector<float> n;
    std::vector<int> Idx_Regulator; // reactant of the regulatory reaction

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetUniqueSubstrateIdx_ReactionList(ReactionSubList);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(ReactionSubList);

    // relevant reaction lists
    std::vector<FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    if (NameSpace_Pathway != "") {
        RegulatoryReactionSubList = Context.FilterByPathway_ReactionList(RegulatoryReactionSubList, NameSpace_Pathway);
    }


    for (auto& Reaction : ReactionSubList) {
        const auto& reaction = dynamic_cast<FStandardReaction *>(Reaction);

        k.push_back(reaction->k);
        krev.push_back(reaction->krev);

        std::vector<int> Idx_Reactant;
        std::vector<int> Idx_Product;

        for (auto& stoich : reaction->Stoichiometry) {
            int Idx = Context.GetIdxByName_MoleculeList(stoich.first);
            if (stoich.second <= 0) {
                Idx_Reactant.push_back(Idx);
            } else if (stoich.second >= 0) {
                Idx_Product.push_back(Idx);
            }
        }
        // Determines the number of substrates to handle (pseudomolecule is implemented to fill in blank spots for now)
        while (Idx_Reactant.size() != N_MoleculesAllowed) {
            Idx_Reactant.push_back(Idx_Pseudo);
        }
        while (Idx_Product.size() != N_MoleculesAllowed) {
            Idx_Product.push_back(Idx_Pseudo);
        }

        for (int i = 0; i < N_MoleculesAllowed; i++) {
            // std::cout << i << " : " << Idx_Reactant[i] << ", " << Idx_Product[i] << endl;
            Idx_Reactants[i].push_back(Idx_Reactant[i]);
            Idx_Products[i].push_back(Idx_Product[i]);
        }

        // Check Idx lengths
        Utils::Assertion((Idx_Reactants.size() == Idx_Products.size()), "# of parsed reactants and products do not match" );

        if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
            for (auto& Reg: RegulatoryReactionSubList) {
                const auto& reg = dynamic_cast<FRegulatoryReaction *>(Reg);

                bool Import = false;
                std::string reactant_reg;
                std::string product_reg;
                for (auto& stoich: reg->Stoichiometry) {
                    if (stoich.second < 0) { reactant_reg = stoich.first; }
                    else if (stoich.second > 0) { product_reg = stoich.first; }
                    // import only if the product of the regulatory reaction targets the current standard reaction
                    if (stoich.first == reaction->Name) {
                        Import = true;
                    }
                }
                if (Import) {
                    Idx_Regulator.push_back(Context.GetIdxByName_MoleculeList(reactant_reg));

                    K.push_back(reg->K);
                    n.push_back(reg->n);
                    break;
                }
            }
        }
    }

    ofs << in+ in+ "# " << TypeText << endl;

    ofs << in+ in+ "self.Const_k_Reactant_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(k) << "])" << endl;
    ofs << in+ in+ "self.Const_k_Product_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(krev) << "])" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
    }
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Product_" << i << "_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ "self.Const_n_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << TypeText << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << TypeText << " = np.asmatrix([" << Utils::Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void FWriter::Initialize_EnzymeReaction(ofstream& ofs, std::string Type, std::string NameSpace_Pathway)
{
    std::string TypeText = Type;
    if (NameSpace_Pathway != "") {
        TypeText += "_" + NameSpace_Pathway;
    }

    ofs << in+ in+ "# " << TypeText << endl;

    ofs << in+ in+ "self.Idx_Enz_" << TypeText << " = None" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Enz_Standard") != string::npos) {
        ofs << in+ in+ "self.Const_k_Reactant_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Const_k_Product_" << TypeText << " = None" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << TypeText << " = None" << endl;
        }
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Product_" << i << "_" << TypeText << " = None" << endl;
        }

    } else if (Type.find("Enz_MichaelisMenten") != string::npos) {
        ofs << in+ in+ "self.Const_kcat_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Const_KM_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Idx_EnzSub_" << TypeText << " = None" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Const_n_" << TypeText << " = None" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << TypeText << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << TypeText << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << TypeText << " = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_EnzymeReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList, std::string NameSpace_Pathway) // to be changed with reaction list
{
    std::string TypeText = Type;
    if (NameSpace_Pathway != "") {
        TypeText += "_" + NameSpace_Pathway;
    }

    if (NameSpace_Pathway != "") {
        ReactionSubList = Context.FilterByPathway_ReactionList(ReactionSubList, NameSpace_Pathway);
    }

    // TODO: Multiplexing for regulation
    std::string RegType = GetRegType(Type);

    // for standard reactions
    std::vector<float> k;
    std::vector<float> krev;

    // for michaelis menten kinetics
    std::vector<float> kcats;
    std::vector<float> KMs;
    std::vector<int> Idx_EnzSub;

    // common variables
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList(Name_Pseudo);
    std::vector<int> Idx_Enz; // for Enzyme where En is not included in the reaction
    std::vector<std::vector<int>> Idx_Reactants(N_MoleculesAllowed);
    std::vector<std::vector<int>> Idx_Products(N_MoleculesAllowed);

    // for regulatory mechanisms
    std::vector<float> K;
    std::vector<float> n;
    std::vector<int> Idx_Regulator;

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetUniqueSubstrateIdx_ReactionList(ReactionSubList);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(ReactionSubList);

    // relevant reaction lists
    std::vector<FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    if (NameSpace_Pathway != "") {
        RegulatoryReactionSubList = Context.FilterByPathway_ReactionList(RegulatoryReactionSubList, NameSpace_Pathway);
    }

    for (auto& Reaction : ReactionSubList) {
        const auto& EnzReaction = dynamic_cast<FEnzymaticReaction *>(Reaction);
        const auto& enzyme = dynamic_cast<FEnzyme *>(Context.GetMolecule_MoleculeList(EnzReaction->Enzyme));

        Idx_Enz.push_back(Context.GetIdxByName_MoleculeList(enzyme->Name));

        // Enz_Standard
        if ((EnzReaction->Type >= 10) & (EnzReaction->Type < 20)) {
            const auto& reaction = dynamic_cast<FEnz_StandardReaction *>(Reaction);

            k.push_back(reaction->k);
            krev.push_back(reaction->krev);

            std::vector<int> Idx_Reactant;
            std::vector<int> Idx_Product;

            for (auto& stoich : reaction->Stoichiometry) {
                int Idx = Context.GetIdxByName_MoleculeList(stoich.first);
                if (stoich.second <= 0) {
                    Idx_Reactant.push_back(Idx);
                } else if (stoich.second >= 0) {
                    Idx_Product.push_back(Idx);
                }
            }
            // Determines the number of substrates to handle (pseudomolecule is implemented to fill in blank spots for now)
            while (Idx_Reactant.size() != N_MoleculesAllowed) {
                Idx_Reactant.push_back(Idx_Pseudo);
            }
            while (Idx_Product.size() != N_MoleculesAllowed) {
                Idx_Product.push_back(Idx_Pseudo);
            }

            for (int i = 0; i < N_MoleculesAllowed; i++) {
                // std::cout << i << " : " << Idx_Reactant[i] << ", " << Idx_Product[i] << endl;
                Idx_Reactants[i].push_back(Idx_Reactant[i]);
                Idx_Products[i].push_back(Idx_Product[i]);
            }

            // Enz_Michaelis Menten
        } else if ((EnzReaction->Type >= 20) & (EnzReaction->Type < 30)) {
            bool bSuccess = false;
            for (auto& kinetics : enzyme->Kinetics) {
                // get substrate, kcat, KM when enzyme's substrate is found as a reactant of the reaction

                if (Reaction->CheckIfReactant(kinetics.first)) {
                    Idx_EnzSub.push_back(Context.GetIdxByName_MoleculeList(kinetics.first));
                    auto& constants = kinetics.second;

                    kcats.push_back(constants[0]);
                    KMs.push_back(constants[1]);
                    bSuccess = true;
                    break;
                }
            }
            Utils::Assertion (bSuccess, "Unable to load kinetic constants for Michaelis Menten Reaction");
        }

        if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
            for (auto& Reg: RegulatoryReactionSubList) {
                const auto& reg = dynamic_cast<FRegulatoryReaction *>(Reg);

                bool Import = false;
                std::string reactant_reg;
                std::string product_reg;
                for (auto& stoich: reg->Stoichiometry) {
                    if (stoich.second < 0) { reactant_reg = stoich.first; }
                    else if (stoich.second > 0) { product_reg = stoich.first; }
                    // import only if the product of the regulatory reaction targets the current standard reaction
                    if (stoich.first == EnzReaction->Name) {
                        Import = true;
                    }
                }
                if (Import) {
                    Idx_Regulator.push_back(Context.GetIdxByName_MoleculeList(reactant_reg));
                    K.push_back(reg->K);
                    n.push_back(reg->n);
                    break;
                }
            }
        }
    }

    // Check Idx lengths
    Utils::Assertion((Idx_Reactants.size() == Idx_Products.size()), "# of parsed reactants and products do not match");

    ofs << in+ in+ "# " << TypeText << endl;
    ofs << in+ in+ "self.Idx_Enz_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Enz) << "])" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Enz_Standard") != string::npos) {
        ofs << in+ in+ "self.Const_k_Reactant_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(k) << "])" << endl;
        ofs << in+ in+ "self.Const_k_Product_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(krev) << "])" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
        }
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Product_" << i << "_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
        }

    } else if (Type.find("Enz_MichaelisMenten") != string::npos) {
        ofs << in+ in+ "self.Const_kcat_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(kcats) << "])" << endl;
        ofs << in+ in+ "self.Const_KM_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(KMs) << "])" << endl;
        ofs << in+ in+ "self.Idx_EnzSub_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_EnzSub) << "])" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (TypeText.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ "self.Const_n_" << TypeText << " = np.array([" << Utils::JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << TypeText << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << TypeText << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << TypeText << " = np.asmatrix([" << Utils::Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void FWriter::Initialize_PolymeraseReaction_Matrix(ofstream& ofs, std::vector<std::vector<FMolecule *>> PolymeraseTypes)
{
    std::string BB_Ch, BB_RNA, BB_Protein;
    std::string Max_Ch, Max_RNA, Max_Protein;
    std::string Len_Ch, Len_RNA, Len_Protein;

    Max_Ch =        "self.MaxLen_NascentChromosomes";
    Max_RNA =       "self.MaxLen_NascentRNAs";
    Max_Protein =   "self.MaxLen_NascentProteins";
    BB_Ch =         "self.Freq_BB_Chromosomes";
    BB_RNA =        "self.Freq_BB_RNAs";
    BB_Protein =    "self.Freq_BB_Proteins";

    Len_Ch =        "self.Len_NascentChromosomes";
    Len_RNA =       "self.Len_NascentRNAs";
    Len_Protein =   "self.Len_NascentProteins";

    // For DNA Replication
    if (!PolymeraseTypes[0].empty() & !PolymeraseTypes[1].empty() & !PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_Ch        << " = None" << endl;
        ofs << in+ in+ BB_RNA       << " = None" << endl;
        ofs << in+ in+ BB_Protein   << " = None" << endl;
        ofs << in+ in+ Max_Ch       << " = None" << endl;
        ofs << in+ in+ Max_RNA      << " = None" << endl;
        ofs << in+ in+ Max_Protein  << " = None" << endl;
        ofs << in+ in+ Len_Ch       << " = None" << endl;
        ofs << in+ in+ Len_RNA      << " = None" << endl;
        ofs << in+ in+ Len_Protein  << " = None" << endl;

    // For Transcription (and Replication)
    } else if (!PolymeraseTypes[1].empty() & PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_Ch        << " = None" << endl;
        ofs << in+ in+ BB_RNA       << " = None" << endl;
        ofs << in+ in+ Max_Ch       << " = None" << endl;
        ofs << in+ in+ Max_RNA      << " = None" << endl;
        ofs << in+ in+ Len_Ch       << " = None" << endl;
        ofs << in+ in+ Len_RNA      << " = None" << endl;

    // For Translation
    } else if (PolymeraseTypes[0].empty() & PolymeraseTypes[1].empty() & !PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_RNA   << " = None" << endl;
        ofs << in+ in+ BB_Protein   << " = None" << endl;
        ofs << in+ in+ Max_RNA      << " = None" << endl;
        ofs << in+ in+ Max_Protein  << " = None" << endl;
        ofs << in+ in+ Len_RNA      << " = None" << endl;
        ofs << in+ in+ Len_Protein  << " = None" << endl;

    // For Replication
    } else if (!PolymeraseTypes[0].empty() & PolymeraseTypes[1].empty() & PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_Ch        << " = None" << endl;
        ofs << in+ in+ Max_Ch       << " = None" << endl;
        ofs << in+ in+ Len_Ch       << " = None" << endl;
    }
}

void FWriter::Initialize_PolymeraseReaction_DNAP(ofstream& ofs, std::vector<FMolecule *> Polymerases)
{
    // Polymerases for future uses

    // std::string InitiationSite
    // std::string TerminationSite

    // initialize common for all polymerase
    Initialize_PolymeraseReaction_Index(ofs, dynamic_cast<FPolymerase *>(Polymerases[0])->Process);
}

void FWriter::Initialize_PolymeraseReaction_RNAP(ofstream& ofs, std::vector<FMolecule *> Polymerases)
{
    // Polymerases for future uses

    // std::string InitiationSite
    // std::string TerminationSite

    // initialize common for all polymerase
    Initialize_PolymeraseReaction_Index(ofs, dynamic_cast<FPolymerase *>(Polymerases[0])->Process);
}

void FWriter::Initialize_PolymeraseReaction_Ribosome(ofstream& ofs, std::vector<FMolecule *> Polymerases)
{
    // Polymerases for future uses

    // std::string InitiationSite
    // std::string TerminationSite

    // initialize common for all polymerase
    Initialize_PolymeraseReaction_Index(ofs, dynamic_cast<FPolymerase *>(Polymerases[0])->Process);
    ofs << in+ in+ "self.Idx_mRNAInRNA = None" << endl;

}

void FWriter::Initialize_PolymeraseReaction_Index(ofstream& ofs, std::string Process)
{
    // Initiation
    ofs << in+ in+ "# " << Process << ": Initialize Initiation Reaction" << endl;
    ofs << in+ in+ "self.Idx_Pol_" << Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset_" << Process << " = None" << endl; // local indexing within the template population for mRNA in RNA for protein translation
    ofs << in+ in+ "self.Weight_" << Process << " = None" << endl;
    ofs << in+ in+ "self.Pol_Threshold_" << Process << " = None" << endl; // here
    ofs << endl;

    // Elongation
    ofs << in+ in+ "# " << Process << ": Initialize Elongation Reaction" << endl;
    ofs << in+ in+ "self.Rate_" << Process << " = None" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Idx_PolSub_" << Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Process << " = None" << endl;
    ofs << endl;
}
void FWriter::SetUp_PolymeraseReaction_Matrix(ofstream& ofs, std::vector<std::vector<FMolecule *>> PolymeraseTypes)
{
    std::string BB_Ch, BB_RNA, BB_Protein;
    std::string Max_Ch, Max_RNA, Max_Protein;
    std::string Len_Ch, Len_RNA, Len_Protein;

    Max_Ch =        "self.MaxLen_NascentChromosomes = np.asmatrix(np.load(r'./Database/Len_ChromosomesInGenome.npy'))";
    Max_RNA =       "self.MaxLen_NascentRNAs = np.asmatrix(np.load(r'./Database/Len_RNAs.npy'))";
    Max_Protein =   "self.MaxLen_NascentProteins = np.asmatrix(np.load(r'./Database/Len_Proteins.npy'))";
    BB_Ch =         "self.Freq_BB_Chromosomes = np.asmatrix(np.load(r'./Database/Freq_NTsInChromosomesInGenome.npy'))";
    BB_RNA =        "self.Freq_BB_RNAs = np.asmatrix(np.load(r'./Database/Freq_NTsInRNAs.npy'))";
    BB_Protein =    "self.Freq_BB_Proteins = np.asmatrix(np.load(r'./Database/Freq_AAsInProteins.npy'))";

    Len_Ch =        "self.Len_NascentChromosomes = np.asmatrix(np.full([10, self.Freq_BB_Chromosomes.shape[0]], -1))";
    Len_RNA =       "self.Len_NascentRNAs = np.asmatrix(np.full([10, self.Freq_BB_RNAs.shape[0]], -1))";
    Len_Protein =   "self.Len_NascentProteins = np.asmatrix(np.full([10, self.Freq_BB_Proteins.shape[0]], -1))";

    // For DNA Replication
    if (!PolymeraseTypes[0].empty() & !PolymeraseTypes[1].empty() & !PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_Ch        << endl;
        ofs << in+ in+ BB_RNA       << endl;
        ofs << in+ in+ BB_Protein   << endl;
        ofs << in+ in+ Max_Ch       << endl;
        ofs << in+ in+ Max_RNA      << endl;
        ofs << in+ in+ Max_Protein  << endl;
        ofs << in+ in+ Len_Ch       << endl;
        ofs << in+ in+ Len_RNA      << endl;
        ofs << in+ in+ Len_Protein  << endl;

        // For Transcription (and Replication)
    } else if (!PolymeraseTypes[1].empty() & PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_Ch        << endl;
        ofs << in+ in+ BB_RNA       << endl;
        ofs << in+ in+ Max_Ch       << endl;
        ofs << in+ in+ Max_RNA      << endl;
        ofs << in+ in+ Len_Ch       << endl;
        ofs << in+ in+ Len_RNA      << endl;

        // For Translation
    } else if (PolymeraseTypes[0].empty() & PolymeraseTypes[1].empty() & !PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_RNA       << endl;
        ofs << in+ in+ BB_Protein   << endl;
        ofs << in+ in+ Max_RNA      << endl;
        ofs << in+ in+ Max_Protein  << endl;
        ofs << in+ in+ Len_RNA      << endl;
        ofs << in+ in+ Len_Protein  << endl;

        // For Replication
    } else if (!PolymeraseTypes[0].empty() & PolymeraseTypes[1].empty() & PolymeraseTypes[2].empty()) {
        ofs << in+ in+ BB_Ch        << endl;
        ofs << in+ in+ Max_Ch       << endl;
        ofs << in+ in+ Len_Ch       << endl;
    }
}

void FWriter::SetUp_PolymeraseReaction_DNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases)
{
    // std::string InitiationSite
    // std::string TerminationSite

    int Threshold = 50;

    SetUp_PolymeraseReaction_Index(ofs, Polymerases, Threshold);

}

void FWriter::SetUp_PolymeraseReaction_RNAP(ofstream& ofs, std::vector<FMolecule*> Polymerases)
{
    // std::string InitiationSite
    // std::string TerminationSite

    int Threshold = 1;

    SetUp_PolymeraseReaction_Index(ofs, Polymerases, Threshold);

}
void FWriter::SetUp_Idx_mRNAInRNA(ofstream& ofs)
{
    std::vector<int> Idx_mRNAInRNA = Context.GetLocalIdxList_MoleculeList(Context.GetSubList_MoleculeList("mRNA"), Context.GetSubList_MoleculeList("RNA"));
    Utils::Assertion(Idx_mRNAInRNA.size() == Context.GetSubList_MoleculeList("Protein").size(), "ERROR: # of mRNA and Proteins do not match. mRNAs: " + to_string(Idx_mRNAInRNA.size()) + " | Proteins: " + to_string(Context.GetSubList_MoleculeList("Protein").size()));
    ofs << in+ in+ "Idx_mRNAInRNA = np.array([" << Utils::JoinInt2Str_Idx(Idx_mRNAInRNA) << "], ndmin=2)" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset_" << "Translation" << " = Idx_mRNAInRNA" << endl; // local indexing within the template population for mRNA in RNA for protein translation

}

void FWriter::SetUp_PolymeraseReaction_Ribosome(ofstream& ofs, std::vector<FMolecule*> Polymerases)
{
    // std::string InitiationSite
    // std::string TerminationSite

    int Threshold = 1;

    SetUp_PolymeraseReaction_Index(ofs, Polymerases, Threshold);
    SetUp_Idx_mRNAInRNA(ofs);

}

void FWriter::SetUp_PolymeraseReaction_Index(ofstream& ofs, std::vector<FMolecule *> Polymerases, int Threshold)
{
    // Polymerase Example
    FPolymerase* Pol = dynamic_cast<FPolymerase *>(Polymerases[0]);
//    std::cout << "Polymerase: " << Pol->Name << " | Process: " << Pol->Process << " | TemplateClass: " << Pol->TemplateClass << " | TargetClass: " << Pol->TargetClass << std::endl;

    std::string MaxLen =            "self.MaxLen_Nascent"   + Pol->TargetClass + "s";
    std::string BuildingBlockFreq = "self.Freq_BB_"         + Pol->TargetClass + "s";
    std::string Len =               "self.Len_Nascent"      + Pol->TargetClass + "s";

    std::vector<float>          Rate;
    std::vector<int>            Idx_Pol,    Idx_Template,       Idx_Target,     Idx_Local,  Idx_PolSub,     Idx_PolBB;
    // For debugging purposes
    std::vector<std::string>    Names_Pol,  Names_Template,     Names_Target,               Names_PolSub,   Names_PolBB;

    // Set Rate array
    for (auto& polymerase : Polymerases) {
        auto Polymerase = dynamic_cast<FPolymerase *>(polymerase);
        Rate.push_back(Polymerase->Rate);
    }

    // Set Pol info
    Names_Pol = Context.GetNameList_MoleculeList(Polymerases);
    Idx_Pol = Context.GetIdxList_MoleculeList(Polymerases);

    // Set Template and targets
    Names_Template = Context.GetNameListByType_MoleculeList(Pol->TemplateClass);
    Idx_Template = Context.GetIdxListByType_MoleculeList(Pol->TemplateClass);

    std::vector<FMolecule *> TargetList = Context.GetSubList_MoleculeList(Pol->TargetClass);
    Names_Target = Context.GetNameListByType_MoleculeList(Pol->TargetClass);
    Idx_Target = Context.GetIdxListByType_MoleculeList(Pol->TargetClass);

    // PolymeraseReaction Example
    auto PolymeraseReaction = Context.GetReactionByPolymeraseName_ReactionList(Polymerases[0]->Name);

    // Substrates
    Idx_PolSub = Context.AddUniqueSubstrateIdxToIdxList(PolymeraseReaction, Idx_PolSub);
    for (auto idx : Idx_PolSub) {
        Names_PolSub.push_back(Context.MoleculeList[idx]->GetName());
    }

    // BuildingBlocks
    Names_PolBB = dynamic_cast<FPolymeraseReaction *>(PolymeraseReaction)->BuildingBlocks;
    Idx_PolBB = Context.GetIdxByStrList_MoleculeList(Names_PolBB);

    // Initiation
    ofs << in+ in+ "# " << Pol->Process << ": Set Up Initiation Reaction" << endl;
    ofs << in+ in+ "self.Idx_Pol_" << Pol->Process << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Pol) << "], ndmin=2)" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Pol->Process << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Template) << "], ndmin=2)" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset_" << Pol->Process << " = np.array([range(" << Idx_Template.size() << ")], ndmin=2)   # Default" << endl;
    ofs << in+ in+ "self.Idx_Target_" << Pol->Process << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Target) << "], ndmin=2)" << endl;
    ofs << in+ in+ "self.Weight_" << Pol->Process << " = np.array([" << "1" << "])" << endl;
    ofs << endl;
    ofs << in+ in+ "# for debugging purposes" << endl;
    ofs << in+ in+ "self.Names_Pol_" << Pol->Process << " = [" << Utils::JoinStr2Str(Names_Pol) << "]" << endl;
    ofs << in+ in+ "self.Names_Template_" << Pol->Process << " = [" << Utils::JoinStr2Str(Names_Template) << "]" << endl;
    ofs << in+ in+ "self.Names_Target_" << Pol->Process << " = [" << Utils::JoinStr2Str(Names_Target) << "]" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Pol_Threshold_" << Pol->Process << " = " << std::to_string(Threshold) << endl;
    ofs << endl;

    // Elongation
    ofs << in+ in+ "# " << Pol->Process << ": Set Up Elongation Reaction" << endl;
    ofs << in+ in+ "self.Rate_" << Pol->Process << " = np.array([" << Utils::JoinFloat2Str(Rate) << "])" << endl;

    ofs << in+ in+ "self.Idx_PolSub_" << Pol->Process << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_PolSub) << "], ndmin=2)" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Pol->Process << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_PolBB) << "], ndmin=2)" << endl;
    ofs << endl;
}

void FWriter::Polymerase_InitiationReaction(ofstream& ofs, std::vector<FMolecule*> Polymerases)
{
    FPolymerase* Pol = dynamic_cast<FPolymerase *>(Polymerases[0]);

//    ofs << in+ in+ "# " << Pol->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Pol->TargetClass << "s = self.Initiation(";
    ofs << "self.State.Len_Nascent" << Pol->TemplateClass << "s, ";
    ofs << "self.State.Len_Nascent" << Pol->TargetClass << "s, ";
    ofs << "self.State.Idx_Pol_" << Pol->Process << ", ";
    ofs << "self.State.Idx_Template_" << Pol->Process << ", ";
    ofs << "self.State.Idx_TemplateSubset_" << Pol->Process << ", ";
    ofs << "self.State.Weight_" << Pol->Process << ", ";
    ofs << "self.State.Pol_Threshold_" << Pol->Process << ") ";
    ofs << endl;
}

void FWriter::Polymerase_ElongationReaction(ofstream& ofs, std::vector<FMolecule*> Polymerases)
{
    FPolymerase* Pol = dynamic_cast<FPolymerase *>(Polymerases[0]);

//    ofs << in+ in+ "# " << Pol->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Pol->TargetClass << "s = self.Elongation(";
    ofs << "self.State.Len_Nascent" << Pol->TargetClass << "s, ";
    ofs << "self.State.MaxLen_Nascent" << Pol->TargetClass << "s, ";
    ofs << "self.State.Rate_" << Pol->Process << ", ";
    ofs << "1, ";
    ofs << "self.State.Freq_BB_" << Pol->TargetClass << "s, ";
    ofs << "self.State.Idx_PolSub_" << Pol->Process << ", ";
    ofs << "self.State.Idx_PolBB_" << Pol->Process << ") ";
    ofs << endl;
}

void FWriter::Polymerase_TerminationReaction(ofstream& ofs, std::vector<FMolecule*> Polymerases)
{
    FPolymerase* Pol = dynamic_cast<FPolymerase *>(Polymerases[0]);

//    ofs << in+ in+ "# " << Pol->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Pol->TargetClass << "s = self.Termination(";
    ofs << "self.State.Len_Nascent" << Pol->TargetClass << "s, ";
    ofs << "self.State.MaxLen_Nascent" << Pol->TargetClass << "s, ";
    ofs << "self.State.Idx_Target_" << Pol->Process << ")";
    ofs << endl;
}

void FWriter::Initialize_TransporterReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ "# " << Type << endl;

    ofs << in+ in+ "self.Idx_Cargo_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Idx_Dist_Cargo_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Idx_Transporter_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_ki_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_ko_" << Type << " = None" << endl;

//    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
//    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
//        ofs << in+ in+ "self.Const_K_" << Type << " = None" << endl;
//        ofs << in+ in+ "self.Const_n_" << Type << " = None" << endl;
//        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = None" << endl;
//    }
//    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_TransporterReaction(ofstream& ofs, std::string Type, std::vector<FReaction *> ReactionSubList)
{
    int Idx_Pseudo = Context.GetIdxByName_MoleculeList(Name_Pseudo);
    std::vector<std::string> Name_Cargo;
    std::vector<int> Idx_Cargo;
    std::vector<int> Idx_Dist_Cargo;
    std::vector<int> Idx_Transporter;
    std::vector<float> ki;
    std::vector<float> ko;

//    // for regulatory mechanisms
//    std::vector<float> K;
//    std::vector<float> n;
//    std::vector<int> Idx_Regulator;

    // for stoichiometry matrix
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetUniqueSubstrateIdx_ReactionList(ReactionSubList);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(ReactionSubList);

//    // relevant reaction lists
//    std::vector<FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    for (auto& Reaction : ReactionSubList) {
        const auto& reaction = dynamic_cast<FTransporterReaction *>(Reaction);

        const auto& molecule = Context.GetMolecule_MoleculeList(reaction->Transporter);
        const auto& Transporter = dynamic_cast<FTransporter *>(molecule);

        Idx_Transporter.push_back(Context.GetIdxByName_MoleculeList(Transporter->Name));
        ki.push_back(Transporter->ki);
        ko.push_back(Transporter->ko);

        for (auto& stoich: reaction->Stoichiometry) {
            std::string CargoName = stoich.first;
            Name_Cargo.push_back(CargoName);

            int Idx = Context.GetIdxByName_MoleculeList(CargoName);
            Idx_Cargo.push_back(Idx);
        }
    }

    // get cargo idx in Dist_All
    std::vector<std::string> MolUniqueNames = Context.GetUniqueNames_LocationList("Molecule");
    for (auto& Name : Name_Cargo) {
        int i = 0;
        for (auto& MolUniqueName : MolUniqueNames) {
            if (MolUniqueName == Name) {
                Idx_Dist_Cargo.push_back(i);
            }
            i++;
        }
    }

    ofs << in+ in+ "# " << Type << endl;
    ofs << in+ in+ "self.Idx_Cargo_" << Type << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Cargo) << "])" << endl;
    ofs << in+ in+ "self.Idx_Dist_Cargo_" << Type << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_Dist_Cargo) << "]).transpose()" << endl;
    ofs << in+ in+ "self.Idx_Transporter_" << Type << " = np.array([" << Utils::JoinInt2Str_Idx(Idx_Transporter) << "])" << endl;
    ofs << in+ in+ "self.Const_ki_" << Type << " = np.array([" << Utils::JoinFloat2Str(ki) << "])" << endl;
    ofs << in+ in+ "self.Const_ko_" << Type << " = np.array([" << Utils::JoinFloat2Str(ko) << "])" << endl;

//    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
//    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
//        ofs << in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
//        ofs << in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
//        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
//    }
//    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Utils::Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void FWriter::Initialize_SpatialSimulation(ofstream& ofs)
{
    // TODO: Take dimension from user input
    int Map_Width = 1200;
    int Map_Height = 800;

    ofs << in+ in+ "# Spatial Simulation" << endl;

    ofs << in+ in+ "self.Dimension_X = " << Map_Width << endl;
    ofs << in+ in+ "self.Dimension_Y = " << Map_Height << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");
//    auto Organisms = Context.GetSubList_ContainerList("Organism");

    ofs << in+ in+ "self.Dist_Names = list()" << endl;
    // TODO: update to 3d array
    ofs << in+ in+ "self.Dist_All = np.zeros((" << MolLoc.size() << ", self.Dimension_X, self.Dimension_Y))" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Pos_Names = list()" << endl;
    ofs << in+ in+ "self.Pos_Name2Idx = dict()" << endl;

    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
    for (auto& UniqueName : ObjUniqueNames) {
        ofs << in+ in+ "self.Idx_Pos_" << UniqueName << " = None" << endl;
    }
    ofs << in+ in+ "self.Pos_X = None" << endl;
    ofs << in+ in+ "self.Pos_Y = None" << endl;
    ofs << in+ in+ "self.Pos_Angle = None" << endl;
    ofs << endl;

    ofs << endl;
    if (!Context.ThresholdList.empty()) {
        ofs << in+ in+ "# Threshold-related" << endl;
        ofs << in+ in+ "self.MolNames_Threshold = list()" << endl;
        ofs << in+ in+ "self.Pos_Threshold = None" << endl;

        ofs << in+ in+ "self.Count_Prev = 0" << endl;
        ofs << in+ in+ "self.Idx_Count_Threshold = None" << endl;
        ofs << in+ in+ "self.Idx_Pos_Threshold = None" << endl;
        ofs << endl;
    }
}

void FWriter::SetUp_SpatialSimulation(ofstream& ofs)
{
    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");

    ofs << in+ in+ "self.Dist_Names = [" << Utils::JoinStr2Str(Context.GetNames_LocationList("Molecule")) << "]" << endl;
    for (int i = 0; i < MolLoc.size(); i++) {
        auto& Coord = MolLoc[i]->Coord;
        float MaxAmount = 0;
        float BasalAmount = 0;
        if (Coord[0] != -1) {
            MaxAmount = Context.GetInitialCountByName_CountList(MolLoc[i]->Name);
        } else {
            BasalAmount = Context.GetInitialCountByName_CountList(MolLoc[i]->Name);
        }

//        std::string Shape, Pattern;
//        int Size;

        ofs << in+ in+ "self.Dist_All[" << i << "] = SimF.InitializeDistribution(self.Dimension_X, self.Dimension_Y, "
            << MolLoc[i]->Coord[0] << ", " << MolLoc[i]->Coord[1] << ", MaxAmount=" << MaxAmount << ", BasalAmount=" << BasalAmount;
//        if (!Shape.empty()) {
//            ofs << ", " << Shape << ", " << Size << ", " << Pattern;
//            }
        ofs << ")" << endl;
    }
    ofs << endl;

    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");

    ofs << in+ in+ "self.Pos_Names = [" << Utils::JoinStr2Str(ObjUniqueNames) << "]" << endl;
    int i = 0;
    for (auto& ObjName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(ObjName));
        ofs << in+ in+ "self.Idx_Pos_" << ObjName << " = np.array([";
        for (int j = 0; j < (Count); j++) {
            ofs << i << ", ";
            i++;
        }
        ofs << "])" << endl;

        // TODO: may not be used
        ofs << in+ in+ "self.Pos_Name2Idx['" << ObjName << "'] = self.Idx_Pos_" << ObjName << endl;
    }
    ofs << endl;
// here
    ofs << in+ in+ "# Currently support X, Y, Angle, Threshold" << endl;
    ofs << in+ in+  "self.Pos_X = np.array([";
    i = 0;
    for (auto& ObjName : ObjUniqueNames) {
        for (auto& objloc : ObjLoc){
            if (objloc->Name == ObjName) {
                int Count = int(objloc->Count->Amount);
                for (int j = 0; j < (Count); j++) {
                    auto& Coord = objloc->Coord;
                    ofs << Coord[0] << ", ";
                    i++;
                }
            }
        }
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Y = np.array([";
    i = 0;
    for (auto& ObjName : ObjUniqueNames) {
        for (auto& objloc : ObjLoc){
            if (objloc->Name == ObjName) {
                int Count = int(objloc->Count->Amount);
                for (int j = 0; j < (Count); j++) {
                    auto& Coord = objloc->Coord;
                    ofs << Coord[1] << ", ";
                    i++;
                }
            }
        }
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Angle = np.array([";
    i = 0;
    for (auto& ObjName : ObjUniqueNames) {
        for (auto& objloc : ObjLoc){
            if (objloc->Name == ObjName) {
                int Count = int(objloc->Count->Amount);
                for (int j = 0; j < (Count); j++) {
                    auto& Coord = objloc->Coord;
                    ofs << "0, ";
                    i++;
                }
            }
        }
    }
    ofs << "]) " << endl;
    ofs << endl;

    if (!Context.ThresholdList.empty()) {

        std::vector<std::string> Names_ThresholdedMolecules;
        std::vector<int> Idx_Count_Threshold;
        std::vector<float> ThresholdValues;

        for (auto& Threshold : Context.ThresholdList) {
            Names_ThresholdedMolecules.push_back(Threshold.first);
            int Idx = Context.GetIdxByName_MoleculeList(Threshold.first);
            Idx_Count_Threshold.push_back(Idx);
            ThresholdValues.push_back(Threshold.second);
        }

        int MatrixSize = 1;
        if (!ObjLoc.empty()) {
            MatrixSize = Context.GetCounts_LocationList("Compartment");
        }

        ofs << in+ in+ "# Threshold-related" << endl;
        ofs << in+ in+ "self.MolNames_Threshold = [" << Utils::JoinStr2Str(Names_ThresholdedMolecules) << "]" << endl;
        ofs << in+ in+ "self.Pos_Threshold = np.tile([" << Utils::JoinFloat2Str(ThresholdValues) << "], [" << MatrixSize << ", 1]).transpose()" << endl;
        ofs << endl;

        // threshold index setting
        ofs << in+ in+ "self.Idx_Count_Threshold = np.array([" << Utils::JoinInt2Str_Idx(Idx_Count_Threshold) <<  "])" << endl;
        ofs << in+ in+ "self.Idx_Pos_Threshold = np.array(range(len(self.MolNames_Threshold)))" << endl;
        ofs << endl;

    }
}
