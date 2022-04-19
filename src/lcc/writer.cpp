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
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(ReactionSubList);
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
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(ReactionSubList);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(ReactionSubList);

    // relevant reaction lists
    std::vector<FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    if (NameSpace_Pathway != "") {
        RegulatoryReactionSubList = Context.FilterByPathway_ReactionList(RegulatoryReactionSubList, NameSpace_Pathway);
    }

    for (auto& Reaction : ReactionSubList) {
        const auto& EnzReaction = dynamic_cast<FEnzymaticReaction *>(Reaction);
        const auto& enzyme = Context.GetEnzyme_EnzymeList(EnzReaction->Enzyme);

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

void FWriter::Initialize_PolymeraseReaction(ofstream& ofs, FPolymerase* Polymerase)
{
    ofs << in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = None" << endl;

    // Initiation
    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset" << Polymerase->Process << " = None" << endl; // local indexing within the template population for mRNA in RNA for protein translation
    ofs << in+ in+ "self.Weight_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Pol_Threshold_" << Polymerase->Process << " = None" << endl;
    ofs << endl;

    // Check if template exists as a target of any polymerases in the system
    bool TemplateExists = false;
//    auto& PolymeraseList = Context.GetList_PolymeraseMoleculeList
    for (auto& PolymeraseInSystem : Context.GetList_Polymerase_MoleculeList()) {
        if (PolymeraseInSystem->Target == Polymerase->Template) {
            TemplateExists = true;
            break;
        }
    }

    if (!TemplateExists) {
        ofs << in+ in+ "self.Len_Nascent" << Polymerase->Template << "s = None" << endl;
    }

    // Elongation
    ofs << in+ in+ "# " << Polymerase->Process << " Initialize ElongationReaction" << endl;
    ofs << in+ in+ "self.Rate_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.MaxLen_Nascent" << Polymerase->Target << "s = None" << endl;
    ofs << in+ in+ "self.Len_Nascent" << Polymerase->Target << "s = None" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = None" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_PolymeraseReaction(ofstream& ofs, FPolymerase* Polymerase, float Rate, std::string FreqBBFileName, std::string MaxLenFileName, int Idx_Pol, std::vector<int> Idx_Template, std::vector<int> Idx_TemplateSubset, std::vector<int> Idx_Target, std::vector<int> Idx_PolSub, std::vector<int> Idx_PolBB, int Threshold)
{
    ofs << in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << FreqBBFileName << "'))" << endl;

    // Initiation
    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Idx_Pol) << "])" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Polymerase->Process << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_Template) << "])" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset_" << Polymerase->Process << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_TemplateSubset) << "])" << endl;
    ofs << in+ in+ "self.Idx_Target_" << Polymerase->Process << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_Target) << "])" << endl;
    ofs << in+ in+ "self.Weight_" << Polymerase->Process << " = np.array([" << "1" << "])" << endl;
    ofs << endl;

    // Check if template exists as a target of any polymerases in the system
    bool TemplateExists = false;
    for (auto& PolymeraseInSystem : Context.GetList_Polymerase_MoleculeList()){
        if (PolymeraseInSystem->Target == Polymerase->Template){
            TemplateExists = true;
            break;
        }
    }

    if (!TemplateExists) {
        ofs << in+ in+ "self.Len_Nascent" << Polymerase->Template << "s = np.asmatrix(np.full([10, self.Freq_BB_" << Polymerase->Target << "s.shape[0]], -1))" << endl;
        ofs << endl;
    }

    int Template_Min, Template_Max, Target_Min, Target_Max;

    if      (Polymerase->Process == "DNAReplication") 	        {Template_Min = 1; Template_Max = 2; Target_Min = 1; Target_Max = 2;}
    else if (Polymerase->Process == "RNATranscription") 	{Template_Min = 1; Template_Max = 2; Target_Min = 0; Target_Max = 1;}
    else if (Polymerase->Process == "ProteinTranslation") 	{Template_Min = 0; Template_Max = 1; Target_Min = 0; Target_Max = 1;}
    else    				                        {Template_Min = 0; Template_Max = 1; Target_Min = 0; Target_Max = 1;} // improve exception handling


    ofs << in+ in+ "self.Pol_Threshold_" << Polymerase->Process << " = " << std::to_string(Threshold) << endl;
    ofs << endl;
    ofs << in+ in+ "Count_Template_" << Polymerase->Process << " = np.random.randint(" << std::to_string(Template_Min) << ", high=" << std::to_string(Template_Max) << ", size=self.Idx_Template_" << Polymerase->Process << ".shape[1])" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Template_" << Polymerase->Process << ", Count_Template_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_Target_" << Polymerase->Process << " = np.random.randint(" << std::to_string(Target_Min) << ", high=" << std::to_string(Target_Max) << ", size=self.Idx_Target_" << Polymerase->Process << ".shape[1])" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Target_" << Polymerase->Process << ", Count_Target_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;


    // Elongation
    ofs << in+ in+ "# " << Polymerase->Process << " SetUp ElongationReaction" << endl;
    ofs << in+ in+ "self.Rate_" << Polymerase->Process << " = " << std::to_string(int(Rate)) << endl;
    ofs << in+ in+ "self.MaxLen_Nascent" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << MaxLenFileName << "'))" << endl;
    ofs << in+ in+ "self.Len_Nascent" << Polymerase->Target << "s = np.asmatrix(np.full([10, self.Freq_BB_" << Polymerase->Target << "s.shape[0]], -1))" << endl;
    ofs << endl;

    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Context.GetIdxByName_MoleculeList(Polymerase->Name)) << "])" << endl;
    ofs << in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_PolSub) << "])" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = np.asmatrix([" << Utils::JoinInt2Str_Idx(Idx_PolBB) << "])" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_Pol_" << Polymerase->Process << " = np.asmatrix(np.full(len(self.Idx_Pol_" << Polymerase->Process << "), 100))" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_Pol_" << Polymerase->Process << ", Count_Pol_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_PolSub_" << Polymerase->Process << " = np.asmatrix(np.full(len(self.Idx_PolSub_" << Polymerase->Process << "), 1000))" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_PolSub_" << Polymerase->Process << ", Count_PolSub_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

    ofs << in+ in+ "Count_PolBB_" << Polymerase->Process << " = np.random.randint(4000, high=5000, size=self.Freq_BB_" << Polymerase->Target << "s.shape[1])" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, self.Idx_PolBB_" << Polymerase->Process << ", Count_PolBB_" << Polymerase->Process << ", axis=1)" << endl;
    ofs << endl;

}

void FWriter::Polymerase_InitiationReaction(ofstream& ofs, FPolymerase* Polymerase)
{
    ofs << in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Initiation(";
    ofs << "self.State.Len_Nascent" << Polymerase->Template << "s, ";
    ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
    ofs << "self.State.Idx_Pol_" << Polymerase->Process << ", ";
    ofs << "self.State.Idx_Template_" << Polymerase->Process << ", ";
    ofs << "self.State.Idx_TemplateSubset_" << Polymerase->Process << ", ";
    ofs << "self.State.Weight_" << Polymerase->Process << ", ";
    ofs << "self.State.Pol_Threshold_" << Polymerase->Process << ") " << endl;
    ofs << endl;
}

void FWriter::Polymerase_ElongationReaction(ofstream& ofs, FPolymerase* Polymerase)
{
    ofs << in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Elongation(";
    ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
    ofs << "self.State.MaxLen_Nascent" << Polymerase->Target << "s, ";
    ofs << "self.State.Rate_" << Polymerase->Process << ", ";
    ofs << "1, ";
    ofs << "self.State.Freq_BB_" << Polymerase->Target << "s, ";
    ofs << "self.State.Idx_PolSub_" << Polymerase->Process << ", ";
    ofs << "self.State.Idx_PolBB_" << Polymerase->Process << ") " << endl;
    ofs << endl;
}

void FWriter::Polymerase_TerminationReaction(ofstream& ofs, FPolymerase* Polymerase)
{
    ofs << in+ in+ "# " << Polymerase->Process << endl;
    ofs << in+ in+ "self.State.Len_Nascent" << Polymerase->Target << "s = self.Termination(";
    ofs << "self.State.Len_Nascent" << Polymerase->Target << "s, ";
    ofs << "self.State.MaxLen_Nascent" << Polymerase->Target << "s, ";
    ofs << "self.State.Idx_Target_" << Polymerase->Process << ")" << endl;
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
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(ReactionSubList);
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
    int Map_Width = 1200;
    int Map_Height = 800;

    ofs << in+ in+ "# Spatial Simulation" << endl;

    ofs << in+ in+ "self.Dimension_X = " << Map_Width << endl;
    ofs << in+ in+ "self.Dimension_Y = " << Map_Height << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");
    auto& Organisms = Context.OrganismList;

    ofs << in+ in+ "self.Dist_Names = list()" << endl;
    // TODO: update to 3d array
    ofs << in+ in+ "self.Dist_All = np.zeros((" << MolLoc.size() << ", self.Dimension_X , self.Dimension_Y))" << endl;
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

    ofs << in+ in+ "self.Pos_Names = [" << Utils::JoinStr2Str(Context.GetNames_LocationList("Compartment")) << "]" << endl;

    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
    for (int i = 0; i < ObjUniqueNames.size(); i++) {
        int Count = int(Context.GetInitialCountByName_CountList(ObjUniqueNames[i]));
        ofs << in+ in+ "self.Idx_Pos_" << ObjUniqueNames[i] << " = np.array([";
        for (int j = i; j < (Count); j++) {
            ofs << j << ", ";
        }
        ofs << "])" << endl;

        // TODO: may not be used
        ofs << in+ in+ "self.Pos_Name2Idx['" << ObjUniqueNames[i] << "'] = self.Idx_Pos_" << ObjUniqueNames[i] << endl;
    }
    ofs << endl;

    ofs << in+ in+ "# Currently support X, Y, Angle, Threshold" << endl;
    ofs << in+ in+  "self.Pos_X = np.array([";
    for (int i = 0; i < ObjUniqueNames.size(); i++) {
        int Count = int(Context.GetInitialCountByName_CountList(ObjUniqueNames[i]));
        for (int j = i; j < (Count); j++) {
            auto& Coord = ObjLoc[j]->Coord;
            ofs << Coord[0] << ", ";
        }
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Y = np.array([";
    for (int i = 0; i < ObjUniqueNames.size(); i++) {
        int Count = int(Context.GetInitialCountByName_CountList(ObjUniqueNames[i]));
        for (int j = i; j < (Count); j++) {
            auto& Coord = ObjLoc[j]->Coord;
            ofs << Coord[1] << ", ";
        }
    }
    ofs << "]) " << endl;

    ofs << in+ in+  "self.Pos_Angle = np.array([";
    for (int i = 0; i < ObjUniqueNames.size(); i++) {
        int Count = int(Context.GetInitialCountByName_CountList(ObjUniqueNames[i]));
        for (int j = i; j < (Count); j++) {
            ofs << "0, ";
//            ofs << Numbers::RandomNumber(0, 1) << " * 2 * np.pi, ";
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

        ofs << in+ in+ "# Threshold-related" << endl;
        ofs << in+ in+ "self.MolNames_Threshold = [" << Utils::JoinStr2Str(Names_ThresholdedMolecules) << "]" << endl;;
        ofs << in+ in+ "self.Pos_Threshold = np.tile([" << Utils::JoinFloat2Str(ThresholdValues) << "], [ " << ObjLoc.size() << ", 1]).transpose()" << endl;
        ofs << endl;

        // threshold index setting
        ofs << in+ in+ "self.Idx_Count_Threshold = np.array([" << Utils::JoinInt2Str_Idx(Idx_Count_Threshold) <<  "])" << endl;
        ofs << in+ in+ "self.Idx_Pos_Threshold = np.array(range(len(self.MolNames_Threshold)))" << endl;
        ofs << endl;

    }
}
