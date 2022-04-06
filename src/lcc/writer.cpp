#include "writer.h"
#include <algorithm>

using namespace std;

// Utils for string expression
std::string JoinStr2Str(std::vector<std::string> StringList)
{ return Utils::JoinStr2Str(StringList); }

std::string JoinInt2Str(std::vector<int> IntList)
{ return Utils::JoinInt2Str(IntList); }

std::string JoinInt2Str_Idx(std::vector<int> IntList)
{ return Utils::JoinInt2Str_Idx(IntList); }

std::string JoinFloat2Str(std::vector<float> FloatList)
{ return Utils::JoinFloat2Str(FloatList); }

std::string Matrix2Str(std::vector<std::vector<int>> Matrix)
{ return Utils::Matrix2Str(Matrix); }

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

void FWriter::Initialize_StandardReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ "# " << Type << endl;

    // Standard vs. MichaelisMenten
    ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_k_Product_" << Type << " = None" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = None" << endl;
    }
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = None" << endl;
    }
    ofs << endl;

    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_StandardReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList)
{
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
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type); // TODO: Add new type
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type); // TODO: Add new type

    // relevant reaction lists
    std::vector<const FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    for (auto Reaction : ReactionSubList) {
        const auto& reaction = dynamic_cast<const FStandardReaction *>(Reaction);

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
            for (auto &Reg: RegulatoryReactionSubList) {
                const auto& reg = dynamic_cast<const FRegulatoryReaction *>(Reg);

                bool Import = false;
                std::string reactant_reg;
                std::string product_reg;
                for (auto &stoich: reg->Stoichiometry) {
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

    ofs << in+ in+ "# " << Type << endl;

    ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = np.array([" << JoinFloat2Str(k) << "])" << endl;
    ofs << in+ in+ "self.Const_k_Product_" << Type << " = np.array([" << JoinFloat2Str(krev) << "])" << endl;
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
    }
    for (int i = 0; i < N_MoleculesAllowed; i++) {
        ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void FWriter::Initialize_EnzymeReaction(ofstream& ofs, std::string Type)
{
    ofs << in+ in+ "# " << Type << endl;

    ofs << in+ in+ "self.Idx_Enz_" << Type << " = None" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Enz_Standard") != string::npos) {
        ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_k_Product_" << Type << " = None" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = None" << endl;
        }
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = None" << endl;
        }

    } else if (Type.find("Enz_MichaelisMenten") != string::npos) {
        ofs << in+ in+ "self.Const_kcat_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_KM_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Idx_EnzSub_" << Type << " = None" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = None" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = None" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = None" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_EnzymeReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList) // to be changed with reaction list
{
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
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type);

    // relevant reaction lists
    std::vector<const FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    for (auto& Reaction : ReactionSubList) {
        const auto& reaction = dynamic_cast<const FEnzymaticReaction *>(Reaction);
        const auto& enzyme = Context.GetEnzyme_EnzymeList(reaction->Enzyme);

        Idx_Enz.push_back(Context.GetIdxByName_MoleculeList(enzyme->Name));

        // Enz_Standard
        if ((reaction->Type >= 10) & (reaction->Type < 20)) {
            const auto& reaction = dynamic_cast<const FEnz_StandardReaction *>(Reaction);

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
        } else if ((reaction->Type >= 20) & (reaction->Type < 30)) {
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
            for (auto &Reg: RegulatoryReactionSubList) {
                const auto &reg = dynamic_cast<const FRegulatoryReaction *>(Reg);

                bool Import = false;
                std::string reactant_reg;
                std::string product_reg;
                for (auto &stoich: reg->Stoichiometry) {
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

    // Check Idx lengths
    Utils::Assertion((Idx_Reactants.size() == Idx_Products.size()), "# of parsed reactants and products do not match");

    ofs << in+ in+ "# " << Type << endl;
    ofs << in+ in+ "self.Idx_Enz_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Enz) << "])" << endl;

    // Standard vs. MichaelisMenten
    if (Type.find("Enz_Standard") != string::npos) {
        ofs << in+ in+ "self.Const_k_Reactant_" << Type << " = np.array([" << JoinFloat2Str(k) << "])" << endl;
        ofs << in+ in+ "self.Const_k_Product_" << Type << " = np.array([" << JoinFloat2Str(krev) << "])" << endl;
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Reactant_" << i << "_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Reactants[i]) << "])" << endl;
        }
        for (int i = 0; i < N_MoleculesAllowed; i++) {
            ofs << in+ in+ "self.Idx_Product_" << i << "_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Products[i]) << "])" << endl;
        }

    } else if (Type.find("Enz_MichaelisMenten") != string::npos) {
        ofs << in+ in+ "self.Const_kcat_" << Type << " = np.array([" << JoinFloat2Str(kcats) << "])" << endl;
        ofs << in+ in+ "self.Const_KM_" << Type << " = np.array([" << JoinFloat2Str(KMs) << "])" << endl;
        ofs << in+ in+ "self.Idx_EnzSub_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_EnzSub) << "])" << endl;
    }

    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
        ofs << in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
        ofs << in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
    ofs << endl;
}

void FWriter::Initialize_PolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase)
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

void FWriter::SetUp_PolymeraseReaction(ofstream& ofs, const FPolymerase* Polymerase, float Rate, std::string FreqBBFileName, std::string MaxLenFileName, int Idx_Pol, std::vector<int> Idx_Template, std::vector<int> Idx_TemplateSubset, std::vector<int> Idx_Target, std::vector<int> Idx_PolSub, std::vector<int> Idx_PolBB, int Threshold)
{
    ofs << in+ in+ "self.Freq_BB_" << Polymerase->Target << "s = np.asmatrix(np.load(r'" << FreqBBFileName << "'))" << endl;

    // Initiation
    ofs << in+ in+ "self.Idx_Pol_" << Polymerase->Process << " = np.asmatrix([" << std::to_string(Idx_Pol) << "])" << endl;
    ofs << in+ in+ "self.Idx_Template_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Template) << "])" << endl;
    ofs << in+ in+ "self.Idx_TemplateSubset_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_TemplateSubset) << "])" << endl;
    ofs << in+ in+ "self.Idx_Target_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Target) << "])" << endl;
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
    ofs << in+ in+ "self.Idx_PolSub_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_PolSub) << "])" << endl;
    ofs << in+ in+ "self.Idx_PolBB_" << Polymerase->Process << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_PolBB) << "])" << endl;
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

void FWriter::Polymerase_InitiationReaction(ofstream& ofs, const FPolymerase* Polymerase)
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

void FWriter::Polymerase_ElongationReaction(ofstream& ofs, const FPolymerase* Polymerase)
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

void FWriter::Polymerase_TerminationReaction(ofstream& ofs, const FPolymerase* Polymerase)
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

void FWriter::SetUp_TransporterReaction(ofstream& ofs, std::string Type, std::vector<const FReaction *> ReactionSubList)
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
    std::vector<int> Idx_Mol_InStoichMatrix = Context.GetIdxForStoichiometryMatrix(Type);
    std::vector<std::vector<int>> StoichMatrix = Context.GetStoichiometryMatrix(Type);

//    // relevant reaction lists
//    std::vector<const FReaction *> RegulatoryReactionSubList = Context.GetSubList_ReactionList(RegType);

    for (auto& Reaction : ReactionSubList) {
        const auto &reaction = dynamic_cast<const FTransporterReaction *>(Reaction);

        const auto& molecule = Context.GetMolecule_MoleculeList(reaction->Transporter);
        const auto& Transporter = dynamic_cast<const FTransporter *>(molecule);

        Idx_Transporter.push_back(Context.GetIdxByName_MoleculeList(Transporter->Name));
        ki.push_back(Transporter->ki);
        ko.push_back(Transporter->ko);

        for (auto &stoich: reaction->Stoichiometry) {
            std::string CargoName = stoich.first;
            Name_Cargo.push_back(CargoName);

            int Idx = Context.GetIdxByName_MoleculeList(CargoName);
            Idx_Cargo.push_back(Idx);
        }
    }

    // get cargo idx in Dist_All
    std::vector<std::string> MolUniqueNames = Context.GetUniqueNames_LocationList("Molecule");
    for (auto Name : Name_Cargo) {
        int i = 0;
        for (auto MolUniqueName : MolUniqueNames) {
            if (MolUniqueName == Name) {
                Idx_Dist_Cargo.push_back(i);
            }
            i++;
        }
    }

    ofs << in+ in+ "# " << Type << endl;
    ofs << in+ in+ "self.Idx_Cargo_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Cargo) << "])" << endl;
    ofs << in+ in+ "self.Idx_Dist_Cargo_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Dist_Cargo) << "]).transpose()" << endl;
    ofs << in+ in+ "self.Idx_Transporter_" << Type << " = np.array([" << JoinInt2Str_Idx(Idx_Transporter) << "])" << endl;
    ofs << in+ in+ "self.Const_ki_" << Type << " = np.array([" << JoinFloat2Str(ki) << "])" << endl;
    ofs << in+ in+ "self.Const_ko_" << Type << " = np.array([" << JoinFloat2Str(ko) << "])" << endl;

//    // Inhibition vs. Activation (improve with number of accepted regulators implemented)
//    if ((Type.find("Inhibition") != string::npos) || (Type.find("Activation") != string::npos)) {
//        ofs << in+ in+ "self.Const_K_" << Type << " = np.array([" << JoinFloat2Str(K) << "])" << endl;
//        ofs << in+ in+ "self.Const_n_" << Type << " = np.array([" << JoinFloat2Str(n) << "])" << endl;
//        ofs << in+ in+ "self.Idx_Regulator_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Regulator) << "])" << endl;
//    }
//    ofs << endl;

    ofs << in+ in+ "self.Idx_Mol_InStoichMatrix_" << Type << " = np.asmatrix([" << JoinInt2Str_Idx(Idx_Mol_InStoichMatrix) << "])" << endl;
    ofs << in+ in+ "self.Const_StoichMatrix_" << Type << " = np.asmatrix([" << Matrix2Str(StoichMatrix) << "])" << endl;
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
    auto Organisms = Context.OrganismList;

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
    ofs << in+ in+ "self.Pos_Threshold = None" << endl;
    ofs << endl;
}

void FWriter::SetUp_SpatialSimulation(ofstream& ofs)
{
    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");

    ofs << in+ in+ "self.Dist_Names = [" << JoinStr2Str(Context.GetNames_LocationList("Molecule")) << "]" << endl;
    for (int i = 0; i < MolLoc.size(); i++) {
        auto Coord = MolLoc[i]->Coord;
        auto Amount = Context.GetInitialCountByName_CountList(MolLoc[i]->Name);
        std::string Shape, Pattern;
        int Size;

        ofs << in + in + "self.Dist_All[" << i << "] = SimF.InitializeDistribution(self.Dimension_X, self.Dimension_Y, "
            << MolLoc[i]->Coord[0] << ", " << MolLoc[i]->Coord[1] << ", " << Amount;
        if (!Shape.empty()) {
            ofs << ", " << Shape << ", " << Size << ", " << Pattern;
            }
        ofs << ")" << endl;
    }
    ofs << endl;

    ofs << in+ in+ "self.Pos_Names = [" << JoinStr2Str(Context.GetNames_LocationList("Compartment")) << "]" << endl;

    int i = 0;
    std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
    for (auto& UniqueName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
        ofs << in+ in+ "self.Idx_Pos_" << UniqueName << " = np.array([";
        for (int j = i; j < (Count); j++) {
            ofs << j << ", ";
        }
        ofs << "])" << endl;

        // TODO: may not be used
        ofs << in+ in+ "self.Pos_Name2Idx['" << UniqueName << "'] = self.Idx_Pos_" << UniqueName << endl;
        i++;
    }
    ofs << endl;

    i = 0;
    ofs << in+ in+ "# Currently support X, Y, Angle, Threshold" << endl;
    ofs << in+ in+  "self.Pos_X = np.array([";
    for (auto& UniqueName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
        for (int j = i; j < (Count); j++) {
            auto Coord = ObjLoc[j]->Coord;
            ofs << Coord[0] << ", ";
        }
    }
    ofs << "]) " << endl;

    i = 0;
    ofs << in+ in+  "self.Pos_Y = np.array([";
    for (auto& UniqueName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
        for (int j = i; j < (Count); j++) {
            auto Coord = ObjLoc[j]->Coord;
            ofs << Coord[1] << ", ";
        }
    }
    ofs << "]) " << endl;

    i = 0;
    ofs << in+ in+  "self.Pos_Angle = np.array([";
    for (auto& UniqueName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
        for (int j = i; j < (Count); j++) {
            ofs << Numbers::RandomNumber(0, 1) << " * 2 * np.pi, ";
        }
    }
    ofs << "]) " << endl;

    i = 0;
    ofs << in+ in+  "self.Pos_Threshold = np.array([";
    for (auto& UniqueName : ObjUniqueNames) {
        int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
        for (int j = i; j < (Count); j++) {
            ofs << "0.0, ";
        }
//        ofs << Numbers::MultiplyByAvogadro(0.983405e-9) << ", ";
    }
    ofs << "]) " << endl;
    ofs << endl;
}

void FWriter::SimIdx() {}

void FWriter::SimModule(int Sim_Steps, int Sim_Resolution)
{
    std::cout << "Generating SimModule..." << std::endl;

    // write SimModule.py
    std::ofstream ofs(Option.SimModuleFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << "import numpy as np" << endl;
    // ofs << "import tensorflow as tf" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import csv" << endl;
    ofs << "import SimFunctions as SimF" << endl;
    // ofs << "import SimIdx as idx" << endl;
    ofs << "import plot" << endl;
    ofs << "from argparse import ArgumentParser" << endl;
    ofs << endl;

    //TODO: Take options from SimModule cmd line
    ofs << "# Temporary global variables" << endl;
    ofs << "N_SimSteps = " << Sim_Steps << endl;
    ofs << "SimStepTimeResolution = " << Sim_Resolution << endl;
    ofs << "DisplayCount = 0   # 0: Display Counts only, 1: Display both Counts & dCounts" << endl;
    ofs << "Unit = 'nM'   # nM or uM supported for now" << endl;
    ofs << endl;

    ofs << "def ListToStr(List, Delimiter):" << endl;
    ofs << in+ "List_Str = ''" << endl;
    ofs << in+ "if List:" << endl;
    ofs << in+ in+ "for item in List:" << endl;
    ofs << in+ in+ in+ "List_Str += item + '_'" << endl;
    ofs << in+ "return List_Str" << endl;
    ofs << endl;

    // class FState
    ofs << "class FState:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.Vol = 0" << endl;
    ofs << endl;

    ofs << in+ in+ "# State Arrays" << endl;

    ofs << in+ in+ "# Temporary legend attribute" << endl;
    ofs << in+ in+ "self.Mol_Names = list()" << endl;
    ofs << in+ in+ "self.Mol_Name2Idx = dict()" << endl;
    ofs << in+ in+ "self.Legends = list()" << endl;
    ofs << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");
    auto ObjLoc = Context.GetSubList_LocationList("Compartment");

    int MatrixSize = 1;
    if (!ObjLoc.empty()) {
        MatrixSize = 0;
        std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
        for (auto UniqueName : ObjUniqueNames) {
            int Count = int(Context.GetInitialCountByName_CountList(UniqueName));
            MatrixSize += Count;
        }
    }
    ofs << in+ in+ "self.Count_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << in+ in+ "self.dCount_All = np.zeros([" << MatrixSize << ", " << Context.MoleculeList.size() << "])" << endl;
    ofs << endl;

    if (!Context.LocationList.empty()) {
        Initialize_SpatialSimulation(ofs);
    }

    // for standard reactions
    std::vector<std::string> ReactionTypes;
    ReactionTypes.emplace_back("Standard_Unregulated");
    ReactionTypes.emplace_back("Standard_Inhibition_Allosteric");
    ReactionTypes.emplace_back("Standard_Activation_Allosteric");

    std::vector<std::string> StandardReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_StandardReaction(ofs, Type);
        }
    }

    // for Enz reactions
    ReactionTypes.clear();
    ReactionTypes.push_back("Enz_Standard_Unregulated");
    ReactionTypes.push_back("Enz_Standard_Inhibition_Allosteric");
    ReactionTypes.push_back("Enz_Standard_Activation_Allosteric");
    ReactionTypes.push_back("Enz_MichaelisMenten_Unregulated");
    ReactionTypes.push_back("Enz_MichaelisMenten_Inhibition_Allosteric");
    ReactionTypes.push_back("Enz_MichaelisMenten_Inhibition_Competitive");
    ReactionTypes.push_back("Enz_MichaelisMenten_Activation_Allosteric");

    std::vector<std::string> EnzReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_EnzymeReaction(ofs, Type);
        }
    }

    // for polymerase reactions (Template-based)
    std::vector<const FPolymerase *> PolymeraseList = Context.GetList_Polymerase_MoleculeList();
//    std::vector<std::string> PolymeraseNames = Context.GetNames_PolymeraseList(PolymeraseList);
    std::vector<const FPolymeraseReaction *> PolymeraseReactionList = Context.GetList_Polymerase_ReactionList();

    if (!PolymeraseList.empty()) {
        for (auto& Polymerase : PolymeraseList) {
            Initialize_PolymeraseReaction(ofs, Polymerase);
        }
    }

    // for Transporter reactions
    ReactionTypes.clear();
    ReactionTypes.push_back("Transporter_Unregulated");
    ReactionTypes.push_back("Transporter_Inhibition_Allosteric");
    ReactionTypes.push_back("Transporter_Inhibition_Competitive");
    ReactionTypes.push_back("Transporter_Activation_Allosteric");

    std::vector<std::string> TransporterReactionTypes = ReactionTypes; // to reuse later

    for (auto& Type : ReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            Initialize_TransporterReaction(ofs, Type);
        }
    }

    ofs << in+ "def Initialize(self):" << endl;
    ofs << in+ in+ "self.Vol = 1" << endl;
    ofs << endl;

    // Print SetUp_SpatialReaction for all spatial simulation (to be updated)
    if (!Context.LocationList.empty()) {
        SetUp_SpatialSimulation(ofs);
    }

    // Print SetUp_StandardReaction for each Reaction Type
    for (auto& Type : StandardReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_StandardReaction(ofs, Type, ReactionSubList);
        }
    }

    // Print SetUp_EnzymeReaction for each Reaction Type
    for (auto& Type : EnzReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_EnzymeReaction(ofs, Type, ReactionSubList);
        }
    }

    if (!PolymeraseList.empty()) {
        for (auto& Polymerase : PolymeraseList) {
            std::string PolymeraseName = Polymerase->Name;
            float Rate = Polymerase->Rate;
            int Idx_Pol = Context.GetIdxByName_MoleculeList(PolymeraseName);
            std::vector<int> Idx_Template;
            std::vector<int> Idx_Target;
            std::vector<int> Idx_TemplateSubset;
            std::vector<int> Idx_PolSub = Context.GetIdx_PolymeraseReactionSubstrate_ByPolymeraseName_MoleculeList(PolymeraseName);
            std::vector<int> Idx_PolBB;
            std::vector<std::string> BuildingBlockNames;
            std::string MaxLenFileName;
            std::string FreqBBFileName;
            int Threshold;

            for (auto& reaction : Context.GetList_Polymerase_ReactionList()){
                if (reaction->Name == PolymeraseName){
                    for (auto& BuildingBlock : reaction->BuildingBlocks){
                        BuildingBlockNames.push_back(BuildingBlock);
                    }
                    Idx_PolBB = Context.GetIdxByStrList_MoleculeList(BuildingBlockNames);
                }
            }

            std::vector<std::vector<int>> StoichMatrix_PolymeraseReaction = Context.GetStoichiometryMatrix_PolymeraseReaction(PolymeraseReactionList);

            if (Polymerase->Process == "DNAReplication") {
                MaxLenFileName = "./Database/Len_ChromosomesInGenome.npy";
                FreqBBFileName = "./Database/Freq_NTsInChromosomesInGenome.npy";
                Idx_Template   = Context.GetIdxListFromMoleculeList("Chromosome");
                Idx_Target     = Context.GetIdxListFromMoleculeList("Chromosome");
                Threshold      = 50;  // The number of polymerases required for a functional unit to initiate the polymerase reaction

                // Idx_TemplateSubset
                for (int i = 0; i < Idx_Template.size(); i++) {
                    Idx_TemplateSubset.push_back(i);
                }
            }

            else if (Polymerase->Process == "RNATranscription") {
                MaxLenFileName = "./Database/Len_RNAs.npy";
                FreqBBFileName = "./Database/Freq_NTsInRNAs.npy";
                Idx_Template   = Context.GetIdxListFromMoleculeList("Gene");
                Idx_Target     = Context.GetIdxListFromMoleculeList("RNA"); // update to link gene to RNA matching when importing and generating this list
                Threshold      = 1;

                // Idx_TemplateSubset
                for (int i = 0; i < Idx_Template.size(); i++) {
                    Idx_TemplateSubset.push_back(i);
                }
            }

            else if (Polymerase->Process == "ProteinTranslation") {
                MaxLenFileName = "./Database/Len_Proteins.npy";
                FreqBBFileName = "./Database/Freq_AAsInProteins.npy";
                Idx_Template   = Context.GetIdxListFromMoleculeList("mRNA");
                Idx_Target     = Context.GetIdxListFromMoleculeList("Protein"); // update to link mRNA to protein matching
                Threshold      = 1;

                Idx_TemplateSubset = Context.GetIdxOfStrListFromStrList(Context.GetNameListFromMoleculeList("mRNA"), Context.GetNameListFromMoleculeList("RNA"));
            }

            else { // add exception handling
            }

            // debugging
            std::cout << Polymerase->Process;
            if (Polymerase->Process == "DNAReplication") {
                std::cout << "\t";
            }
            std::cout << "\t | Idx_Template.size(): " << Idx_Template.size() << "\t Idx_Target.size(): " << Idx_Target.size() << endl;
            Utils::Assertion((Idx_Template.size() == Idx_Target.size()), "# of Target and Target do not match for Polymerase: " + PolymeraseName);

            SetUp_PolymeraseReaction(ofs, Polymerase, Rate, FreqBBFileName, MaxLenFileName, Idx_Pol, Idx_Template, Idx_TemplateSubset, Idx_Target, Idx_PolSub, Idx_PolBB, Threshold);
        }
    }

    // Print SetUp_TransporterReaction for each Reaction Type
    for (auto& Type : TransporterReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            SetUp_TransporterReaction(ofs, Type, ReactionSubList);
        }
    }

    // for legends
    std::vector<std::string> MolNames = Context.GetNames_MoleculeList();
    ofs << in+ in+ "self.Mol_Names = [" << JoinStr2Str(MolNames) << "]" << endl;
    ofs << in+ in+ "self.Mol_Name2Idx = {";
    int i = 0;
    for (auto MolName : MolNames) {
        ofs << "'" << MolName << "' : np.array([" << i << "]), ";
        i++;
    } ofs << "}" << endl;
    ofs << in+ in+ "self.Legends = ['SimStep', 'Vol', " << JoinStr2Str(MolNames) << "]" << endl;
    ofs << endl;

    // Initialize all molecule counts
    std::vector<int> Idx_Mol = Context.GetIdxListFromMoleculeList("Molecule");
    std::vector<float> InitialCount_Molecules;
    std::vector<float> MolarityFactor_Molecules;

    for (auto& molecule : Context.MoleculeList) {
        float Count = Context.GetInitialCountByName_CountList(molecule->Name);
        float MolarityFactor = Context.GetMolarityFactorByName_CountList(molecule->Name);
        InitialCount_Molecules.push_back(Count);
        MolarityFactor_Molecules.push_back(MolarityFactor);
    }

    ofs << in+ in+ "Idx_Mol = np.array([";
    if (ObjLoc.empty()) {
        ofs << JoinInt2Str_Idx(Idx_Mol);
    } else {
        for (auto objLoc : ObjLoc) {
            ofs << "[" << JoinInt2Str_Idx(Idx_Mol) << "], ";
        }
    }
    ofs << "], ndmin=2)" << endl;

    ofs << in+ in+ "Count_Mol = np.array([";
    if (ObjLoc.empty()) {
        ofs << JoinFloat2Str(InitialCount_Molecules);
    } else {
        for (auto objLoc : ObjLoc) {
            ofs << "[" << JoinFloat2Str(InitialCount_Molecules) << "], ";
        }
    }
    ofs << "])" << endl;

    ofs << in+ in+ "MolarityFactor_Mol = np.array([" << JoinFloat2Str(MolarityFactor_Molecules) << "])" << endl;
    ofs << in+ in+ "MolarityFactor_Mol = np.where(MolarityFactor_Mol == 1, self.Vol, 1)" << endl;
    ofs << in+ in+ "Count_Mol *= MolarityFactor_Mol" << endl;
    ofs << in+ in+ "np.put_along_axis(self.Count_All, Idx_Mol, Count_Mol, axis=1)" << endl;
    ofs << endl;

    ofs << in+ "def GetMolNames(self):" << endl;
    ofs << in+ in+ "return self.Mol_Names" << endl;
    ofs << endl;

    ofs << in+ "def GetDistNames(self):" << endl;
    ofs << in+ in+ "return self.Dist_Names" << endl;
    ofs << endl;

    ofs << in+ "def GetPosNames(self):" << endl;
    ofs << in+ in+ "return self.Pos_Names" << endl;
    ofs << endl;

    ofs << in+ "def ExportLegend(self):" << endl;
    ofs << in+ in+ "return self.Legends" << endl;
    ofs << endl;

    ofs << in+ "def ExportData(self, Time):" << endl;
    ofs << in+ in+ "Data = np.asmatrix(np.zeros(2 + " << MolNames.size() << "))" << endl;
    i = 0;
    int i_SimStep = i + 1;
    ofs << in+ in+ "Data[:, " << i << ":" << i_SimStep << "] = Time" << endl;

    int i_Vol = i_SimStep + 1;
    ofs << in+ in+ "Data[:, " << i_SimStep << ":" << i_Vol << "] = self.Vol" << endl;

    int i_Count_Mol = i_Vol + MolNames.size();
    ofs << in+ in+ "Data[:, " << i_Vol << ":" << i_Count_Mol << "] = self.Count_All" << endl;

    ofs << in+ in+ "return Data" << endl;
    ofs << endl;

    // class FDataset
    ofs << "class FDataset:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.Data = None" << endl;
    ofs << endl;

    ofs << in+ "def PrintLegend(self):" << endl;
    ofs << in+ in+ "print(self.Legend)" << endl;
    ofs << endl;

    ofs << in+ "def PrintData(self):" << endl;
    ofs << in+ in+ "print(self.Data)" << endl;
    ofs << endl;

    ofs << in+ "def ExportDataAsList(self):" << endl;
    ofs << in+ in+ "return np.array(self.Data).flatten().tolist()" << endl;
    ofs << endl;

    // class FSimulation
    ofs << "class FSimulation:" << endl;
    ofs << in+ "def __init__(self, InState, InDataset, InDataManager):" << endl;
    ofs << in+ in+ "self.N_SimSteps = 0" << endl;
    ofs << in+ in+ "self.SimStep = 0" << endl;
    ofs << in+ in+ "self.SimTimeResolutionPerSecond = 0" << endl;
    ofs << endl;
    ofs << in+ in+ "self.State = InState" << endl;
    ofs << in+ in+ "self.Dataset = InDataset" << endl;
    ofs << in+ in+ "self.DataManager = InDataManager" << endl;
    ofs << endl;


    // Distribution to Count
    // Update spatial distribution to a specific coord values --> to LocationList method

    if (!MolLoc.empty()) {
        for (auto molLoc : MolLoc) {
            ofs << in+ in+ "self.Idx_DistToCoord_" << molLoc->Name << " = None" << endl;
        }
    }

    // Restore --> to CountList method
    std::vector<std::string> Names;
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if (count->End == -1) { // indicator for fixed
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                ofs << in+ in+ "self.Idx_Restore_" << Name << " = None" << endl;
                Names.push_back(Name);
            }
        }
    }

    // Event --> to CountList method
    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if ((count->End >= 0) & (count->Begin != count->End)) {
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                ofs << in+ in+ "self.Idx_Event_" << Name << " = None" << endl;
                Names.push_back(Name);
            }
        }
    }

    // Homeostasis --> to SimModule cmd option
    ofs << in+ in+ "self.Count_Prev = 0" << endl;
    ofs << in+ in+ "self.Idx_Count_Homeostasis = None" << endl;
    ofs << in+ in+ "self.Idx_Pos_Homeostasis = None" << endl;
    ofs << endl;

//    old
//    // Restore
//    for (auto& molecule : Context.MoleculeList) {
//        auto& count = molecule->Count;
//        if (count.Fixed) {
//            ofs << in+ in+ "self.Idx_Restore_" << molecule->Name << " = None" << endl;
//        }
//    }
//
//    // Event
//    for (auto& count : Context.CountList) {
//        if (Begin != 0) {
//            }
//        }
//    }

//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                ofs << in+ in+ "self.Idx_Restore_" << Pathway.Name << " = None" << endl;
//            }
//        }
//    }

    ofs << in+ in+ "# Debugging" << endl;
    ofs << in+ in+ "self.Debug_Idx_Molecules = list()" << endl;
    ofs << in+ in+ "self.UnitTxt = ''" << endl;
    ofs << in+ in+ "self.Unit = 0" << endl;
    ofs << endl;

    ofs << in+ in+ "# Debugging_Spatial" << endl;
    ofs << in+ in+ "self.Debug_Idx_Pos = list()" << endl;
    ofs << in+ in+ "self.Debug_Idx_Dist = list()" << endl;
    ofs << endl;

    ofs << in+ "def Initialize(self, InN_SimSteps=1000, InTimeResolution=100):" << endl;
    ofs << in+ in+ "print('Simulation Initialized...')" << endl;
    ofs << in+ in+ "self.N_SimSteps = np.array([InN_SimSteps])" << endl;
    ofs << in+ in+ "self.SimTimeResolutionPerSecond = InTimeResolution" << endl;
    ofs << endl;

    // Distribution To Count
    i = 0;
    if (!MolLoc.empty()) {
        for (auto molLoc : MolLoc) {
            int Idx = Context.GetIdxByName_MoleculeList(molLoc->Name);
            ofs << in+ in+ "self.Idx_DistToCoord_" << molLoc->Name << " = np.array([" << std::to_string(Idx) << "])" << endl;
            i++;
        }
    }

    // Restore --> to CountList method
    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if (count->End == -1) { // indicator for fixed
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                int Idx = Context.GetIdxByName_MoleculeList(Name);
                ofs << in+ in+ "self.Idx_Restore_" << Name << " = np.array([" << std::to_string(Idx) << "])" << endl;
                Names.push_back(Name);
            }
        }
    }

    // Event --> to CountList method
    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;
        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }
        if ((count->End >= 0) & (count->Begin != count->End)) { //
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                int Idx = Context.GetIdxByName_MoleculeList(Name);
                ofs << in+ in+ "self.Idx_Event_" << Name << " = np.array([" << std::to_string(Idx) << "], ndmin=2)" << endl;
                Names.push_back(Name);
            }
        }
    }
    ofs << endl;

//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                std::string MoleculeToRestore;
//                int Idx;
//
//                MoleculeToRestore = "acetyl-CoA";
//                Idx = Context.GetIdxByName_MoleculeList(MoleculeToRestore);
//                ofs << in+ in+ "self.Idx_Restore_" << Pathway.Name << " = np.asmatrix([" << Idx << "]) # " << MoleculeToRestore << endl;
//
//	                MoleculeToRestore.clear();
//                MoleculeToRestore = "malate";
//                Idx = Context.GetIdxByName_MoleculeList(MoleculeToRestore);
//                ofs << in+ in+ "self.Idx_Restore_" << Pathway.Name << " = np.asmatrix([" << Idx << "]) # " << MoleculeToRestore << endl;
//            }
//        }
//    }
    ofs << in+ in+ "self.State.Initialize()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Legend Export" << endl;
    ofs << in+ in+ "self.Dataset.Legend = self.State.ExportLegend()" << endl;
    ofs << in+ in+ "self.DataManager.SetLegend(self.Dataset.Legend)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Data Export" << endl;
    ofs << in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Debugging" << endl;
    ofs << in+ in+ "self.Debug_SetIdxMoleculesToTrack()" << endl;
    ofs << in+ in+ "self.Debug_SetIdxDistAndPosToTrack()" << endl;
    ofs << in+ in+ "self.Debug_SetUnit(Unit)" << endl;
    ofs << endl;

    // homeostasis--threshold index setting
    ofs << in+ "def SetIdxForHomeostasis(self, MoleculeNameList):" << endl;
    ofs << in+ in+ "if MoleculeNameList:" << endl;
    ofs << in+ in+ in+ "self.Idx_Count_Homeostasis = np.array([self.GetMolIdx(Name)[0] for Name in MoleculeNameList])" << endl;
    ofs << in+ in+ in+ "self.Idx_Pos_Homeostasis = np.array(range(len(MoleculeNameList)))" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "self.Idx_Count_Homeostasis = np.array(range(self.State.Count_All.shape[1]))" << endl;
    ofs << in+ in+ in+ "self.Idx_Pos_Homeostasis = np.array(range(self.State.Count_All.shape[1]))" << endl;
    ofs << endl;

    // regular simloop
    ofs << in+ "def SimLoop_WithSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#"; }
//    ofs << "self.Debug_PrintCounts(DisplayCount)" << endl;
//    ofs << in+ in;
////    if (!Option.bDebug) { ofs << "#"; }
//    ofs << "self.Debug_PrintDistributions()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Run Spatial Simulation" << endl;
    ofs << in+ in+ "self.SpatialSimulation()" << endl; // TODO: Take delta, instead of updating directly
    ofs << endl;

    ofs << in+ in+ "# Run Reactions" << endl;
    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;


    ofs << in+ in+ "# Update Substrate Count" << endl;
    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << endl;

    ofs << in+ in+ "# temporary: Update Threshold" << endl;
    ofs << in+ in+ "self.UpdateThreshold()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Update Spatially Distributed Molecules On Count" << endl;
    ofs << in+ in+ "self.DistributionToCount()" << endl;
    ofs << in+ in+ "self.CountToDistribution()" << endl;
    ofs << endl;

    ofs << in+ in+ "# Restore Substrate Count for Sustained Substrate InTransporter" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ "def SimLoop_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << "self.Debug_PrintCounts(DisplayCount)" << endl;
//    ofs << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << endl;

    // simloop for receptivity
    ofs << in+ "def SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution(self):" << endl;
    ofs << in+ in+ "self.IncrementSimStep()" << endl;

//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
//    ofs << in+ in;
//    if (!Option.bDebug) { ofs << "#";}
//    ofs << "self.Debug_PrintCounts(DisplayCount)" << endl;
//    ofs << endl;

    ofs << in+ in+ "self.NonSpatialSimulation()" << endl;

    ofs << in+ in+ "self.UpdateCounts()" << endl;
    ofs << in+ in+ "self.RestoreMoleculeCount()" << endl;
    ofs << in+ in+ "self.DistributionToCount()" << endl;
    ofs << in+ in+ "self.CountToDistribution()" << endl;
    ofs << endl;

    ofs << in+ "def Run(self, Spatial=0):" << endl;
    ofs << in+ in+ "if Spatial:" << endl;
    ofs << in+ in+ in+ "self.Run_WithSpatialSimulation()" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "self.Run_WithoutSpatialSimulation()" << endl;
    ofs << endl;

    ofs << in+ "def Run_WithSpatialSimulation(self):" << endl;
    ofs << in+ in+ "print('Simulation Run Begins...')" << endl;
    ofs << endl;

    ofs << in+ in+ "while self.SimStep < self.N_SimSteps:" << endl;
    ofs << in+ in+ in+ "self.SimLoop_WithSpatialSimulation()" << endl;

    ofs << in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ "def Run_WithoutSpatialSimulation(self):" << endl;
    ofs << in+ in+ "print('Simulation Run_WithoutSpatialSimulation Begins...')" << endl;
    ofs << endl;

    ofs << in+ in+ "while self.SimStep < self.N_SimSteps:" << endl;
    ofs << in+ in+ in+ "self.SimLoop_WithoutSpatialSimulation()" << endl;

    ofs << in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;

    ofs << in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    ofs << in+ in+ "print('Simulation Run_WithoutSpatialSimulation Completed')" << endl;
    ofs << endl;

    // TODO: scale up exporting
    ofs << in+ "def ExportData(self):" << endl;
    if (Context.LocationList.empty()) {
        ofs << in+ in+ "self.Dataset.Data = self.State.ExportData(self.SimStep/self.SimTimeResolutionPerSecond)" << endl;
        ofs << in+ in+ "self.DataManager.Add(self.Dataset.Data)" << endl;
    } else { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def SaveState(self, FileName):" << endl;
    ofs << in+ in+ "self.DataManager.SaveToNpyFile(self.GetCount_All(), FileName)" << endl;
    ofs << endl;

    ofs << in+ "def ApplySimTimeResolution(self, Rate):" << endl;
    ofs << in+ in+ "return Rate / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    // Distribution to Count

    // get coord for 'e'
    // get 'L' count from distribution by 'e' coord
    // get idx for 'L'
    // Update count at idx


    // TODO: make it conditional for specially treated molecules in the system
    ofs << in+ "def DistributionToCount(self):" << endl;
    bool PassSwitch = true;
    if (!MolLoc.empty() & !ObjLoc.empty()) {
        for (auto molLoc : MolLoc) {
            std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
            for (auto UniqueName : ObjUniqueNames) {
                if (molLoc->Name == "L") {
                    ofs << in + in + "Count = self.GetCountFromDistributionByNameAndPos('" << molLoc->Name << "', " << "'" << UniqueName << "')" << endl;
                    ofs << in + in + "self.State.Count_All[:, self.Idx_DistToCoord_" << molLoc->Name
                        << "] = Count.reshape(-1, 1)" << endl;
                    PassSwitch = false;
                }
            }
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "def CountToDistribution(self):" << endl;
    {
//    if (!MolLoc.empty() & !ObjLoc.empty()) {
//        for (auto molLoc : MolLoc) {
//            std::vector<std::string> ObjUniqueNames = Context.GetUniqueNames_LocationList("Compartment");
//            for (auto UniqueName : ObjUniqueNames) {
//                ofs << in + in + "Count = self.GetCountFromDistributionByNameAndPos('" << molLoc->Name << "', " << "'" << UniqueName << "')" << endl;
//                ofs << in + in + "self.State.Count_All[:, self.Idx_DistToCoord_" << molLoc->Name
//                    << "] = Count.transpose()" << endl;
//            }
//        }
//    } else {
        ofs << in+ in+ "pass" << endl;
    }
    ofs << endl;

    // Restore
    ofs << in+ "def RestoreMoleculeCount(self):" << endl;

    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;

        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }

        if (count->End == -1) { // indicator for fixed
            if (std::find(Names.begin(), Names.end(), Name) == Names.end()) {
                std::string Amount;
                float MolarityFactor = Context.GetMolarityFactorByName_CountList(Name);

                if (MolarityFactor)      { Amount = Utils::SciFloat2Str(count->Amount) + " * self.State.Vol"; }
                else                     { Amount = Utils::SciFloat2Str(count->Amount); }

                ofs << in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << Name << ".reshape(1, -1), " << Amount << ", axis=1)" << endl;
                Names.push_back(Name);
            }
        }
    }

    if (Names.empty()) {
        ofs << in+ in+ "pass" << endl;
    }


//    if (!Context.PathwayList.empty()){
//        for (auto& Pathway : Context.PathwayList) {
//            if (Pathway.Name == "TCA") {
//                std::string MoleculeToRestore = "acetyl-CoA";
//                int Count = 446331;
//                ofs << in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Restore_" << Pathway.Name << ", " << std::to_string(Count) << ", axis=1)   # " << MoleculeToRestore << endl;
//                Pass = false;
//            }
//        }
//    }
    ofs << endl;

    // Event
    ofs << in+ "def TriggerEventMoleculeCount(self):" << endl;
    ofs << in+ in+ "Time = self.SimStep / self.SimTimeResolutionPerSecond" << endl;
    bool ElseSwitch = false;

    Names.clear();
    for (auto& count : Context.CountList) {
        std::string Name = count->Name;

        // ignore if not relevant to the system
        if (std::find(MolNames.begin(), MolNames.end(), Name) == MolNames.end()) {
            continue;
        }

        std::string Amount;
        float MolarityFactor = Context.GetMolarityFactorByName_CountList(Name);

        if (MolarityFactor) { Amount = Utils::SciFloat2Str(count->Amount) + " * self.State.Vol"; }
        else                { Amount = Utils::SciFloat2Str(count->Amount); }


        if ((count->End >= 0) & (count->Begin != count->End)) { //
            ofs << in+ in+ "if (Time >= " << Utils::SciFloat2Str(count->Begin) << ") & (Time < " << Utils::SciFloat2Str(count->End) << "):" << endl;
            ofs << in+ in+ in+ "np.put_along_axis(self.State.Count_All, self.Idx_Event_" << Name << ", " << Amount << ", axis=1)" << endl;
            Names.push_back(Name);
        }
    }
    ofs << endl;

    // temporary code
    ofs << in + "def GetNewThresholdValues(self):" << endl;
    float ThresholdFactor = 0.999;
    ofs << in+ in+ "return self.GetCount_All()[:, self.Idx_Count_Homeostasis].transpose() * " << Utils::SciFloat2Str(ThresholdFactor) << endl;
    ofs << endl;

    ofs << in + "def SetThreshold(self):" << endl;
    ofs << in+ in+ "self.State.Pos_Threshold = self.GetNewThresholdValues()" << endl;
    ofs << endl;

    ofs << in + "def UpdateThreshold(self):" << endl;
    float HomeostasisFactor = 1e-7;
    ofs << in+ in+ "Count_Now = self.GetCount_All()[:, self.Idx_Count_Homeostasis].transpose()" << endl;
    ofs << in+ in+ "bSteadyState = abs(Count_Now - self.Count_Prev) / Count_Now < " << Utils::SciFloat2Str(HomeostasisFactor) << endl;
    ofs << in+ in+ "Idx_SteadyState = np.where(bSteadyState)" << endl;
    ofs << in+ in+ "self.State.Pos_Threshold[Idx_SteadyState] = Count_Now[Idx_SteadyState]" << endl;
//    ofs << in+ in+ "if np.any(Idx_SteadyState):" << endl;
//    ofs << in+ in+ in+ "print('[Homeostasis] Steady state has been updated.')" << endl;
    ofs << in+ in+ "self.Count_Prev = Count_Now" << endl;
    ofs << endl;

    ofs << in + "def Homeostasis(self, MoleculeNames=None):" << endl;
    // TODO: Get homeostasis info from FPathway,
    ofs << in+ in+ "bNotHomeostasis = True" << endl;
    ofs << in+ in+ "self.SetIdxForHomeostasis(MoleculeNames)" << endl;
    ofs << in+ in+ "MoleculeNames_Str = ListToStr(MoleculeNames, '_')" << endl;

    ofs << in+ in+ "# Temporary homeostasis save file" << endl;
    ofs << in+ in+ "FileName_Homeostasis = 'SimSave_Homeostasis_%s" << int(Numbers::RandomNumber(0, 1000)) << "' % MoleculeNames_Str" << endl;
    ofs << endl;
    ofs << in+ in+ "if os.path.exists('%s.npy' % FileName_Homeostasis):" << endl;
    ofs << in+ in+ in+ "print('[Homeostasis] Loading previous steady state...')" << endl;
    ofs << in+ in+ in+ "self.State.Count_All = np.load('%s.npy' % FileName_Homeostasis)" << endl;
    ofs << in+ in+ in+ "print('[Homeostasis] Steady state has been achieved!')" << endl;
    ofs << in+ in+ in+ "self.SpecialRedistribution()" << endl;
    ofs << in+ in+ in+ "if MoleculeNames:" << endl;
    ofs << in+ in+ in+ in+ "self.SetThreshold()" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "print('[Homeostasis] Running simulation to achieve steady state...')" << endl;

    ofs << in+ in+ in+ "while (bNotHomeostasis):" << endl;
    ofs << in+ in+ in+ in+ "self.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()" << endl;
    ofs << in+ in+ in+ in+ "Count_Now = self.GetCount_All()[:, self.Idx_Count_Homeostasis].transpose()" << endl;
    ofs << in+ in+ in+ in+ "if np.all(Count_Now > 0) and np.all(abs(Count_Now - self.Count_Prev) / Count_Now < " << Utils::SciFloat2Str(HomeostasisFactor) << "):" << endl;
    ofs << in+ in+ in+ in+ in+ "bNotHomeostasis = False" << endl;
    ofs << in+ in+ in+ in+ in+ "print('[Homeostasis] Steady state has been achieved.')" << endl;
    ofs << in+ in+ in+ in+ in+ "self.SaveState(FileName_Homeostasis)" << endl;
    ofs << in+ in+ in+ in+ in+ "print('[Homeostasis] Steady state has been saved: %s' % FileName_Homeostasis)" << endl;
    ofs << in+ in+ in+ in+ in+ "if MoleculeNames:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "self.SetThreshold()" << endl;
    ofs << in+ in+ in+ in+ in+ "self.SpecialRedistribution()" << endl;
//    ofs << in+ in+ in+ in+ "print('[Homeostasis] achieved  : " << name << " @', self.Debug_ApplyUnit(Current), self.UnitTxt)" << endl;
//    ofs << in+ in+ in+ in+ "print('[Homeostasis] threshold : " << name << " @', self.Debug_ApplyUnit(Threshold), self.UnitTxt)" << endl;
    //ofs << in+ in+ in+ in+ "print('Homeostasis achieved for      : " << name << " @ {:.010f}'.format(self.Debug_ApplyUnit(" << Now << ")), self.UnitTxt)" << endl;
    //ofs << in+ in+ in+ in+ "print('Homeostasis threshold set for : " << name << " @ {:.010f}'.format(self.Debug_ApplyUnit(" << Now << " * " << Utils::SciFloat2Str(ThresholdFactor) << ")), self.UnitTxt)" << endl;
    ofs << in+ in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ in+ "self.Count_Prev = Count_Now" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ in+ "# Trigger Event on Substrate Count" << endl;
    ofs << in+ in+ in+ in+ in+ "self.TriggerEventMoleculeCount()" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ in+ "# Save and Export Data" << endl;
    ofs << in+ in+ in+ in+ in+ "self.ExportData()" << endl;
    ofs << endl;

    // Special Distribution Reset for glucose 'L' distribution
    ofs << in+ "def SpecialRedistribution(self):" << endl;
    PassSwitch = true;
//    for (int i = 0; i < MolLoc.size(); i++) {
//        auto Coord = MolLoc[i]->Coord;
//        auto Amount = Context.GetInitialCountByName_CountList(MolLoc[i]->Name);
//        std::string Shape, Pattern;
//        int Size;
//
//        if (MolLoc[i]->Name == "L") {
//            Shape = "'circle'"; Size = 500; Pattern = "'diffuse'";
//
//            ofs << in + in + "self.State.Dist_All[" << i << "] = SimF.InitializeDistribution(self.State.Dimension_X, self.State.Dimension_Y, "
//                << MolLoc[i]->Coord[0] << ", " << MolLoc[i]->Coord[1] << ", " << Amount << ", " << 0 ;
//            if (!Shape.empty()) {
//                ofs << ", " << Shape << ", " << Size << ", " << Pattern;
//                }
//            ofs << ")" << endl;
//            PassSwitch = false;
//        }
//    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;


    if (!Context.LocationList.empty()) {
        ofs << in+ "# Spatial Simulation related routines" << endl;
        ofs << in + "def SpatialSimulation(self):" << endl;
        ofs << in + in + "self.SpatialDiffusion()" << endl;
        if (!TransporterReactionTypes.empty()) {
            ofs << in + in + "self.TransporterReactions()" << endl;
        }

            ofs << endl;
            ofs << in+ in+ "self.Debug_PrintSimStepTime()" << endl;
            ofs << in+ in+ "print()" << endl;
            ofs << in+ in+ "self.Debug_PrintCountsAndDistributions() # Temporary placement. Update with implementing delta for spatial simulation" << endl;
            ofs << endl;

        ofs << in + in + "self.SpatialLocation()" << endl;
        ofs << endl;
    }

    ofs << in+ "def SpatialDiffusion(self):" << endl;
    PassSwitch = true;
    if (!MolLoc.empty()) {
        int N_Dist = Context.GetNames_LocationList("Molecule").size();
        // TODO: update to 3d array
        for (i = 0; i < N_Dist; i++) {
            if (MolLoc[i]->Name == "L") { continue; }
            else {
                ofs << in+ in+ "self.State.Dist_All[" << i << "] = SimF.DiffuseDistribution_4Cell(self.State.Dist_All[" << i << "])" << endl;
                PassSwitch = false;
            }
            // Future implementation
            // ofs << in+ in+ "self.State.Dist_All[" << i << "] = SimF.DiffuseDistribution_8Cell(self.State.Dist_All[" << i << "])" << endl;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    // Print StandardReaction for each Reaction Type
    ofs << in+ "def TransporterReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : TransporterReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.TransporterReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    // TransporterReaction function
    for (auto& Type : TransporterReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            bool bMolaritySys = Context.CheckMolarityFactorTrueForAny_CountList();

            std::string AmountTextStr;
            if (bMolaritySys) { AmountTextStr = "Conc_"; }
            else              { AmountTextStr = "Count_"; }

            ofs << in+ "def TransporterReaction_" << Type << "(self):" << endl;

//            if (Option.bDebug) {
                ofs << endl;
                ofs << in+ in+ "# Debugging tool" << endl;
                ofs << in+ in+ "# self.Debug_PrintNames(PUT_IDX_HERE)" << endl;
                ofs << endl;
//            }

            std::string Idx_Cargo = "self.State.Idx_Cargo_" + Type;
            std::string Idx_Dist_Cargo = "self.State.Idx_Dist_Cargo_" + Type;
            std::string Idx_Transporter = "self.State.Idx_Transporter_" + Type;

            std::string Amount_Cargo_Inside = AmountTextStr + "Cargo_Inside";
            std::string Amount_Cargo_Outside = AmountTextStr + "Cargo_Outside";
            std::string Amount_Transporter = AmountTextStr + "Transporter";

            ofs << in+ in+ Amount_Transporter << " = ";
            std::string Line = "self.State.Count_All[:, " + Idx_Transporter + "]";
            if (bMolaritySys) { Line = "SimF.CountToConc(" + Line + ", self.State.Vol)"; }
            ofs << Line << endl;

            // Get Concentrations
            ofs << in+ in+ Amount_Cargo_Inside << " = ";
            Line = "self.State.Count_All[:, " + Idx_Cargo + "]";
            if (bMolaritySys) { Line = "SimF.CountToConc(" + Line + ", self.State.Vol)"; }
            ofs << Line << endl;

            ofs << in+ in+ Amount_Cargo_Outside << " = ";
            Line = "self.GetCountFromDistribution(" + Idx_Dist_Cargo + ", self.State.Pos_X, self.State.Pos_Y).transpose()";
            if (bMolaritySys) { Line = "SimF.CountToConc(" + Line + ", self.State.Vol)"; }
            ofs << Line << endl;

//            // Regulators
//            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
//                ofs << in+ in+ AmountTextStr << "Regulator = ";
//                if (bMolaritySys) {
//                    ofs << "SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
//                } else {
//                    ofs << "self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "]" << endl;
//                }
//            }

            // Calculate Rate
            ofs << in+ in+ "Rate = SimF.Eqn_" << Type << "(" <<
            Amount_Transporter << ", " <<
            Amount_Cargo_Inside << ", " <<
            Amount_Cargo_Outside << ", " <<
            "self.State.Const_ki_" << Type << ", " <<
            "self.State.Const_ko_" << Type << ")" << endl;

            //            // Regulators
//            if (Typing[substrate] != StandardReactionTypes[0]) {
//                ofs << AmountTextStr << "Regulator, ";
//            }
//            // Reaction constants
//            ofs << "self.State.Const_k_" << substrate << "_" << Type;
//            if (Typing[substrate] != StandardReactionTypes[0]) {
//                ofs << ", self.State.Const_K_" << Type << ", ";
//                ofs << "self.State.Const_n_" << Type;
//            }
//            ofs << ")" << endl;

            // Apply Time Resolution
            ofs << in+ in+ "Rate = self.ApplySimTimeResolution(Rate)" << endl;

            if (bMolaritySys) {
                // Apply stoichiometry
                ofs << in+ in+ "dConc_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
                // Convert to counts
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            } else {
                // Apply stoichiometry
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            }
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << in+ in+ "self.AddToDist(" << Idx_Dist_Cargo << ", self.State.Pos_X, self.State.Pos_Y, -dCount_Mol_InStoichMatrix.transpose())" << endl;
            ofs << endl;
        }
    }

    ofs << in+ "def SpatialLocation(self):" << endl;
    ofs << in+ in+ "HomeostasisMolecule = self.GetCount(self.Idx_Count_Homeostasis).transpose()[0]" << endl; // TODO:set up dynamic indexing
    ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = SimF.BacterialChemotaxis(np.array(HomeostasisMolecule), self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle, self.State.Pos_Threshold)" << endl;
    ofs << in+ in+ "self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle = SimF.CorrectOutOfBounds(self.State.Pos_X, self.State.Pos_Y, self.State.Pos_Angle, self.State.Dimension_X, self.State.Dimension_Y)" << endl;
    ofs << endl;

    ofs << in+ "def NonSpatialSimulation(self):" << endl;
    PassSwitch = true;
    if (!StandardReactionTypes.empty()) {
        for (auto &Type: StandardReactionTypes) {
            std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
            if (!ReactionSubList.empty()) {
                ofs << in + in + "self.StandardReactions()" << endl;
                ofs << endl;
                PassSwitch = false;
                break;
            }
        }
    }

    if (!EnzReactionTypes.empty()) {
        for (auto &Type: EnzReactionTypes) {
            std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
            if (!ReactionSubList.empty()) {
                ofs << in + in + "self.EnzymaticReactions()" << endl;
                ofs << endl;
                PassSwitch = false;
                break;
            }
        }
    }

    // TODO: encapsulate this part with each polymerase to allow more process-specific customization
    if (!PolymeraseList.empty()){
        ofs << in+ in+ "self.Polymerase_InitiationReactions()" << endl;
        ofs << in+ in+ "self.Polymerase_ElongationReactions()" << endl;
        ofs << in+ in+ "self.Polymerase_TerminationReactions()" << endl;
        PassSwitch = false;
    }

    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    ofs << in+ "# Biochemical Reaction related routines" << endl;
    // StandardReaction function
    for (auto& Type : StandardReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            // Typing set up
            std::map<std::string, std::string> Typing;
            std::vector<std::string> SubstrateTypes { "Reactant", "Product" };
            Typing.emplace("Reactant", Type);
            Typing.emplace("Product", StandardReactionTypes[0]); // default type

            bool bMolaritySys = Context.CheckMolarityFactorTrueForAny_CountList();
            std::string AmountTextStr;
            if (bMolaritySys) { AmountTextStr = "Conc_"; }
            else              { AmountTextStr = "Count_"; }

            ofs << in+ "def StandardReaction_" << Type << "(self):" << endl;

//            if (Option.bDebug) {
                ofs << endl;
                ofs << in+ in+ "# Debugging tool" << endl;
                ofs << in+ in+ "# self.Debug_PrintNames(PUT_IDX_HERE)" << endl;
                ofs << endl;
//            }

            // Get Concentrations
            // Reactants and products
            for (auto& substrate : SubstrateTypes) {
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << in+ in+ AmountTextStr << substrate << "_" << i << " = ";
                    if (bMolaritySys) {
                        ofs << "SimF.CountToConc(self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "], self.State.Vol)" << endl;
                    } else {
                        ofs << "self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "]" << endl;
                    }
                }
            }

            // Regulators
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << in+ in+ AmountTextStr << "Regulator = ";
                if (bMolaritySys) {
                    ofs << "SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
                } else {
                    ofs << "self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "]" << endl;
                }
            }

            // Calculate Rate
            for (auto& substrate : SubstrateTypes) {
                ofs << in+ in+ "Rate_" << substrate << " = SimF.Eqn_" << Typing[substrate] << "_" << N_MoleculesAllowed << "(";

                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << AmountTextStr << substrate << "_" << i << ", ";
                }

                // Regulators
                if (Typing[substrate] != StandardReactionTypes[0]) {
                    ofs << AmountTextStr << "Regulator, ";
                }
                // Reaction constants
                ofs << "self.State.Const_k_" << substrate << "_" << Type;
                if (Typing[substrate] != StandardReactionTypes[0]) {
                    ofs << ", self.State.Const_K_" << Type << ", ";
                    ofs << "self.State.Const_n_" << Type;
                }
                ofs << ")" << endl;

                // Apply Time Resolution
                ofs << in+ in+ "Rate_" << substrate << " = self.ApplySimTimeResolution(Rate_" << substrate << ")" << endl;

                // compare conc
                ofs << in+ in+ "Rate_" << substrate << " = SimF.CheckRateAndConc_" << N_MoleculesAllowed << "(Rate_" << substrate;
                for (int i = 0; i < N_MoleculesAllowed; i++) {
                    ofs << ", " << AmountTextStr << substrate << "_" << i;
                }
                ofs << ")" << endl;
            }

            // Tally Rates
            ofs << in+ in+ "Rate = Rate_Reactant - Rate_Product" << endl;

            if (bMolaritySys) {
                // Apply stoichiometry
                ofs << in+ in+ "dConc_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
                // Convert to counts
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            } else {
                // Apply stoichiometry
                ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            }
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    // EnzReaction function // TODO: Add bMolaritySys option as above
    for (auto& Type : EnzReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {

            // Typing set up
            std::map<std::string, std::string> Typing;
            std::vector<std::string> SubstrateTypes { "Reactant", "Product" };
            Typing.emplace("Reactant", Type);
            Typing.emplace("Product", EnzReactionTypes[0]); // default type

            ofs << in+ "def EnzymaticReaction_" << Type << "(self):" << endl;

//            if (Option.bDebug) {
                ofs << endl;
                ofs << in+ in+ "# Debugging tool" << endl;
                ofs << in+ in+ "# self.Debug_PrintNames(PUT_IDX_HERE)" << endl;
                ofs << endl;
//            }

            // Get Concentrations
            // Enzyme
            ofs << in+ in+ "Conc_Enz = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Enz_" << Type << "], self.State.Vol)" << endl;
            // Reactants and products or EnzSubstrate
            if (Type.find("Enz_Standard") != std::string::npos) {
                for (auto& substrate : SubstrateTypes) {
                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << in+ in+ "Conc_" << substrate << "_" << i << " = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_" << substrate << "_" << i << "_" << Type << "], self.State.Vol)" << endl;
                    }
                }
            } else if (Type.find("Enz_MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ "Conc_EnzSub = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_EnzSub_" << Type << "], self.State.Vol)" << endl;
            } else {
                Utils::Assertion(false, "Unsupported Enz Reaction Type: " + Type);
            }

            // Regulators
            if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                ofs << in+ in+ "Conc_Regulator = SimF.CountToConc(self.State.Count_All[:, self.State.Idx_Regulator_" << Type << "], self.State.Vol)" << endl;
            }

            // Calculate Rate
            if (Type.find("Enz_Standard") != std::string::npos) {
                for (auto& substrate : SubstrateTypes) {
                    ofs << in+ in+ "Rate_" << substrate << " = SimF.Eqn_" << Typing[substrate] << "_" << N_MoleculesAllowed << "(Conc_Enz, ";

                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << "Conc_" << substrate << "_" << i << ", ";
                    }

                    // Regulators
                    if (Typing[substrate] != EnzReactionTypes[0]) {
                        ofs << "Conc_Regulator, ";
                    }
                    // Reaction constants
                    ofs << "self.State.Const_k_" << substrate << "_" << Type;
                    if (Typing[substrate] != EnzReactionTypes[0]) {
                        ofs << ", self.State.Const_K_" << Type << ", ";
                        ofs << "self.State.Const_n_" << Type;
                    }
                    ofs << ")" << endl;

                    // Apply Time Resolution
                    ofs << in+ in+ "Rate_" << substrate << " = self.ApplySimTimeResolution(Rate_" << substrate << ")" << endl;

                    // compare conc
                    ofs << in+ in+ "Rate_" << substrate << " = SimF.CheckRateAndConc_" << N_MoleculesAllowed << "(Rate_" << substrate;
                    for (int i = 0; i < N_MoleculesAllowed; i++) {
                        ofs << ", Conc_" << substrate << "_" << i;
                    }
                    ofs << ")" << endl;
                }

                // Tally Rates
                ofs << in+ in+ "Rate = Rate_Reactant - Rate_Product" << endl;

            } else if (Type.find("Enz_MichaelisMenten") != std::string::npos) {
                ofs << in+ in+ "Rate = SimF.Eqn_" << Type << "(Conc_Enz, Conc_EnzSub";

                // Regulators
                if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                    ofs << ", Conc_Regulator";
                }

                // Reaction constants
                ofs << ", self.State.Const_kcat_" << Type << ", self.State.Const_KM_" << Type;

                if ((Type.find("Inhibition") != std::string::npos) || (Type.find("Activation") != std::string::npos)) {
                    ofs << ", self.State.Const_K_" << Type;
                    if (Type.find("Allosteric") != std::string::npos) {
                        ofs << ", self.State.Const_n_" << Type;
                    }
                }
                ofs << ")" << endl;

                // Apply TimeResolution to the rate
                ofs << in+ in+ "Rate = self.ApplySimTimeResolution(Rate)" << endl;
                // compare conc for michaelis menten?
            }
            // Apply stoichiometry
            ofs << in+ in+ "dConc_Mol_InStoichMatrix = SimF.GetDerivativeFromStoichiometryMatrix(self.State.Const_StoichMatrix_" << Type <<", Rate)" << endl;
            // Convert to counts
            ofs << in+ in+ "dCount_Mol_InStoichMatrix = SimF.ConcToCount(dConc_Mol_InStoichMatrix, self.State.Vol)" << endl;
            // Apply delta counts for molecules in the stoichiometry matrix
            ofs << in+ in+ "self.AddTodCount(self.State.Idx_Mol_InStoichMatrix_" << Type << ", dCount_Mol_InStoichMatrix)" << endl;
            ofs << endl;
        }
    }

    // Print StandardReaction for each Reaction Type
    ofs << in+ "def StandardReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : StandardReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.StandardReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;


    // Print EnzymeReaction for each Reaction Type
    ofs << in+ "def EnzymaticReactions(self):" << endl;
    PassSwitch = true;
    for (auto& Type : EnzReactionTypes) {
        std::vector<const FReaction *> ReactionSubList = Context.GetSubList_ReactionList(Type);
        if (!ReactionSubList.empty()) {
            ofs << in+ in+ "self.EnzymaticReaction_" << Type << "()" << endl;
            PassSwitch = false;
        }
    }
    if (PassSwitch) { ofs << in+ in+ "pass" << endl; }
    ofs << endl;

    if (!PolymeraseList.empty()) {

        ofs << in+ "def Polymerase_InitiationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {

            Polymerase_InitiationReaction(ofs, Polymerase);
        }


        ofs << in+ "def Polymerase_ElongationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {

            Polymerase_ElongationReaction(ofs, Polymerase);
        }

        ofs << in+ "def Polymerase_TerminationReactions(self):" << endl;
        for (auto& Polymerase : PolymeraseList) {

            Polymerase_TerminationReaction(ofs, Polymerase);
        }
    }

//    if (EnzymeList.empty() & PolymeraseList.empty()) {
//            ofs << in+ in+ "pass" << endl;
//            ofs << endl;
//    }

    ofs << in+ "# Useful routines" << endl;

    ofs << in+ "def GetSimStep(self):" << endl;
    ofs << in+ in+ "return self.SimStep" << endl;
    ofs << endl;

    ofs << in+ "def GetSimTime(self):" << endl;
    ofs << in+ in+ "return self.SimStep / self.SimTimeResolutionPerSecond" << endl;
    ofs << endl;

    ofs << in+ "def IncrementSimStep(self):" << endl;
    ofs << in+ in+ "self.SimStep += 1" << endl;
    ofs << endl;

    ofs << in+ "def GetCount(self, Idx):" << endl;
    ofs << in+ in+ "return self.State.Count_All[:, Idx]" << endl;
    ofs << endl;

    ofs << in+ "def GetCount_All(self):" << endl;
    ofs << in+ in+ "return np.copy(self.State.Count_All)" << endl;
    ofs << endl;

    ofs << in+ "def GetConcentration(self, Idx):" << endl;
    ofs << in+ in+ "return SimF.CountToConc(self.State.Count_All[:, Idx])" << endl;
    ofs << endl;

    ofs << in+ "def GetDistribution(self, Idx):" << endl;
    ofs << in+ in+ "return self.State.Dist_All[Idx]" << endl;
    ofs << endl;

    ofs << in+ "def GetDistWidth(self):" << endl;
    ofs << in+ in+ "return self.State.Dimension_X" << endl;
    ofs << endl;

    ofs << in+ "def GetDistHeight(self):" << endl;
    ofs << in+ in+ "return self.State.Dimension_Y" << endl;
    ofs << endl;

    ofs << in+ "def AddTodCount(self, Idx, Values):" << endl;
    ofs << in+ in+ "dCountToAdd = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << in+ in+ "np.put_along_axis(dCountToAdd, Idx, Values, axis=1)" << endl;
    ofs << in+ in+ "dCount_All_New = self.State.dCount_All + dCountToAdd" << endl;
    ofs << in+ in+ "ZeroTest = dCount_All_New + self.State.Count_All" << endl;
    ofs << in+ in+ "self.State.dCount_All = np.where(ZeroTest < 0, dCount_All_New - ZeroTest, dCount_All_New)" << endl;
    ofs << endl;

    ofs << in+ "def UpdateCounts(self):" << endl;
    ofs << in+ in+ "self.State.Count_All += self.State.dCount_All" << endl;
    ofs << in+ in+ "self.CleardCounts()" << endl;
    ofs << endl;

    ofs << in+ "def CleardCounts(self):" << endl;
    ofs << in+ in+ "self.State.dCount_All = np.zeros_like(self.State.dCount_All)" << endl;
    ofs << endl;

    ofs << in+ "def AddToDist(self, Idx, X, Y, Values):" << endl;
    ofs << in+ in+ "self.State.Dist_All[Idx, X.astype(int), Y.astype(int)] += Values" << endl;
    ofs << endl;

    // TODO: Use SimIdx in the future
    ofs << in+ "# Temporary routines" << endl;

    ofs << in+ "def GetMolIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.Mol_Name2Idx[Name]" << endl;
    ofs << endl;

    ofs << in+ "def GetCountByName(self, Name):" << endl;
    ofs << in+ in+ "return self.GetCount(self.GetMolIdx(Name))" << endl;
    ofs << endl;

    ofs << in+ "def GetConcentrationByName(self, Name):" << endl;
    ofs << in+ in+ "return self.GetConcentration(self.GetMolIdx(Name), self.State.Vol)" << endl;
    ofs << endl;

    ofs << in+ "def GetPosIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.Pos_Name2Idx[Name]" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXY(self, Idx):" << endl;
    ofs << in+ in+ "return np.take(self.State.Pos_X, Idx), np.take(self.State.Pos_Y, Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXYByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetPosIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetPositionXY(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXYAngle(self, Idx):" << endl;
    ofs << in+ in+ "return np.take(self.State.Pos_X, Idx), np.take(self.State.Pos_Y, Idx), np.take(self.State.Pos_Angle, Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetPositionXYAngleByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetPosIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetPositionXYAngle(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetDistIdx(self, Name):" << endl;
    ofs << in+ in+ "return self.State.GetDistNames().index(Name)" << endl;
    ofs << endl;

    ofs << in+ "def GetDistributionByName(self, Name):" << endl;
    ofs << in+ in+ "Idx = self.GetDistIdx(Name)" << endl;
    ofs << in+ in+ "return self.GetDistribution(Idx)" << endl;
    ofs << endl;

    ofs << in+ "def GetCountFromDistribution(self, Dist_Idx, X, Y):" << endl;
    ofs << in+ in+ "return self.State.Dist_All[Dist_Idx, X.astype(int), Y.astype(int)]" << endl;
    ofs << endl;

    ofs << in+ "def GetCountFromDistributionByNameAndPos(self, NameOfDist, NameOfPos):" << endl;
    // temporary code
    ofs << in+ in+ "X, Y = self.GetPositionXYByName(NameOfPos)" << endl;
    ofs << in+ in+ "Dist_Idx = self.GetDistIdx(NameOfDist)" << endl;
    ofs << in+ in+ "return self.GetCountFromDistribution(Dist_Idx, X, Y)" << endl;
    ofs << endl;

    ofs << in+ "def ApplyCountToDistribution(self, Dist, X, Y, Count):" << endl;
    ofs << in+ in+ "Dist[X.astype(int), Y.astype(int)] += Count" << endl;
    ofs << endl;

    ofs << in+ "def ApplyCountToDistributionByNameAndPos(self, NameOfDist, NameOfPos):" << endl;
    ofs << in+ in+ "Count = self.GetCountByName" << endl;
    ofs << endl;


    ofs << in+ "# Temporary routines" << endl;

    ofs << in+ "def OverElongationCorrection(self, Len_Elongated, Max):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ "Len_Over = np.where(Len_Elongated > Max, Len_Elongated - Max, 0)" << endl;
    ofs << in+ in+ "return Len_Elongated - Len_Over" << endl;
    ofs << endl;

    ofs << in+ "def BuildingBlockConsumption(self, Freq, N_Elongated_PerSpecies):" << endl;
    ofs << in+ in+ "Raw = SimF.DetermineAmountOfBuildingBlocks(Freq, N_Elongated_PerSpecies)" << endl;
    ofs << in+ in+ "Rounded = np.around(Raw)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Discrepancy handling" << endl;
    ofs << in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << in+ in+ "Discrepancy = np.sum(Rounded) - N_Elongated" << endl;

    ofs << in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ "Sets, Remainder = np.divmod(Discrepancy, NUniq_BuildingBlocks)" << endl;
    ofs << in+ in+ "return Rounded + np.ones(NUniq_BuildingBlocks) * np.int32(Sets) + np.concatenate((np.ones(np.int32(np.round(Remainder))), np.zeros(np.int32(np.around(NUniq_BuildingBlocks - Remainder)))))" << endl;
    ofs << endl;

    ofs << in+ "# Polymerase Reaction related" << endl;
    ofs << in+ "def Initiation(self, Len_Template, Len_Target, Idx_Pol, Idx_Template, Idx_TemplateSubset, Weight, PolThreshold):" << endl;
    ofs << in+ in+ "# Get available, active polymerase count - TO BE UPDATED with more regulatory algorithms" << endl;
    ofs << in+ in+ "Count_Pol = self.GetCount(Idx_Pol)" << endl;
    ofs << in+ in+ "Count_Pol_Active = np.floor_divide(Count_Pol, 2).astype(int)" << endl;
    ofs << in+ in+ "Count_Pol_Occupied = np.count_nonzero(np.where(Len_Target != -1, 1, 0)) * PolThreshold" << endl;
    ofs << in+ in+ "Count_Pol_Avail = Count_Pol_Active - Count_Pol_Occupied" << endl;
    ofs << in+ in+ "Count_Pol_FunctionalUnit = np.floor_divide(Count_Pol_Avail, PolThreshold)" << endl;
    ofs << in+ in+ "Count_Pol_Avail = np.where(Count_Pol_FunctionalUnit > 0, Count_Pol_FunctionalUnit, 0)[0, 0]" << endl;
    ofs << endl;

    ofs << in+ in+ "# Get final initiation weight by applying initiation site count" << endl;
    ofs << in+ in+ "Count_Template_Complete = self.GetCount(Idx_Template)" << endl;
    ofs << in+ in+ "Count_Template_Nascent = np.count_nonzero(np.where(Len_Template != -1, 1, 0), axis=0)   # Assumption: each nascent template has one highly efficient initiation site" << endl;
    ofs << in+ in+ "Count_TemplateSubset_Nascent = np.take(Count_Template_Nascent, Idx_TemplateSubset)   # May need to update np.take with fancy indexing [:, idx]" << endl;
    ofs << in+ in+ "Count_InitiationSite = Count_Template_Complete + Count_TemplateSubset_Nascent" << endl;
    ofs << in+ in+ "Weight_Initiation = Count_InitiationSite * Weight " << endl;
    ofs << endl;

    ofs << in+ in+ "# Get randomly selected target indices" << endl;
    ofs << in+ in+ "Idx_Selected = SimF.PickRandomIdx(Count_Pol_Avail, Idx_Template, Weight_Initiation)" << endl;
    ofs << in+ in+ "Len_Target_Initiated = SimF.InsertZeroIntoNegOneElementInLenMatrix(Len_Target, Idx_Selected)" << endl;
    ofs << in+ in+ "# Export Data" << endl;
    ofs << in+ in+ "# N_Initiated" << endl;
    ofs << endl;

    ofs << in+ in+ "return Len_Target_Initiated" << endl;
    ofs << endl;

    ofs << in+ "def Elongation(self, Len, Max, Rate, Weight, Freq, Idx_PolSub, Idx_BB):" << endl;
    ofs << in+ in+ "NUniq_BuildingBlocks = Freq.shape[1]" << endl;
    ofs << in+ in+ "NUniq_Species = Freq.shape[0]" << endl;
    ofs << endl;

//    ofs << in+ in+ "dLength = np.matmul(SMatrix,Rate)
    ofs << in+ in+ "dLength = self.ApplySimTimeResolution(Rate)   # this is not necessarily true based on the reaction input" << endl;
    ofs << in+ in+ "Len_Elongated = np.where(Len >= 0, Len + dLength, Len)" << endl;
    ofs << endl;

    ofs << in+ in+ "Len_Trimmed = self.OverElongationCorrection(Len_Elongated, Max)" << endl;

    ofs << in+ in+ "N_Elongated_PerSpecies = np.asmatrix(np.sum(Len_Trimmed - Len, axis=0))   # This step loses shape for some reason, hence apply matrix again" << endl;
    ofs << in+ in+ "N_Elongated = np.sum(N_Elongated_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ "Consumed_BB = self.BuildingBlockConsumption(Freq, N_Elongated_PerSpecies)" << endl;

    ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
    ofs << in+ in+ "self.AddTodCount(Idx_BB, -Consumed_BB)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Update dCount for Polymerase Reaction Substrates (To be updated by the reaction matrix form" << endl;
    ofs << in+ in+ "self.AddTodCount(Idx_PolSub, N_Elongated)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Export Data" << endl;
    ofs << in+ in+ "# N_Elongated" << endl;

    ofs << in+ in+ "return Len_Trimmed" << endl;
    ofs << endl;

    ofs << in+ "def Termination(self, Len, Max, Idx_Target):   # Some polymerization process may not have max" << endl;
    ofs << in+ in+ "Bool_Completed = (Len == Max)" << endl;
    ofs << in+ in+ "N_Completed_PerSpecies = np.sum(Bool_Completed, axis=0)" << endl;
    ofs << in+ in+ "N_Completed = np.sum(N_Completed_PerSpecies)" << endl;
    ofs << in+ in+ "Len_Completed = np.where(Bool_Completed, -1, Len)" << endl;

    ofs << in+ in+ "# Update dCount for BuildingBlocks" << endl;
    ofs << in+ in+ "self.AddTodCount(Idx_Target, N_Completed_PerSpecies)" << endl;
    ofs << endl;

    ofs << in+ in+ "# Export Data" << endl;
    ofs << in+ in+ "# N_Completed" << endl;

    ofs << in+ in+ "return Len_Completed" << endl;
    ofs << endl;

    ofs << in+ "# Debugging tools" << endl;
    ofs << in+ "def Debug_PrintNames(self, IdxArray):" << endl;
    ofs << in+ in+ "ListOfNames = list()" << endl;
    ofs << in+ in+ "if IdxArray.dim == 1:" << endl;
    ofs << in+ in+ in+ "IdxArray = IdxArray.tolist()" << endl;
    ofs << in+ in+ "if IdxArray.dim == 2:" << endl;
    ofs << in+ in+ in+ "IdxArray = IdxArray.tolist()[0]" << endl;
    ofs << in+ in+ "for Idx in IdxArray.tolist()[0]:" << endl;
    ofs << in+ in+ in+ "ListOfNames.append(self.State.GetMolNames()[Idx])" << endl;
    ofs << in+ in+ "print(ListOfNames)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_SetIdxMoleculesToTrack(self):" << endl;
    ofs << in+ in+ "# Add a list of molecules to track for debugging every simulation step" << endl;
    ofs << in+ in+ "#Debug_Names_Molecules = []" << endl; // TODO: take input from command line
    ofs << in+ in+ "Debug_Names_Molecules = ['Am', 'AmL', 'L', 'qL', 'pc_qL']" << endl; // TODO: take input from command line
    ofs << endl;
    ofs << in+ in+ "if Debug_Names_Molecules == []:" << endl;
    ofs << in+ in+ in+ "Debug_Names_Molecules = self.State.GetMolNames()" << endl;
    ofs << endl;
    ofs << in+ in+ "for Name in Debug_Names_Molecules:" << endl;
    ofs << in+ in+ in+ "self.Debug_Idx_Molecules.append(self.State.GetMolNames().index(Name))" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCounts(self, Switch, Idx_Pos=0):" << endl;
    ofs << in+ in+ "for Idx_Pos in self.Debug_Idx_Pos:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCountsForSinglePosition(DisplayCount, Idx_Pos)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCountsForSinglePosition(self, Switch, Idx_Pos=0):" << endl;
    ofs << in+ in+ "for Idx_Mol in self.Debug_Idx_Molecules:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCount(Idx_Mol, Idx_Pos)" << endl;
    ofs << in+ in+ "print()" << endl;
    ofs << in+ in+ "if Switch:" << endl;
    ofs << in+ in+ in+ "for Idx_Mol in self.Debug_Idx_Molecules:" << endl;
    ofs << in+ in+ in+ in+ "self.Debug_PrintdCount(Idx_Mol, Idx_Pos)" << endl;
    ofs << in+ in+ in+ "print()" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintSimStepTime(self):" << endl;
    ofs << in+ in+ "Time = self.GetSimTime()" << endl;
    ofs << in+ in+ "print(self.SimStep, '(', round(Time,3), 's)', end='\t| ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCount(self, Idx_Mol, Idx_Pos=0):" << endl;
    ofs << in+ in+ "print(' ' + self.State.GetMolNames()[Idx_Mol], end=': ')" << endl;
    ofs << in+ in+ "print('{:010e}'.format(self.Debug_ApplyUnit(self.State.Count_All[Idx_Pos][Idx_Mol])), self.UnitTxt, end=' | ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintdCount(self, Idx_Mol, Idx_Pos=0):" << endl;
    ofs << in+ in+ "print('d' + self.State.GetMolNames()[Idx_Mol], end=': ')" << endl;
    ofs << in+ in+ "print('{:010e}'.format(self.Debug_ApplyUnit(self.State.dCount_All[Idx_Pos][Idx_Mol])), self.UnitTxt, end=' | ')" << endl;
    ofs << endl;

    ofs << in+ "def Debug_SetUnit(self, Input):" << endl;
    ofs << in+ in+ "self.UnitTxt = Input" << endl;
    ofs << in+ in+ "if Unit == 'nM':" << endl;
    ofs << in+ in+ in+ "self.Unit = 1e-9 * " << Numbers::GetAvogadroStr() << endl;
    ofs << in+ in+ "elif Unit == 'uM':" << endl;
    ofs << in+ in+ in+ "self.Unit = 1e-6 * " << Numbers::GetAvogadroStr() << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "self.Unit = 1" << endl;
    ofs << endl;

    ofs << in+ "def Debug_ApplyUnit(self, Value):" << endl;
    ofs << in+ in+ "return Value / self.Unit" << endl;
    ofs << endl;

    // Debugging for spatial simulation
    ofs << in+ "def Debug_SetIdxDistAndPosToTrack(self):" << endl;
    ofs << in+ in+ "# Add a list of distributions at and around selected positions to track for debugging every simulation step" << endl;
    ofs << in+ in+ "Debug_Names_Positions = []" << endl; // TODO: take input from command line
    ofs << in+ in+ "if Debug_Names_Positions == []:" << endl;
    ofs << in+ in+ in+ "self.Debug_Idx_Pos = list(range(len(self.State.GetPosNames())))" << endl;
    ofs << in+ in+ "else:   # does not properly cover all cases due to redundant names allowed through for loops" << endl;
    ofs << in+ in+ in+ "for Name in Debug_Names_Positions:" << endl;
    ofs << in+ in+ in+ in+ "self.Debug_Idx_Pos.append(self.State.GetPosNames().index(Name))" << endl;
    ofs << endl;
    ofs << in+ in+ "Debug_Names_Distributions = []" << endl; // TODO: take input from command line
    ofs << in+ in+ "if Debug_Names_Distributions == []:" << endl;
    ofs << in+ in+ in+ "Debug_Names_Distributions = self.State.GetDistNames()" << endl;
    ofs << in+ in+ "for Name in Debug_Names_Distributions:" << endl;
    ofs << in+ in+ in+ "self.Debug_Idx_Dist.append(self.State.GetDistNames().index(Name))" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintDistributions(self):" << endl;
    ofs << in+ in+ "for Idx_Pos in self.Debug_Idx_Pos:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintDistributionsForSinglePosition(Idx_Pos)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintCountsAndDistributions(self):" << endl;
    ofs << in+ in+ "for Idx_Pos in self.Debug_Idx_Pos:" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintDistributionsForSinglePosition(Idx_Pos)" << endl;
    ofs << in+ in+ in+ "self.Debug_PrintCountsForSinglePosition(DisplayCount, Idx_Pos)" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintDistributionsForSinglePosition(self, Idx_Pos):" << endl;
    ofs << in+ in+ "X = self.State.Pos_X[Idx_Pos].astype(int)" << endl;
    ofs << in+ in+ "Y = self.State.Pos_Y[Idx_Pos].astype(int)" << endl;
    ofs << in+ in+ "print('Position #' + str(Idx_Pos) + ' Name: ' + self.State.GetPosNames()[Idx_Pos] + ' @ ({}, {})'.format(X, Y))" << endl;
    ofs << in+ in+ "DistNames = ''" << endl;
    ofs << in+ in+ "DistValues = dict()" << endl;
    ofs << in+ in+ "DistValues[0] = ''" << endl;
    ofs << in+ in+ "DistValues[1] = ''" << endl;
    ofs << in+ in+ "DistValues[2] = ''" << endl;
    ofs << in+ in+ "for Idx_Dist in self.Debug_Idx_Dist:" << endl;
    ofs << in+ in+ in+ "DistNames, DistValues = self.Debug_PrintDistribution(DistNames, DistValues, Idx_Dist, X, Y)" << endl;
    ofs << in+ in+ "print(DistNames, '\\n', DistValues[0], '\\n', DistValues[1], '\\n', DistValues[2])" << endl;
    ofs << endl;

    ofs << in+ "def Debug_PrintDistribution(self, DistNames, DistValues, Idx_Dist, X, Y):" << endl;
    ofs << in+ in+ "DistNames += self.State.GetDistNames()[Idx_Dist] + ' (' + self.UnitTxt + ')' + '\\t\\t\\t\\t\\t\\t\\t\\t'" << endl;
    ofs << in+ in+ "DistValues[0] += self.Debug_GetArrayInString(self.Debug_ApplyUnit(self.State.Dist_All[Idx_Dist][X-1:X+2, Y-1])) + '\\t'" << endl;
    ofs << in+ in+ "DistValues[1] += self.Debug_GetArrayInString(self.Debug_ApplyUnit(self.State.Dist_All[Idx_Dist][X-1:X+2, Y]))   + '\\t'" << endl;
    ofs << in+ in+ "DistValues[2] += self.Debug_GetArrayInString(self.Debug_ApplyUnit(self.State.Dist_All[Idx_Dist][X-1:X+2, Y+1])) + '\\t'" << endl;
    ofs << in+ in+ "return DistNames, DistValues" << endl;
    ofs << endl;

    ofs << in+ "def Debug_GetArrayInString(self, Array_1D):" << endl;
    ofs << in+ in+ "Array_1D_Str = ''" << endl;
    ofs << in+ in+ "for Idx in range(Array_1D.shape[0]):" << endl;
    ofs << in+ in+ in+ "Array_1D_Str += SimF.SciFloat(Array_1D[Idx]) + ' '" << endl;
    ofs << in+ in+ "return Array_1D_Str" << endl;
    ofs << endl;



    // class FDataManager
    ofs << "class FDataManager:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "# SimOut" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.DataBuffer = list()" << endl;

    ofs << in+ "def SetLegend(self, InLegend):" << endl;
    ofs << in+ in+ "self.Legend = InLegend" << endl;
    ofs << endl;

    ofs << in+ "def Add(self, InData):" << endl;
    ofs << in+ in+ "self.DataBuffer.append(InData)" << endl;
    ofs << endl;

    ofs << in+ "def SaveToTsvFile(self, Legend, Data, InFileName):" << endl;
    ofs << in+ in+ "with open(InFileName, 'w', newline='', encoding='utf-8') as OutFile:" << endl;
    ofs << in+ in+ in+ "TsvWriter = csv.writer(OutFile, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ "if Legend:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(Legend)" << endl;
    ofs << in+ in+ in+ "for Row in Data:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(np.array(Row).flatten().tolist())" << endl;
    ofs << endl;

    ofs << in+ "def SaveToNpyFile(self, Data, InFileName):" << endl;
    ofs << in+ in+ "np.save('%s.npy' % InFileName, Data)" << endl;
    ofs << endl;

    ofs << in+ "def SaveCountToFile(self, InFileName):" << endl;
    ofs << in+ in+ "self.SaveToTsvFile(self.Legend, self.DataBuffer, InFileName)" << endl;
    ofs << endl;

    // BODY
    ofs << "def main():   # add verbose" << endl;
    ofs << endl;

    // Instantiate Objects
    ofs << in+ "State = FState()" << endl;
    ofs << in+ "Data = FDataset()" << endl;
    ofs << in+ "DataManager = FDataManager()" << endl;
    ofs << in+ "Simulation = FSimulation(State, Data, DataManager)" << endl;
    ofs << endl;

    // Simulation Module
    ofs << in+ "Simulation.Initialize(N_SimSteps, SimStepTimeResolution)" << endl;
    ofs << in+ "Simulation.Run(Spatial=0) # 0: WithoutSpatialSimulation, 1: WithSpatialSimulation" << endl;
    ofs << endl;
    ofs << in+ "DataManager.SaveCountToFile('" << Option.SimResultFile.c_str() << "')" << endl;
    ofs << endl;

    // MAIN
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in + "parser = ArgumentParser()" << endl;
    ofs << in + "parser.add_argument('--save-fig', dest='save_fname', type=str, help='Save figure to file')" << endl;
    ofs << in + "args = parser.parse_args()" << endl;
    ofs << in + "main()" << endl;
    ofs << in + "plot.main(args.save_fname)" << endl;
    ofs << endl;

    std::cout << "  Simulation program has been generated: ";
}

void FWriter::SimVis2D()
{
    std::cout << "Generating SimVis2D..." << std::endl;

    // write SimVis2D.py
    std::ofstream ofs(Option.SimVis2DFile.c_str());
    std::string endl = "\n";

    ofs << "import sys" << endl;
    ofs << "import pygame" << endl;
    ofs << "import random" << endl;
    ofs << "from datetime import datetime" << endl;
    ofs << "import SimModule" << endl;
    ofs << "import SimFunctions as SimF" << endl;
    ofs << "import numpy as np" << endl;
    ofs << "import math" << endl;
    ofs << endl;
    ofs << "# Colors" << endl;
    ofs << "BLACK = (0, 0, 0)" << endl;
    ofs << "WHITE = (255, 255, 255)" << endl;
    ofs << "GRAY1 = (230, 230, 230)" << endl;
    ofs << "GRAY2 = (210, 210, 210)" << endl;
    ofs << "GRAY3 = (150, 150, 150)" << endl;
    ofs << "GRAY4 = (100, 100, 100)" << endl;
    ofs << "YELLOW_FAINT = (200, 200, 150)" << endl;
    ofs << "RED = (255, 0, 0)" << endl;
    ofs << "GREEN = (0, 255, 0)" << endl;
    ofs << "YELLOW = (255, 255, 0)" << endl;
    ofs << "BLUE = (0, 0, 255)" << endl;
    ofs << "MAGENTA = (255, 0, 255)" << endl;
    ofs << "CYAN = (0, 255, 255)" << endl;
    ofs << endl;
    ofs << "# Global variables" << endl;
    ofs << "NA = 6.0221409e+23" << endl;
    ofs << "pi = np.pi" << endl;
    ofs << "uM = 1e-6" << endl;
    ofs << "nM = 1e-9" << endl;
    ofs << "UnitTxt = ''" << endl;
    ofs << endl;
    ofs << "# Determine global unit" << endl;
    ofs << "Unit = nM" << endl;
    ofs << "if Unit == nM:" << endl;
    ofs << in+ "UnitTxt = 'nM'" << endl;
    ofs << "elif Unit == uM:" << endl;
    ofs << in+ "UnitTxt = 'uM'" << endl;
    ofs << endl;

    auto MolLoc = Context.GetSubList_LocationList("Molecule");

    // Define indices for spatially distributed molecules and containers
    ofs << "GlucoseName = 'L'" << endl;

    // Pathway dependent key molecules to monitor
    ofs << "HomeostasisMolName = [";
    for (auto location : MolLoc) {
        if      (location->Name == "L")  { ofs << "'" << "Am" << "', "; }
        else if (location->Name == "qL") { ofs << "'" << "qAm" << "', "; }
    } ofs << "]" << endl;
    ofs << endl;

    ofs << "# Utilities" << endl;
    ofs << "def GetColorGradient(Fade, baseColor=None):" << endl;
    ofs << in+ "assert Fade <= 255, 'ERROR: Fade factor is higher than 255!'" << endl;
    ofs << in+ "if baseColor == 'Blue':" << endl;
    ofs << in+ in+ "return (Fade, Fade, 255)" << endl;
    ofs << in+ "else:  # Assumes White" << endl;
    ofs << in+ in+ "return (Fade, Fade, Fade)" << endl;
    ofs << endl;
    ofs << "# pygame" << endl;
    ofs << "pygame.init()" << endl;
    ofs << endl;
    ofs << "# Load model" << endl;
    ofs << "State = SimModule.FState()" << endl;
    ofs << "Data = SimModule.FDataset()" << endl;
    ofs << "DataManager = SimModule.FDataManager()" << endl;
    ofs << "SimM = SimModule.FSimulation(State, Data, DataManager)" << endl;
    ofs << endl;
    ofs << "# Initialize model" << endl;
    ofs << "SimM.Initialize()" << endl;
    ofs << endl;
    ofs << "Screen_Size = W_S, H_S = " << "SimM.GetDistWidth()" << "," << "SimM.GetDistHeight()" << endl;
    ofs << "Screen = pygame.display.set_mode(Screen_Size)" << endl;
    ofs << endl;
    ofs << "LEFT = 0" << endl;
    ofs << "MID_X = W_S / 2" << endl;
    ofs << "RIGHT = W_S" << endl;
    ofs << endl;
    ofs << "TOP = 0" << endl;
    ofs << "MID_Y = H_S / 2" << endl;
    ofs << "BOTTOM = H_S" << endl;
    ofs << endl;
    ofs << "CenterTop = (MID_X, TOP)" << endl;
    ofs << "Center = (MID_X, MID_Y)" << endl;
    ofs << endl;
    ofs << "def AddPos(A, B):" << endl;
    ofs << in+ "X_A, Y_A = A" << endl;
    ofs << in+ "X_B, Y_B = B" << endl;
    ofs << in+ "return (X_A + X_B, Y_A + Y_B)" << endl;
    ofs << endl;
    ofs << "def GetMidPoint(A, B):" << endl;
    ofs << in+ "X_A, Y_A = A" << endl;
    ofs << in+ "X_B, Y_B = B" << endl;
    ofs << in+ "return ((X_A + X_B) / 2, (Y_A + Y_B) / 2)" << endl;
    ofs << endl;
    ofs << "# # Transparent control board" << endl;
    ofs << "# ControlBoard = pygame.Surface(Screen_Size, pygame.SRCALPHA)" << endl;
    ofs << "# ControlBoard.fill((0, 0, 0, 255))" << endl;
    ofs << endl;
    ofs << "Title = 'Vis2D'" << endl;
    ofs << "pygame.display.set_caption(Title)" << endl;
    ofs << endl;
    ofs << "Font_Sans = pygame.font.Font('freesansbold.ttf', 20)" << endl;
    ofs << "Font_Monospace = pygame.font.SysFont('monospace', 18, True)" << endl;
    ofs << endl;

    ofs << "class FEnvironment:" << endl;
    ofs << in+ "def __init__(self, InX=W_S/2, InY=H_S/2, InShape='circle', InRadius=W_S*0.5, InThickness=10):" << endl;
    ofs << in+ in+ "self.X = InX" << endl;
    ofs << in+ in+ "self.Y = InY" << endl;
    ofs << in+ in+ "self.Shape = InShape" << endl;
    ofs << in+ in+ "self.Radius = int(InRadius)" << endl;
    ofs << in+ in+ "self.Thickness = InThickness" << endl;
    ofs << in+ in+ "self.TransparentCircleArea = None" << endl;
    ofs << endl;
    ofs << in+ "def Draw(self, shape=None):" << endl;
    ofs << in+ in+ "if shape == 'circle':" << endl;
    ofs << in+ in+ in+ "pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)" << endl;
    ofs << in+ in+ "elif shape == 'lining':" << endl;
    ofs << in+ in+ in+ "# pygame.draw.circle(Screen, YELLOW_FAINT, (self.X, self.Y), self.Radius)" << endl;
    ofs << in+ in+ in+ "pygame.draw.circle(Screen, GRAY4, (self.X, self.Y), self.Radius, self.Thickness)" << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "pygame.draw.rect(Screen, YELLOW_FAINT, ((0, 0), (W_S, H_S)))" << endl;
    ofs << endl;
    ofs << in+ "def DrawTransparentArea(self):" << endl;
    ofs << in+ in+ "self.TransparentCircleArea = pygame.Surface((self.Radius * 2, self.Radius * 2), pygame.SRCALPHA)" << endl;
    ofs << in+ in+ "self.TransparentCircleArea.fill((255, 255, 255, 255))" << endl;
    ofs << in+ in+ "pygame.draw.circle(self.TransparentCircleArea, (0, 0, 0, 0), (self.Radius, self.Radius), self.Radius)" << endl;
    ofs << endl;
    ofs << "class FOrganism:" << endl;
    ofs << in+ "def __init__(self, InName, InSpecies):" << endl;
    ofs << in+ in+ "self.Name = InName" << endl;
    ofs << in+ in+ "self.Species = InSpecies" << endl;
    ofs << in+ in+ "self.X = None" << endl;
    ofs << in+ in+ "self.Y = None" << endl;
    ofs << in+ in+ "self.Angle = None" << endl;
    ofs << endl;
    ofs << in+ in+ "self.BodyLength = 20" << endl;
    ofs << in+ in+ "self.BodyThickness = 10" << endl;
    ofs << in+ in+ "self.FlagellaLength_Fold = 2" << endl;
    ofs << in+ in+ "self.FlagellaThickness = 3" << endl;
    ofs << endl;
    ofs << in+ in+ "# Memory for display & debugging" << endl;
    ofs << in+ in+ "self.Ligand_Prev = 0" << endl;
    ofs << in+ in+ "self.Am = 0" << endl; // TODO: hardcoded
    ofs << in+ in+ "self.SimCount = 0" << endl;
    ofs << endl;
    ofs << in+ in+ "# Trajectory" << endl;
    ofs << in+ in+ "self.TrajectorySwitch = True" << endl;
    ofs << in+ in+ "self.Trajectory = dict()" << endl;
    ofs << in+ in+ "self.TrajectoryColor = list()" << endl;
    ofs << endl;
    ofs << in+ "def Draw(self):" << endl;
    ofs << in+ in+ "if self.Species == 'Ecoli':" << endl;
    ofs << in+ in+ in+ "Color = ''" << endl;
    ofs << in+ in+ in+ "# pygame.draw.circle(Screen, Color, (self.X, self.Y), 5)" << endl;
    ofs << in+ in+ in+ "dX =  np.cos(self.Angle) * -self.BodyLength" << endl;
    ofs << in+ in+ in+ "dY = -np.sin(self.Angle) * -self.BodyLength" << endl;
    ofs << in+ in+ in+ "X_BodyEnd = self.X + dX" << endl;
    ofs << in+ in+ in+ "Y_BodyEnd = self.Y + dY" << endl;
    ofs << in+ in+ in+ "X_TailEnd = self.X + self.FlagellaLength_Fold * dX" << endl;
    ofs << in+ in+ in+ "Y_TailEnd = self.Y + self.FlagellaLength_Fold * dY" << endl;
    ofs << in+ in+ in+ "for i in range(self.X.size):" << endl;
    ofs << in+ in+ in+ in+ "if i == 0:" << endl;
    ofs << in+ in+ in+ in+ in+ "Color = RED" << endl;
    ofs << in+ in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ in+ "Color = YELLOW" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.line(Screen, Color, (self.X[i], self.Y[i]), (X_BodyEnd[i], Y_BodyEnd[i]), self.BodyThickness)" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.line(Screen, Color, (self.X[i], self.Y[i]), (X_TailEnd[i], Y_TailEnd[i]), self.FlagellaThickness)" << endl;
    ofs << endl;
    ofs << in+ "def SetPosition(self, Position):" << endl;
    ofs << in+ in+ "self.X = Position[0]" << endl;
    ofs << in+ in+ "self.Y = Position[1]" << endl;
    ofs << in+ in+ "self.Angle = Position[2]" << endl;
    ofs << in+ in+ "# Threshold value (Position[3]) is not used" << endl;
    ofs << endl;
    ofs << in+ "def Receptivity(self, N_SimulationsToPass=50):" << endl;
    ofs << in+ in+ "for _ in range(N_SimulationsToPass):" << endl;
    ofs << in+ in+ in+ "SimM.SimLoop_WithoutSpatialSimulation_WithMoleculeDistribution()" << endl;
    ofs << in+ in+ "SimM.SimLoop_WithSpatialSimulation()" << endl;
    ofs << endl;
    ofs << in+ "def ReportStatus(self):" << endl;
    ofs << in+ in+ "# for debugging" << endl;
    ofs << in+ in+ "self.Am = SimM.GetCountByName(HomeostasisMolName[0])" << endl;
//    ofs << in+ in+ "self.qAm = SimM.GetCountByName(HomeostasisMolName[1])" << endl;
    ofs << in+ in+ "self.SimCount += 1" << endl;
    ofs << in+ in+ "Ligand_Now = SimM.GetCountFromDistributionByNameAndPos(GlucoseName, self.Name)[0]" << endl; // TODO: HARDCODED
    ofs << in+ in+ "SimStep = SimM.GetSimStep()" << endl;
    ofs << in+ in+ "Delta = 0" << endl;
    ofs << in+ in+ "if Ligand_Now != 0:" << endl;
    ofs << in+ in+ in+ "Delta = (Ligand_Now - self.Ligand_Prev[0]) / Ligand_Now * 100" << endl;
    ofs << in+ in+ "# print('SimStep {:06d} [Chemotaxis  {:06d}] Ligand:{:.6f} {} ({}{:.4f}%) Am:{:.6f} {} (X:{:.2f}, Y:{:.2f}, {:3.1f} degree)'.format" << endl;
    ofs << in+ in+ "#in+ '   (SimStep, self.SimCount, Ligand_Now / Unit / NA, UnitTxt, ('+' if Delta >= 0 else ''), Delta, self.Am / Unit / NA, UnitTxt , self.X, self.Y, self.Angle / pi * 180))" << endl;
    ofs << endl;
    ofs << in+ "def Homeostasis(self, MolName=[]):" << endl;
    ofs << in+ in+ "SimM.Homeostasis(MolName)   # Input 'Am', 'qAm' here" << endl;
    ofs << in+ in+ "self.Ligand_Prev = SimM.GetCountFromDistributionByNameAndPos(GlucoseName, self.Name)" << endl; // TODO: HARDCODED
    ofs << in+ in+ "self.SetPosition(SimM.GetPositionXYAngleByName(self.Name))" << endl;
    ofs << in+ in+ "self.InitializeTrajectory()" << endl;
    ofs << in+ in+ "self.ReportStatus()" << endl;
    ofs << endl;
    ofs << in+ "def InitializeTrajectory(self):" << endl;
    ofs << in+ in+ "for i in range(self.X.size):" << endl;
    ofs << in+ in+ in+ "self.Trajectory[i] = [(self.X[i], self.Y[i])]" << endl;
    ofs << in+ in+ in+ "self.TrajectoryColor.append(tuple(np.random.randint(0, 255, 3)))" << endl;
    ofs << in+ in+ "self.TrajectoryColor[0] = MAGENTA" << endl;
    ofs << endl;
    ofs << in+ "def AddToTrajectory(self):" << endl;
    ofs << in+ in+ "for i in range(self.X.size):" << endl;
    ofs << in+ in+ in+ "self.Trajectory[i].append((self.X[i], self.Y[i]))" << endl;
    ofs << endl;
    ofs << in+ "def DrawTrajectory(self):" << endl;
    ofs << in+ in+ "for i in range(len(self.Trajectory)):" << endl;
    ofs << in+ in+ in+ "pygame.draw.aalines(Screen, self.TrajectoryColor[i], False, self.Trajectory[i])" << endl;
    ofs << endl;
    ofs << "class FMolecule:" << endl;
    ofs << in+ "def __init__(self, InName, InX, InY):" << endl;
    ofs << in+ in+ "self.Name = InName" << endl;
    ofs << in+ in+ "self.X = InX" << endl;
    ofs << in+ in+ "self.Y = InY" << endl;
    ofs << in+ in+ "self.DiffusionFactor = 2" << endl;
    ofs << in+ in+ "self.SpaceFactor = 50" << endl;
    ofs << in+ in+ "self.Color = ''" << endl;
    ofs << in+ in+ "self.Pattern = ''" << endl;
    ofs << in+ in+ "self.MaxBrightness = None" << endl;
    ofs << in+ in+ "self.NormalizationType = ''" << endl;
    ofs << endl;
    ofs << in+ in+ "# Heatmap Drawing" << endl;
    ofs << in+ in+ "self.ReductionFactor = 5" << endl;
    ofs << in+ in+ "self.Max = 0" << endl;
    ofs << in+ in+ "self.InitializeHeatmapMax()" << endl;
    ofs << endl;
    ofs << in+ in+ "# Particle Drawing" << endl;
    ofs << in+ in+ "self.Particle_N = 300" << endl;
    ofs << in+ in+ "self.Particle_PerLayer = 3" << endl;
    ofs << in+ in+ "self.Particle_Radius = 2" << endl;
    ofs << in+ in+ "self.Particle_SpreadFactor = 1.2" << endl;
    ofs << in+ in+ "self.Particle_XY_Static = []" << endl;
    ofs << in+ in+ "self.InitializeStaticParticles()" << endl;
    ofs << endl;
    ofs << in+ "def InitializeHeatmapMax(self):" << endl;
    ofs << in+ in+ "self.Max = np.max(SimM.GetDistributionByName(self.Name))" << endl;
    ofs << endl;
    ofs << in+ "def InitializeStaticParticles(self):" << endl;
    ofs << in+ in+ "for i in range(int(self.Particle_N / self.Particle_PerLayer)):" << endl;
    ofs << in+ in+ in+ "for j in range(self.Particle_PerLayer):" << endl;
    ofs << in+ in+ in+ in+ "X = self.X + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))" << endl;
    ofs << in+ in+ in+ in+ "Y = self.Y + random.randint(-int(i ** self.Particle_SpreadFactor), int(i ** self.Particle_SpreadFactor))" << endl;
    ofs << in+ in+ in+ in+ "self.Particle_XY_Static.append((X, Y))" << endl;
    ofs << endl;
    ofs << in+ "def Reposition(self):" << endl;
    ofs << in+ in+ "self.X = random.randint(W_S * 2 / 5, W_S * 3 / 5)" << endl;
    ofs << in+ in+ "self.Y = random.randint(H_S * 2 / 5, H_S * 3 / 5)" << endl;
    ofs << in+ in+ "self.Particle_XY_Static = []" << endl;
    ofs << in+ in+ "self.InitializeStaticParticles()" << endl;
    ofs << endl;
    ofs << in+ "def Move(self, dX, dY):" << endl;
    ofs << in+ in+ "New_Particle_XY_Static = []" << endl;
    ofs << in+ in+ "for (X, Y) in self.Particle_XY_Static:" << endl;
    ofs << in+ in+ in+ "X += dX" << endl;
    ofs << in+ in+ in+ "Y += dY" << endl;
    ofs << in+ in+ in+ "New_Particle_XY_Static.append((X, Y))" << endl;
    ofs << in+ in+ "self.Particle_XY_Static = New_Particle_XY_Static" << endl;
    ofs << endl;
    ofs << in+ "def Draw(self, Data, threshold=0):" << endl;
    ofs << in+ in+ "if self.Pattern == 'spots':" << endl;
    ofs << in+ in+ in+ "Max = np.max(Data)" << endl;
    ofs << in+ in+ in+ "if Max == 0:" << endl;
    ofs << in+ in+ in+ in+ "return"<< endl;
    ofs << in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ "Data = SimF.Normalize_Linear(Data)" << endl;
    ofs << in+ in+ in+ in+ "CoordsToDraw = np.where(Data > threshold)" << endl;
    ofs << in+ in+ in+ in+ "for X, Y, Value in zip(CoordsToDraw[0], CoordsToDraw[1], Data[CoordsToDraw]):" << endl;
    ofs << in+ in+ in+ in+ in+ "intensity = math.floor(Value * self.MaxBrightness)" << endl;
    ofs << in+ in+ in+ in+ in+ "if intensity > self.MaxBrightness:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "intensity = self.MaxBrightness" << endl;
    ofs << in+ in+ in+ in+ in+ "color = self.GetColor(intensity)"<< endl;
    ofs << in+ in+ in+ in+ in+ "pygame.draw.circle(Screen, color, (X, Y), self.Particle_Radius)" << endl;
    ofs << endl;
    ofs << in+ in+ "elif self.Pattern == 'heatmap':" << endl;
    ofs << in+ in+ in+ "Max = np.max(Data)" << endl;
    ofs << in+ in+ in+ "if Max == 0:" << endl;
    ofs << in+ in+ in+ in+ "return"<< endl;
    ofs << in+ in+ in+ "else:" << endl;
    ofs << in+ in+ in+ in+ "if self.NormalizationType == 'linear':" << endl;
    ofs << in+ in+ in+ in+ in+ "Data = SimF.Normalize_Linear(Data)" << endl;
    ofs << in+ in+ in+ in+ "if self.NormalizationType == 'log':" << endl;
    ofs << in+ in+ in+ in+ in+ "Data = SimF.Normalize_P1Log(Data)" << endl;
    ofs << in+ in+ in+ in+ "for x in range(0, Data.shape[0], self.ReductionFactor):" << endl;
    ofs << in+ in+ in+ in+ in+ "for y in range(0, Data.shape[1], self.ReductionFactor):" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "PercentMolLevel = Data[x][y]" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "if PercentMolLevel < threshold or PercentMolLevel == 0:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ in+ "continue" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "intensity = math.floor(PercentMolLevel * self.MaxBrightness)" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "if intensity > self.MaxBrightness:" << endl;
    ofs << in+ in+ in+ in+ in+ in+ in+ "intensity = self.MaxBrightness" << endl;
    ofs << in+ in+ in+ in+ in+ in+ "color = self.GetColor(intensity)"<< endl;
    ofs << in+ in+ in+ in+ in+ in+ "pygame.draw.rect(Screen, color, ((x, y), (self.ReductionFactor, self.ReductionFactor)))" << endl;
    ofs << endl;
    ofs << in+ in+ "elif self.Pattern == 'particle':" << endl;
    ofs << in+ in+ in+ "for XY in self.Particle_XY_Static:" << endl;
    ofs << in+ in+ in+ in+ "pygame.draw.circle(Screen, BLUE, XY, self.Particle_Radius)" << endl;
    ofs << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "assert True, 'Unsupported molecule distribution pattern for drawing: %s' % self.Pattern" << endl;
    ofs << endl;

    ofs << in+ "def GetColor(self, Intensity):" << endl;
    ofs << in+ in+ "if self.Color == 'Yellow':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness, self.MaxBrightness, self.MaxBrightness - Intensity)" << endl;
    ofs << in+ in+ "if self.Color == 'Blue':" << endl;
    ofs << in+ in+ in+ "return (self.MaxBrightness - Intensity, self.MaxBrightness - Intensity, self.MaxBrightness)" << endl;
    ofs << endl;

    ofs << in+ "def SetColor(self, Color, MaxBrightness):" << endl;
    ofs << in+ in+ "self.Color = Color" << endl;
    ofs << in+ in+ "self.MaxBrightness = MaxBrightness" << endl;
    ofs << endl;

    ofs << in+ "def SetPattern(self, Pattern, NormalizationType, ReductionFactor):" << endl;
    ofs << in+ in+ "self.Pattern = Pattern" << endl;
    ofs << in+ in+ "self.NormalizationType = NormalizationType" << endl;
    ofs << in+ in+ "self.ReductionFactor = ReductionFactor" << endl;
    ofs << endl;

    ofs << "class FControl:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "# self.FPS = 30" << endl;
    ofs << in+ in+ "self.MovementResolution = 30" << endl;
    ofs << in+ in+ "self.MessageWelcome = 'Welcome to Bacterial Chemotaxis!'" << endl;
    ofs << in+ in+ "self.Pos_Welcome = GetMidPoint(Center, CenterTop)" << endl;
    ofs << in+ in+ "self.Message = ''" << endl;
    ofs << in+ in+ "self.MessageTimer = 3000" << endl;
    ofs << endl;
    ofs << in+ in+ "self.InstructionSwitch = False" << endl;
    ofs << in+ in+ "self.Instructions = {" << endl;
    ofs << in+ in+ in+ "'T' : 'Display Trajectory Switch'," << endl;
    ofs << in+ in+ in+ "'I' : 'Display Instruction Switch'," << endl;
    ofs << in+ in+ in+ "'S' : 'Display Score Switch'," << endl;
    ofs << in+ in+ in+ "'A' : 'Display Status Switch'," << endl;
    ofs << in+ in+ in+ "'P' : 'Pause Visualization'," << endl;
    ofs << in+ in+ "}" << endl;
    ofs << in+ in+ "self.InstructionText = ''" << endl;
    ofs << in+ in+ "self.SetInstructionText()" << endl;
    ofs << endl;
    ofs << in+ in+ "self.MoleculeGradientText = ''" << endl;
    ofs << in+ in+ "self.MoleculeGradientColor = list()" << endl;
    ofs << endl;
    ofs << in+ in+ "self.Score = 0" << endl;
    ofs << in+ in+ "self.ScoreSwitch = False" << endl;
    ofs << endl;
    ofs << in+ in+ "self.Time = 0" << endl;
    ofs << in+ in+ "self.ScoreSwitch = False" << endl;
    ofs << endl;
    ofs << in+ in+ "# DK - debugging purposes" << endl;
    ofs << in+ in+ "self.StatusSwitch = True" << endl;
    ofs << in+ in+ "# self.StatusSwitch = False" << endl;
    ofs << endl;
    ofs << in+ in+ "# Pause" << endl;
    ofs << in+ in+ "self.MessagePause = 'PAUSE'" << endl;
    ofs << in+ in+ "self.PauseSwitch = False" << endl;
    ofs << endl;
    ofs << in+ "def SetInstructionText(self):" << endl;
    ofs << in+ in+ "self.InstructionText = 'Instructions \\" << "n'" << endl;
    ofs << in+ in+ "for Key, Value in self.Instructions.items():" << endl;
    ofs << in+ in+ in+ "Space = ' ' * (6 - len(Key))" << endl;
    ofs << in+ in+ in+ "self.InstructionText = self.InstructionText + '  ' + Key + Space + ': ' + Value + '\\" << "n'" << endl;
    ofs << endl;
    ofs << in+ "def DisplayInstructions(self):" << endl;
    ofs << in+ in+ "TextLines = self.InstructionText.splitlines()" << endl;
    ofs << in+ in+ "Height = Font_Monospace.get_linesize() + 2" << endl;
    ofs << in+ in+ "X, Y = Screen.get_rect().topleft" << endl;
    ofs << in+ in+ "Color = None" << endl;
    ofs << in+ in+ "for i, TextLine in enumerate(TextLines):" << endl;
    ofs << in+ in+ in+ "Color = BLACK" << endl;
    ofs << in+ in+ in+ "Text = Font_Monospace.render(TextLine, True, Color)" << endl;
    ofs << in+ in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ in+ "Text_Rect.topleft = (X, Y + Height * i)" << endl;
    ofs << in+ in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayWelcome(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render(self.MessageWelcome, True, BLACK, WHITE)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.center = self.Pos_Welcome" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayInput(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render(self.Message, True, BLACK)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.bottomleft = Screen.get_rect().bottomleft" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << in+ in+ "self.MessageTimer -= 1" << endl;
    ofs << endl;
    ofs << in+ "def SetMessage(self, Key):" << endl;
    ofs << in+ in+ "assert Key in self.Instructions" << endl;
    ofs << in+ in+ "return 'Input [' + Key + '] : ' + self.Instructions[Key] + '   '" << endl;
    ofs << endl;
    ofs << in+ "def DisplayScore(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render('Score: ' + str(round(self.Score)), True, RED)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.midtop = tuple(map(lambda i, j: i + j, Screen.get_rect().midtop, (0, Text.get_height())))" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayTime(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render('Simulation Time: ' + str(round(self.Time)), True, BLACK)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.midtop = Screen.get_rect().midtop" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "# TODO: Update display status" << endl;
    ofs << in+ "# def DisplayStatus(self, Ligand_Total, Ligand_Local, Ligand_Prev_Local, Am):" << endl;
    ofs << in+ "def DisplayStatus(self, Ligand_Local, Ligand_Prev_Local, Am_Local):" << endl;
    ofs << in+ in+ "Ligand_Now = Ligand_Local[0]" << endl;
    ofs << in+ in+ "Ligand_Prev = Ligand_Prev_Local[0]" << endl;
    ofs << in+ in+ "Am = Am_Local[-1, 0]" << endl;
    ofs << in+ in+ "dLigand = 0" << endl;
    ofs << in+ in+ "if Ligand_Now != 0:" << endl;
    ofs << in+ in+ in+ "dLigand = (Ligand_Now - Ligand_Prev) / Ligand_Now * 100" << endl;
    ofs << endl;
    ofs << in+ in+ "StatusText = 'Ligand @ RED :' + '{:.2f} '.format(Ligand_Now/ Unit / NA) + UnitTxt + '\\n' \\" << endl;
    ofs << in+ in+ in+ in+ in+ " + 'dLigand @ RED : ' + ('+' if dLigand >= 0 else '') + '{:.5f}'.format(dLigand) + ' % \\n' \\" << endl;
    ofs << in+ in+ in+ in+ in+ " + 'Am level of RED : ' + '{:.5f} '.format(Am / Unit / NA) + UnitTxt   # Get the last E coli's info " << endl;
    ofs << endl;
    ofs << in+ in+ "TextLines = StatusText.splitlines()" << endl;
    ofs << in+ in+ "Height = Font_Monospace.get_linesize() + 2" << endl;
    ofs << in+ in+ "X, Y = Screen.get_rect().topright" << endl;
    ofs << in+ in+ "Color = BLACK" << endl;
    ofs << in+ in+ "for i, TextLine in enumerate(TextLines):" << endl;
    ofs << in+ in+ in+ "if 'RED' in TextLine:" << endl;
    ofs << in+ in+ in+ in+ "Color = RED" << endl;
    ofs << in+ in+ in+ "Text = Font_Monospace.render(TextLine, True, Color)" << endl;
    ofs << in+ in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ in+ "Text_Rect.topright = (X, Y + Height * i)" << endl;
    ofs << in+ in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << in+ "def DisplayPause(self):" << endl;
    ofs << in+ in+ "Text = Font_Sans.render(self.MessagePause, True, BLACK, WHITE)" << endl;
    ofs << in+ in+ "Text_Rect = Text.get_rect()" << endl;
    ofs << in+ in+ "Text_Rect.center = self.Pos_Welcome" << endl;
    ofs << in+ in+ "Screen.blit(Text, Text_Rect)" << endl;
    ofs << endl;
    ofs << endl;
    ofs << "def main():" << endl;
    ofs << in+ "Control = FControl()" << endl;
    ofs << endl;
    ofs << in+ "SimUnitTime = 0.25" << endl;
    ofs << endl;
    ofs << in+ "PetriDish = FEnvironment()" << endl;
    ofs << endl;
    ofs << in+ "# TODO: Communicate to initialize in Sim" << endl;

    //                                           L,             qL
    std::vector<std::string> Color              {"Yellow",      "Blue",     "Green"     };
    std::vector<int> MaxBrightness              {200,           170,        50          };

    std::vector<std::string> Pattern            {"heatmap",     "heatmap",  "particles" };
    std::vector<std::string> NormalizationType  {"linear",      "log",      "particles" };
    std::vector<int> ReductionFactor            {5,             2,          1           };

    int i = 0;

    // Instantiate Molecules for Distribution
    // auto MolLoc = Context.GetSubList_LocationList("Molecule");
    for (auto Mol : MolLoc) {
        ofs << in+ Mol->Name << " = FMolecule('" << Mol->Name << "', ";
        ofs << Mol->Coord[0] << ", " << Mol->Coord[1] << ")" << endl;
        ofs << in+ Mol->Name << ".SetColor('" << Color[i] << "', " << MaxBrightness[i] << ")" << endl;
        ofs << in+ Mol->Name << ".SetPattern('" << Pattern[i] << "', '" << NormalizationType[i] << "', " << ReductionFactor[i] << ",)" << endl;
        ofs << endl;
        i++;
    }

    // Instantiate Organisms
    auto OrgNames = Context.GetUniqueNames_LocationList("Organism");
    for (auto OrganismName : OrgNames) {
        ofs << in+ OrganismName << " = FOrganism('" << OrganismName << "', " << "'Ecoli'" << ")" << endl; // TODO: Get Species later
        ofs << endl;
    }

    // Run Homeostasis for "HomeostasisMolName"
    for (auto OrganismName : OrgNames) {
        ofs << in+ OrganismName << ".Homeostasis(HomeostasisMolName)" << endl; // TODO: HARDCODED
    }
    ofs << endl;
    ofs << in+ "ElapsedTime = 0" << endl;
    ofs << in+ "PrevTime = datetime.now()" << endl;
    ofs << endl;
    ofs << in+ "SimState = True" << endl;
    ofs << in+ "while SimState:" << endl;
    ofs << endl;
    ofs << in+ in+ "if not Control.PauseSwitch:" << endl;
    ofs << in+ in+ in+ "CurrTime = datetime.now()" << endl;
    ofs << in+ in+ in+ "ElapsedTime += (CurrTime - PrevTime).total_seconds()" << endl;
    ofs << in+ in+ in+ "PrevTime = CurrTime" << endl;
    ofs << endl;
    ofs << in+ in+ "for event in pygame.event.get():" << endl;
    ofs << in+ in+ in+ "if event.type == pygame.QUIT:" << endl;
    ofs << in+ in+ in+ in+ "SimState = False" << endl;
    ofs << in+ in+ in+ "elif event.type == pygame.KEYDOWN:" << endl;
    ofs << in+ in+ in+ in+ "Control.MessageTimer = 5000" << endl;
    ofs << in+ in+ in+ in+ "if event.key == pygame.K_x:" << endl;
    ofs << in+ in+ in+ in+ in+ "SimState = False" << endl;
    ofs << endl;

    ofs << in+ in+ in+ in+ "# Organism Control" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_t:" << endl;
    for (auto OrganismName : OrgNames) {
        ofs << in + in + in + in + in + OrganismName << ".TrajectorySwitch = not " << OrganismName << ".TrajectorySwitch" << endl;
    }
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('T')" << endl;
    ofs << endl;
    ofs << in+ in+ in+ in+ "# Control panel" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_i:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.InstructionSwitch = not Control.InstructionSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('I')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_s:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.ScoreSwitch = not Control.ScoreSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('S')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_a:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.StatusSwitch = not Control.StatusSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('A')" << endl;
    ofs << in+ in+ in+ in+ "elif event.key == pygame.K_p:" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.PauseSwitch = not Control.PauseSwitch" << endl;
    ofs << in+ in+ in+ in+ in+ "Control.Message = Control.SetMessage('P')" << endl;
    ofs << endl;
    ofs << in+ in+ "Screen.fill(GRAY1)" << endl;
    ofs << endl;
    ofs << in+ in+ "# PetriDish.Draw()" << endl;
    ofs << in+ in+ "PetriDish.Draw(shape='circle')" << endl;
    ofs << endl;

    for (auto Mol : MolLoc) {
        ofs << in+ in+ Mol->Name << ".Draw(SimM.GetDistributionByName('" << Mol->Name << "'))" << endl;
    }
    ofs << endl;
    ofs << in+ in+ "if Control.PauseSwitch:" << endl;
    ofs << in+ in+ in+ "Control.DisplayPause()" << endl;
    ofs << endl;
    ofs << in+ in+ "else:" << endl;
    ofs << in+ in+ in+ "while ElapsedTime >= SimUnitTime:" << endl;
    ofs << endl;

    for (auto OrganismName : OrgNames) {
        ofs << in+ in+ in+ in+ OrganismName << ".Receptivity(50)" << endl;
        ofs << in+ in+ in+ in+ OrganismName << ".SetPosition(SimM.GetPositionXYAngleByName('" << OrganismName << "'))" << endl;
        ofs << in+ in+ in+ in+ OrganismName << ".ReportStatus()" << endl;
    }
    ofs << endl;

    ofs << in+ in+ in+ in+ "ElapsedTime -= SimUnitTime" << endl;
    ofs << in+ in+ in+ in+ "Control.Time += 1" << endl;
    ofs << endl;

    for (auto OrganismName : OrgNames) {
        ofs << in+ in+ OrganismName << ".AddToTrajectory()" << endl;
        ofs << in+ in+ "if " << OrganismName << ".TrajectorySwitch:" << endl;
        ofs << in+ in+ in+ OrganismName << ".DrawTrajectory()" << endl;
    }
    ofs << endl;

    ofs << endl;

    for (auto OrganismName : OrgNames) {
        ofs << in+ in+ OrganismName << ".Draw()" << endl;
    }
    ofs << in+ in+ "PetriDish.Draw(shape='lining')" << endl;
    ofs << endl;

    // Show first molecule only
    std::string MolName = MolLoc[0]->Name;
    std::string OrgName = OrgNames[0];

    ofs << in+ in+ MolName << "_Now = SimM.GetCountFromDistributionByNameAndPos('" << MolName << "', '" << OrgName << "')" << endl;

    ofs << endl;
    ofs << in+ in+ "if Control.Time < 50:" << endl;
    ofs << in+ in+ in+ "Control.DisplayWelcome()" << endl;
    ofs << in+ in+ "if Control.MessageTimer > 0:" << endl;
    ofs << in+ in+ in+ "Control.DisplayInput()" << endl;
    ofs << in+ in+ "if Control.InstructionSwitch:" << endl;
    ofs << in+ in+ in+ "Control.DisplayInstructions()" << endl;
    ofs << in+ in+ "if Control.ScoreSwitch:" << endl;
    ofs << in+ in+ in+ "Control.Score += " << MolName << "_Now / 10" << endl;
    ofs << in+ in+ in+ "Control.DisplayScore()" << endl;
    ofs << endl;
    ofs << in+ in+ "if Control.StatusSwitch:" << endl;
    ofs << in+ in+ in+ "Control.DisplayStatus(" << MolName << "_Now, " << OrgName << ".Ligand_Prev, " << OrgName << ".Am)" << endl;
    ofs << endl;
    ofs << in+ in+ "# Update Ligand Prev" << endl;
    ofs << in+ in+ OrgName << ".Ligand_Prev = " << MolName << "_Now" << endl;
    ofs << endl;
    ofs << in+ in+ "Control.DisplayTime()" << endl;
    ofs << in+ in+ "# Control.DisplayMoleculeGradient()" << endl;
    ofs << endl;
    ofs << in+ in+ "pygame.display.update()" << endl;
    ofs << endl;
    ofs << in+ "pygame.quit()" << endl;
    ofs << in+ "sys.exit()" << endl;
    ofs << endl;
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "main()" << endl;
    ofs << endl;

    std::cout << "  Visualization_2D program has been generated: ";
}
