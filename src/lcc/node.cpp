#include <iostream>
#include "node.h"
#include "number.h"

using namespace std;

void NExpression::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "NExpression" << std::endl;
}
void NStatement::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "NStatement" << std::endl;
}

void NIdentifier::Visit(FTraversalContext& Context) const {
//    std::cerr << "NIdentifier(" << Name << ")" << std::endl;
}

void NMoleculeIdentifier::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "NMoleculeIdentifier(" << Name << ", " << Id << ")"<< std::endl;
}

void NSubstrate::Visit(FTraversalContext& Context) const {
}

void NBlock::Visit(FTraversalContext& Context) const {
    for (auto& stmt: Statements) {
        Context.Queue.push(static_cast<const NNode *>(stmt.get())); // Fixme
    }
}

void NMoleculeReaction::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Reaction(" << bBiDirection << ")" << std::endl;

    for (auto& reactant: Reactants) {
        Context.Queue.push(static_cast<const NNode *>(reactant.get())); // Fixme
    }

    for (auto& product: Products) {
        Context.Queue.push(static_cast<const NNode *>(product.get())); // Fixme
    }
}

void NReaction::Visit(FTraversalContext &Context) const {

}

void NReactionDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Reaction Declaration(" << Id.Name << ")" << std::endl;

    if (OverallReaction) {
        Context.Queue.push(OverallReaction);
    }
}

void NProteinDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Protein Declaration(" << Id.Name << ")" << std::endl;

    if (OverallReaction) {
        Context.Queue.push(OverallReaction);
    }
}

void NProcessDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Process Declaration(" << Id.Name << ")" << std::endl;

    Context.Queue.push(Block.get());
}



void NPathwayExpression::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "PathwayExpression(" << Type << ")" << std::endl;

    Context.Queue.push(Lhs.get());
    Context.Queue.push(Rhs.get());

}

void NPathwayDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "PathwayDeclaration(" << Id.Name << ")" << std::endl;

    if (PathwayExpression) Context.Queue.push(PathwayExpression.get());
    if (Block) Context.Queue.push(Block.get());
}


void NDescriptionStatement::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Description(" << Description << ")" << std::endl;
}

void NPathwayReactionIDStatement::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "PathwayReactionID(" << Id.Name << ")" << std::endl;
}

void NPathwayReactionStatement::Visit(FTraversalContext &Context) const {
    Context.Queue.push(static_cast<const NNode *>(&PathwayExpression));
}

void NOrganismDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Organism(" << Id.Name << ", " << Description << ")" << std::endl;
}

void NExperimentDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Experiment(" << Id.Name << ", " << Description << ")" << std::endl;
    if (Block) Context.Queue.push(Block.get());
}

void NPropertyStatement::Visit(FTraversalContext& Context) const {
}

void NProteinCofactorStatement::Visit(FTraversalContext& Context) const {
}

void NProteinDomainStatement::Visit(FTraversalContext &Context) const {
}

void NUsingStatement::Visit(FTraversalContext& Context) const {
}

void NPolymeraseDeclaration::Visit(FTraversalContext& Context) const {
    for (auto& stmt: Statements) {
        Context.Queue.push(static_cast<const NNode *>(stmt.get())); // Fixme
    }
}

void NRibosomeDeclaration::Visit(FTraversalContext& Context) const {
    for (auto& stmt: Statements) {
        Context.Queue.push(static_cast<const NNode *>(stmt.get())); // Fixme
    }
}

void NDummyDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Dummy(" << StringLiteral << ")" << std::endl;
}

void NIfStatement::Visit(FTraversalContext& Context) const {
    if (CondExpression) Context.Queue.push(CondExpression.get());
    if (Body) Context.Queue.push(Body.get());
    if (ElseBody) Context.Queue.push(ElseBody.get());
}

void NLoopStatement::Visit(FTraversalContext& Context) const {
    if (Body) Context.Queue.push(Body.get());
}

void NExpressionStatement::Visit(FTraversalContext& Context) const {
    if (Expression) Context.Queue.push(Expression.get());

}

// parsing routines
//std::vector<std::pair<std::string, int>> NReaction::GetStoichFromReaction(bool bProductIsReaction=false) {
//
//    map<string, int> Stoichiometry;
//		string Location = Reaction->Location.Name;
//
//    std::string Name;
//    int Coefficient;
//    for (const auto& reactant : Reaction->Reactants) {
//        Name = reactant->Id.Name;
//        Coefficient = -reactant->Coeff; // update when coeff is fully implemented in parser
//        std::cout << "    Reactant: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
//        Stoichiometry[Name] = Coefficient;
//
//    }
//
//    for (const auto& product : Reaction->Products) {
//        Name = product->Id.Name;
//        Coefficient = product->Coeff; // update when coeff is fully implemented in parser
//        std::cout << "    Product: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
//        if (Stoichiometry.count(Name) > 0) {
//            Coefficient += Stoichiometry[Name];
//            std::cout << "    Updated Stoichiometry: " << "(" << Coefficient << ") " << Name << ", " << std::endl;
//        }
//        Stoichiometry[Name] = Coefficient;
//    }
//
//		if (!Location.empty()) {
//        std::cout << "    Location: " << Location << endl;
//		}
//
//    // convert stoichiometry map to vector of pairs and add new molecules to the system
//    std::vector<std::pair<std::string, int>> Stoichiometry_Ordered;
//
//    for (auto& stoich : Stoichiometry) {
//        std::pair<std::string, int> Stoich(stoich.first, stoich.second);
//        Stoichiometry_Ordered.push_back(Stoich);
//
//        // Do not add new molecule if the product is a reaction name
//        if ((stoich.second > 0) & (bProductIsReaction)) {
//            continue;
//        }
//        FMolecule * NewMolecule = new FMolecule(stoich.first);
//        if (Option.bDebug) { NewMolecule->Print(os); }
//        Context.AddToMoleculeList(NewMolecule);
//    } 
//
//    return Stoichiometry_Ordered;
//}
//            
//void NReaction::AddEnzReaction(std::string ReactionName, std::string EnzymeName)
//{
//    float k = Float_Init;
//    float krev = Float_Init;
//
//    // properties
//    const auto& propertylist = Reaction->Property;
//    for (auto& property :propertylist) {
//        auto& Key = property->Key;
//        // auto Value = property->Value->Evaluate();
//        const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
//        auto Value = Value_Exp->EvaluateValueAndPrefix();
//
//        if (Key == "k") {
//            k = std::stof(Value);
//        } else if (Key == "krev") {
//            krev = std::stof(Value);
//        }
//    }
//
//    std::vector<std::pair<std::string, int>> Stoichiometry = GetStoichFromReaction(Reaction);
//    string Location = Reaction->Location.Name;
//
//    // for Enz_Standard Reaction
//    // Fill in presumably irreversible reaction kinetic values 
//    if      ((k != Float_Init) & (krev == Float_Init))   { krev = 0; }
//    else if ((k == Float_Init) & (krev != Float_Init))   { k = 0; }
//
//    // add new enzymatic reaction to the system
//    if ((k >= 0) & (krev >= 0)) {
//        FEnz_StandardReaction *NewReaction = new FEnz_StandardReaction(ReactionName, Stoichiometry, EnzymeName, k, krev);
//        if (Option.bDebug) { NewReaction->Print(os); }
//        Context.AddToReactionList(NewReaction);
//
//    } else {
//        FEnzymaticReaction *NewReaction = new FEnzymaticReaction(ReactionName, Stoichiometry, EnzymeName);
//        if (Option.bDebug) { NewReaction->Print(os); }
//        Context.AddToReactionList(NewReaction);
//    }
//}
//
//std::pair<std::string, std::vector<float>> NReaction::GetEnzKinetics(std::string EnzymeName)
//{
//    std::pair<std::string, std::vector<float>> SubConstPair;
//    std::string Substrate;
//    float kcat = Float_Init;
//    float KM = Float_Init;                
//
//    // properties
//    const auto& propertylist = Reaction->Property;
//    for (auto& property :propertylist) {
//        auto& Key = property->Key;
//
//        if (Key == "Substrate") {
//            auto Value = property->Value->Evaluate();
//            Substrate = std::stof(Value);
//
//        } else if ((Key == "kcat") || (Key == "kCat")) {
//            const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
//            auto Value = Value_Exp->EvaluateValueAndPrefix();
//            kcat = std::stof(Value);
//
//        } else if ((Key == "KM") || (Key == "kM") || (Key == "km")) {
//            const auto Value_Exp = dynamic_pointer_cast<const NConstantExpression>(property->Value);
//            auto Value = Value_Exp->EvaluateValueAndPrefix();
//            KM = std::stof(Value);
//        }
//    }
//     
//
//    // if Substrate not defined by user input, search the database            
//    if (Substrate.empty()) {
//        Substrate = Context.QueryTable(EnzymeName, "Substrate", Context.EnzymeTable);
//        if (!Substrate.empty()) {
//            os << "  Substrate imported from database: " << Substrate << endl;
//        }
//    }
//
//    // import constants from database
//    if (kcat == Float_Init) {
//        string kcat_Database = Context.QueryTable(EnzymeName, "kcat", Context.EnzymeTable);
//        if (!kcat_Database.empty()) {
//            kcat = std::stof(kcat_Database);
//            os << "  kcat imported from database: " << kcat_Database << endl;
//        }
//    }
//    if (KM == Float_Init) {
//        string KM_Database = Context.QueryTable(EnzymeName, "KM", Context.EnzymeTable);
//        if (!KM_Database.empty()) {
//            KM = std::stof(KM_Database);
//            os << "  KM imported from database: " << KM_Database << endl;
//        }
//    }
//
//    // Use first reactant as substrate for MichaelisMenten if still not assigned
//    // This may not always work with Michaelis Menten without database. Excluding common molecules will improve a chance.
//    if (Substrate.empty()) {
//        for (const auto& reactant : Reaction->Reactants) {
//            if (reactant->Id.Name == EnzymeName) {
//                continue;
//            }
//            Substrate = reactant->Id.Name;
//            break;
//        }
//    }
//
//    Utils::Assertion(((kcat >= 0) & (KM >= 0)), "Kinetic constants are not properly extracted: " + EnzymeName);
// 
//    std::vector<float> KConstants{ kcat, KM };
//    std::cout << EnzymeName << ", " << Substrate << ", " << kcat << ", " << KM << std::endl;
//    SubConstPair = {Substrate, KConstants};
//    
//    return SubConstPair;
//}

std::string NVariableExpression::GetName() const {
    std::string Name_Out;

    if (Utils::is_class_of<const NIdentifier, const NExpression>(Variable.get())) {
        const auto Id = dynamic_pointer_cast<const NIdentifier>(Variable);
        Name_Out = Id->Name;
    } else if (Utils::is_class_of<const NVariableExpression, const NExpression>(Variable.get())) {
        const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(Variable);
        Name_Out = VarExp->GetName();
    } else if (Utils::is_class_of<NFunctionCallExpression, NExpression>(Variable.get())) {
        const auto FCExp = dynamic_pointer_cast<const NFunctionCallExpression>(Variable);
        Name_Out = FCExp->GetName();
    }

    return Name_Out;
}

std::vector<float> NRangeExpression::GetBeginEndStep() const {
    float fBegin, fEnd, fStep = Numbers::GetFloatDefault();

    if (Begin) {
        if (Utils::is_class_of<NConstantExpression, NExpression>(Begin.get())) {
            fBegin = std::stof(dynamic_pointer_cast<const NConstantExpression>(Begin)->Evaluate());
            Utils::Assertion(fBegin >= 0, "Range input (beginning) cannot be negative.");
        }
    }
    if (End) {
        if (Utils::is_class_of<NConstantExpression, NExpression>(End.get())) {
            fEnd = std::stof(dynamic_pointer_cast<const NConstantExpression>(End)->Evaluate());
            Utils::Assertion(fEnd > 0, "Range input (end) cannot be 0 or negative.");
        }
    }
    if (Step) {
        if (Utils::is_class_of<NConstantExpression, NExpression>(Step.get())) {
            fStep = std::stof(dynamic_pointer_cast<const NConstantExpression>(Step)->Evaluate());
            Utils::Assertion(fStep > 0, "Range input (step) cannot not be 0 or negative.");
        }
    }

    if (fStep < 0) {
        fStep  = 0;   // default value
    }

    // categorizing
    if ((fBegin <  0) & (fEnd < 0))      { fBegin = 0; fEnd = -1;} //    os << "Begin<0 & End<0: " << Name << endl;} // fixed amount
    else if ((fBegin >= 0) & (fEnd < 0)) {             fEnd = fBegin;} //    os << "Begin>=0 & End<0: " << Name << endl; } // single step event treated the same as range for now

    // package into a vector and return
    std::vector<float> Range = {fBegin, fEnd, fStep};
 
    return Range;
}

std::string NFunctionCallExpression::GetName() const {
    std::string Name_Out;

    if (Utils::is_class_of<const NVariableExpression, const NExpression>(Name.get())){
        const auto VarExp = dynamic_pointer_cast<const NVariableExpression>(Name);
        Name_Out = VarExp->GetName();
    }
    return Name_Out;
}

std::vector<float> NFunctionCallExpression::GetParameters(std::string Type = "") const {
    std::vector<float> Parameters;
    bool bLocation = false;
    if (Type == "Location") {
        bLocation = true;
    }     

    for (const auto& arg: * Args) {

        if (Utils::is_class_of<NConstantExpression, NExpression>(arg.get())) {
            const auto& parameter = dynamic_cast<const NConstantExpression *>(arg.get());
            Parameters.push_back(std::stof(parameter->Evaluate()));
        // Other NExpression classes may be added for parsing
//        } else if () {
        }
    }

    if (bLocation) {
        Utils::Assertion((Parameters.size() == 3), "Location information requires three (X, Y, Z) coordinates");
    }

    return Parameters;
}


