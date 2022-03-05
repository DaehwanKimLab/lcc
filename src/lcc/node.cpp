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


