#include <iostream>
#include "node.h"

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

void NProteinDeclaration::Visit(FTraversalContext &Context) const {
//    Context.OutStream << "Protein Declaration(" << Id.Name << ")" << std::endl;

    Context.Queue.push(&OverallReaction);
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
