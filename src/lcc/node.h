#ifndef LCC_NODE_H
#define LCC_NODE_H
#include <iostream>
#include <vector>
#include <queue>
#include <memory>

//#include <llvm/IR/Value.h>


class NStatement;
class NMoleculeIdentifier;
class NIdentifier;
class NPathwayExpression;

typedef std::vector<std::shared_ptr<NStatement>> StatementList;
typedef std::vector<std::shared_ptr<NMoleculeIdentifier>> MoleculeList;
typedef std::vector<std::shared_ptr<NIdentifier>> IdentifierList;
typedef std::vector<std::shared_ptr<NPathwayExpression>> PathwayExprList;

class NNode;
class FTraversalContext {
public:
    std::queue<const NNode *> Queue;
    std::ostream& OutStream;

    FTraversalContext(std::ostream& InOutStream) : OutStream(InOutStream) {};
};

class NNode {
public:
	virtual ~NNode() {}

	virtual void Print(std::ostream& os) const {
		os << "Node";
	}

    virtual void Visit(FTraversalContext& Context) const {
        Context.OutStream << "Visit Node " << this << std::endl;
    }
};

class NExpression : public NNode {
    virtual void Visit(FTraversalContext& Context) const override;
};

class NStatement : public NNode {
    virtual void Visit(FTraversalContext& Context) const override;
};


class NIdentifier : public NExpression {
public:
    const std::string Name;

    NIdentifier() {}

    NIdentifier(const std::string &InName) : Name(InName) {}

    virtual void Print(std::ostream &os) const override {
        os << "NIdentifier(" << Name << ")";
    }

    virtual void Visit(FTraversalContext &Context) const override;
};

class NMoleculeIdentifier : public NIdentifier {
public:
    const int Id;

    NMoleculeIdentifier(int InId, const std::string& InName)
        : Id(InId), NIdentifier(InName) {}

    virtual void Print(std::ostream& os) const override {
        os << "NMoleculeIdentifier(" << Name << ", " << Id << ")";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NBlock : public NExpression {
public:
    StatementList Statements;

    NBlock() {}

    virtual void Print(std::ostream& os) const override {
        os << "NBlock(" << std::endl;
        for (const auto stmt : Statements) {
            os << "\t"; stmt->Print(os); os << std::endl;
        }
        os << ")" << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NReaction : public NExpression {
public:
    MoleculeList Reactants;
    MoleculeList Products;
    bool bBiDirection;

    NReaction() : bBiDirection(false) {}
    NReaction(const MoleculeList& InReactants, const MoleculeList& InProducts, bool bInBiDirection)
        : Reactants(InReactants), Products(InProducts), bBiDirection(bInBiDirection) {}

    virtual void Print(std::ostream& os) const override {
        os << "NReaction(" << std::endl;

        for (const auto& reactant : Reactants) {
            reactant->Print(os); os << ", ";
        }

        for (const auto& product : Products) {
            product->Print(os); os << ", ";
        }

        os << "BiDirection: " << bBiDirection;

        os << std::endl;

    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NGeneralReaction : public NExpression {
public:
    IdentifierList Reactants;
    IdentifierList Products;
    bool bBiDirection;

    NGeneralReaction() : bBiDirection(false) {}
    NGeneralReaction(const IdentifierList& InReactants, const IdentifierList& InProducts, bool bInBiDirection)
    : Reactants(InReactants), Products(InProducts), bBiDirection(bInBiDirection) {}

    virtual void Print(std::ostream& os) const override {
        os << "NGeneralReaction(" << std::endl;
        os << "  ";
        for (const auto& item: Reactants) {
             os << item->Name << ", ";
        }
        os << std::endl;
        os << "  ";
        for (const auto& item: Products) {
            os << item->Name << ", ";
        }
        os << std::endl;
        os << "), " << "BiDirection: " << bBiDirection << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};


class NProteinDeclaration : public NStatement {
public:
    const NIdentifier Id;
    NReaction Reaction;
    std::shared_ptr<NBlock> Block;

    NProteinDeclaration(const NIdentifier& InId, const NReaction& InReaction, NBlock* InBlock)
        : Id(InId), Reaction(InReaction), Block(InBlock) {}

    NProteinDeclaration(const NIdentifier &InId, const NReaction& InReaction)
            : Id(InId), Reaction(InReaction) {}

    NProteinDeclaration(const NIdentifier &InId)
            : Id(InId) {}

    virtual void Print(std::ostream &os) const override {
        os << "NProteinDeclaration: ";
        Id.Print(os);
        os << ", ";
        Reaction.Print(os);
        os << std::endl;
    }

    virtual void Visit(FTraversalContext &Context) const override;
};

class NPathwayExpression : public NExpression {
public:
    int Type;
    std::shared_ptr<NExpression> Lhs;
    std::shared_ptr<NExpression> Rhs;

    NPathwayExpression(NPathwayExpression* InLhs, NPathwayExpression* InRhs, int InType)
        : Lhs(InLhs), Rhs(InRhs), Type(InType) {}

    virtual void Print(std::ostream& os) const override {
        os << "NPathwayExpression(Type: " << Type << ", ";
        Lhs->Print(os); os << ", ";
        Rhs->Print(os); os << ")";
        os << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPathwayDeclaration : public NStatement {
public:
    const NIdentifier Id;

    NPathwayExpression* PathwayExpression;
    NBlock* Block;

    NPathwayDeclaration(const NIdentifier& InId, NPathwayExpression* InPathwayExpression)
        : Id(InId), PathwayExpression(InPathwayExpression), Block(nullptr) {}

    /* New */
    NPathwayDeclaration(const NIdentifier& InId, NBlock* InBlock)
        : Id(InId), Block(InBlock), PathwayExpression(nullptr) {}

    virtual ~NPathwayDeclaration() {
        if (Block) delete Block;
        if (PathwayExpression) delete PathwayExpression;
    }


    virtual void Print(std::ostream& os) const override {
        os << "NPathwayDeclaration("; Id.Print(os); os << std::endl;
        if (PathwayExpression) PathwayExpression->Print(os);
        if (Block) Block->Print(os);
        os << ")" << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NDescriptionStatement : public NStatement {
public:
    const std::string Description;

    NDescriptionStatement(const std::string& InDescription)
        : Description(InDescription) {}

    virtual void Print(std::ostream& os) const override {
        os << "Description: " << Description;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPathwayReactionIDStatement : public NStatement {
public:
    const NIdentifier Id;

    NPathwayReactionIDStatement(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "NPathwayReactionIDStatement("; Id.Print(os); os << ")" << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPathwayReactionStatement : public NStatement {
public:
    NPathwayExpression PathwayExpression;

    NPathwayReactionStatement(NPathwayExpression& InPathwayExpression)
        : PathwayExpression(InPathwayExpression) {}


    virtual void Print(std::ostream& os) const override {
        os << "NPathwayExpression("; PathwayExpression.Print(os); os << ")" << std::endl;
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NOrganismDeclaration : public NStatement {
public:
    const NIdentifier Id;
    const std::string Description;

    NOrganismDeclaration(const NIdentifier& InId, const std::string& InDescription)
        : Id(InId), Description(InDescription) {}
    NOrganismDeclaration(const NIdentifier& InId)
        : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "NOrganismDeclaration(";
        Id.Print(os); os << ", " << Description << ")" << std::endl;
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NExperimentDeclaration : public NStatement {
public:
    const NIdentifier Id;

    const std::string Description;
    std::shared_ptr<NBlock> Block;

    NExperimentDeclaration(const NIdentifier& InId, const std::string& InDescription)
        : Id(InId), Description(InDescription), Block(nullptr) {}
    NExperimentDeclaration(const NIdentifier& InId)
        : Id(InId), Block(nullptr) {}
    /* New */
    NExperimentDeclaration(const NIdentifier& InId, NBlock* InBlock)
        : Id(InId), Block(InBlock) {}

    virtual void Print(std::ostream& os) const override {
        os << "NExperimentDeclaration(";
        Id.Print(os); os << ", " << Description << ")" << std::endl;
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NPropertyStatement : public NStatement {
public:
    const std::string Key;
    const std::string Value;

    NPropertyStatement(const std::string& InKey, const std::string& InValue)
        : Key(InKey), Value(InValue) {}


    virtual void Print(std::ostream& os) const override {
        os << "Property(" << Key << "): " << Value;
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NProteinCofactorStatement : public NStatement {
public:
    const NIdentifier Id;
    IdentifierList CofactorList;

    NProteinCofactorStatement(const NIdentifier& InId, const IdentifierList& InCofactorList) : Id(InId), CofactorList(InCofactorList) {};

    virtual void Print(std::ostream& os) const override {
        os << "Cofactor(" << Id.Name << "): ";
        for (const auto& arg: CofactorList) {
            os << arg->Name << ", ";
        }
        os << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NProteinDomainStatement : public NStatement {
public:
    const NIdentifier Id;
    IdentifierList DomainList;

    NProteinDomainStatement(const NIdentifier& InId, const IdentifierList& InDomainList) : Id(InId), DomainList(InDomainList) {};

    virtual void Print(std::ostream& os) const override {
        os << "Domain(" << Id.Name << "): ";
        for (const auto& domain: DomainList) {
            os << domain->Name << ", ";
        }
        os << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NProteinStepStatement : public NStatement {
public:
    const NIdentifier Id;
    NGeneralReaction Reaction;

    NProteinStepStatement(const NIdentifier& InId, const NGeneralReaction& InReaction)
    : Id(InId), Reaction(InReaction) {}

    virtual void Print(std::ostream& os) const override {
        os << "NProteinStepStatement(" << Id.Name << "): ";
        Reaction.Print(os);
    }

    virtual void Visit(FTraversalContext &Context) const override {};
};

class NProteinSequenceStatement : public NStatement {
public:
    const NIdentifier Id;
    IdentifierList StepList;

    NProteinSequenceStatement(const NIdentifier& InId, const IdentifierList& InStepList)
    : Id(InId), StepList(InStepList) {}

    virtual void Print(std::ostream& os) const override {
        os << "Protein Process Sequence(" << Id.Name << "): ";
        for (const auto& step: StepList) {
            os << step->Name << ", ";
        }
        os << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NUsingStatement : public NStatement {
public:
    int Type;
    const NIdentifier Id;

    NUsingStatement(int InType, const NIdentifier& InId) : Type(InType), Id(InId) {};

    virtual void Print(std::ostream& os) const override {
        os << "Using(" << Type << "): "; Id.Print(os); os << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NDummyDeclaration : public NStatement {
public:
    const std::string StringLiteral;

    NDummyDeclaration(const std::string& InStringLiteral) : StringLiteral(InStringLiteral) {}

    virtual void Print(std::ostream& os) const override {
        os << "NDummyDeclaration(" << StringLiteral << ")" << std::endl;
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

#endif /* LCC_NODE_H */
