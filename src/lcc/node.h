#include <iostream>
#include <vector>

//#include <llvm/IR/Value.h>


class NStatement;
class NMoleculeIdentifer;
class NIdentifier;
class NPathwayExpression;

typedef std::vector<NStatement*> StatementList;
typedef std::vector<NMoleculeIdentifer*> MoleculeList;
typedef std::vector<NIdentifier*> IdentifierList;
typedef std::vector<NPathwayExpression*> PathwayExprList;


class NNode {
public:
	virtual ~NNode() {}

	virtual void Print(std::ostream& os) const {
		os << "Node";
	}
};

class NExpression : public NNode {
};

class NStatement : public NNode {
};


class NIdentifier : public NExpression {
public:
    std::string Name;
    NIdentifier() {}
    NIdentifier(const std::string& InName) : Name(InName) {}

    virtual void Print(std::ostream& os) const override {
        os << "NIdentifier(" << Name << ")";
    }
};

class NMoleculeIdentifer : public NExpression {
public:
    const int Id;
    const std::string Name;

    NMoleculeIdentifer(int InId, const std::string& InName)
        : Id(InId), Name(InName) {}

    virtual void Print(std::ostream& os) const override {
        os << "NMoleculeIdentifer(" << Name << ", " << Id << ")";
    }

};

class NBlock : public NExpression {
public:
    StatementList Statements;

    NBlock() {}
    virtual ~NBlock() {
        for (auto& stmt : Statements) {
            delete stmt;
        }
    }

    virtual void Print(std::ostream& os) const override {
        os << "NBlock(" << std::endl;
        for (const auto stmt : Statements) {
            os << "\t"; stmt->Print(os); os << std::endl;
        }
        os << ")" << std::endl;
    }

};

class NReaction : public NExpression {
public:
    MoleculeList Reactants;
    MoleculeList Products;
    bool bBiDirection;

    NReaction() : bBiDirection(false) {}
    NReaction(const MoleculeList& InReactants, const MoleculeList& InProducts, bool bInBiDirection)
        : Reactants(InReactants), Products(InProducts), bBiDirection(bInBiDirection) {}

    virtual ~NReaction() {
        for(auto& mol : Reactants) {
            delete mol;
        }
        for(auto& mol : Products) {
            delete mol;
        }
    }

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

};



class NProteinDeclaration : public NStatement {
public:
    const NIdentifier& Id;
    NReaction Reaction;

    NProteinDeclaration(const NIdentifier& InId, const NReaction& InReaction) 
        : Id(InId), Reaction(InReaction) {}
    NProteinDeclaration(const NIdentifier& InId)
        : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "NProteinDeclaration: " ; Id.Print(os); os << ", ";
        Reaction.Print(os);
        os << std::endl;
    }
};

class NPathwayExpression : public NExpression {
public:
    int Type;
    NExpression& Lhs;
    NExpression& Rhs;

    NPathwayExpression(NPathwayExpression& InLhs, NPathwayExpression& InRhs, int InType) 
        : Lhs(InLhs), Rhs(InRhs), Type(InType) {}

    virtual void Print(std::ostream& os) const override {
        os << "NPathwayExpression(Type: " << Type << ", ";
        Lhs.Print(os); os << ", ";
        Rhs.Print(os); os << ")";
        os << std::endl;
    }
};

class NPathwayDeclaration : public NStatement {
public:
    const NIdentifier& Id;

    NPathwayExpression& PathwayExpression;

    NPathwayDeclaration(const NIdentifier& InId, NPathwayExpression& InPathwayExpression) 
        : Id(InId), PathwayExpression(InPathwayExpression) {}


    virtual void Print(std::ostream& os) const override {
        os << "NPathwayDeclaration("; Id.Print(os); os << std::endl;
        PathwayExpression.Print(os);
        os << ")" << std::endl;
    }
};

