#ifndef LCC_NODE_H
#define LCC_NODE_H
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <memory>
#include <cassert>

//#include <llvm/IR/Value.h>


class NStatement;
class NMoleculeIdentifier;
class NIdentifier;
class NReactionDeclaration;
class NPathwayExpression;
class NProteinDeclaration;
class NPolymeraseDeclaration;
class NPropertyStatement;
class NSubstrate;
class NExpression;
class NDeclaraionStatement;

typedef std::vector<std::shared_ptr<NStatement>> StatementList;
typedef std::vector<std::shared_ptr<NMoleculeIdentifier>> MoleculeList;
typedef std::vector<std::shared_ptr<NIdentifier>> IdentifierList;
typedef std::vector<std::shared_ptr<NReactionDeclaration>> ReactionDeclList;
typedef std::vector<std::shared_ptr<NPathwayExpression>> PathwayExprList;
typedef std::vector<std::shared_ptr<NProteinDeclaration>> ProteinDeclList;
typedef std::vector<std::shared_ptr<NPolymeraseDeclaration>> PolymeraseDeclList;
typedef std::vector<std::shared_ptr<NPropertyStatement>> PropertyList;
typedef std::vector<std::shared_ptr<NSubstrate>> SubstrateList;
typedef std::vector<std::shared_ptr<NExpression>> ExpressionList;

class NNodeUtil {
public:
    static StatementList* InitStatementList(NStatement* InStatementPtr = nullptr)
    {
        StatementList* ptr = new StatementList();
		if (InStatementPtr) {
			ptr->emplace_back(InStatementPtr);
		}
        return ptr;
    }
};


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
	}

    virtual void Visit(FTraversalContext& Context) const {
        Context.OutStream << "Visit Node " << this << std::endl;
    }
};

class NExpression : public NNode {
public:
    NExpression() {};
    virtual void Visit(FTraversalContext& Context) const override;

    virtual std::string Evaluate() const {
        /* Temporary */
        return "";
    }
};

class NStatement : public NNode {
public:
	NStatement() {};

    virtual void Visit(FTraversalContext& Context) const override;
};


class NIdentifier : public NExpression {
public:
    std::string Name;

    NIdentifier() {}

    NIdentifier(const std::string &InName) : Name(InName) {}

    virtual void Print(std::ostream &os) const override {
        os << "Id: " << Name;
    }

    virtual std::string Evaluate() const override {
        return Name;
    }

    virtual void Visit(FTraversalContext &Context) const override;
};

class NMoleculeIdentifier : public NIdentifier {
public:
    const int Id;

    NMoleculeIdentifier(int InId, const std::string& InName)
        : Id(InId), NIdentifier(InName) {}

    virtual void Print(std::ostream& os) const override {
        os << "MolIdentifier: {";
        os << Name << ": " << Id;
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NSubstrate : public NExpression {
public:
    int Coeff;
    NIdentifier Id;

    NSubstrate(NIdentifier& InId)
        : Id(InId), Coeff(1) {};

    NSubstrate(NIdentifier& InId, int InCoeff)
        : Id(InId), Coeff(InCoeff) {};

    NSubstrate(NIdentifier& InId, const std::string& InCoeff)
        : Id(InId)
    {
        Coeff = std::stoi(InCoeff);
    }

    virtual void Print(std::ostream& os) const override {
        os << "Substrate: {";
        os << Id.Name << ": " << Coeff;
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NBlock : public NExpression {
public:
    StatementList Statements;

    NBlock() {}


    void AddStatment(NStatement* InStatement) {
        Statements.emplace_back(InStatement);
    }

    void AddStatment(const std::vector<std::shared_ptr<NStatement>>* InStatements) {
        Statements.insert(Statements.end(), InStatements->begin(), InStatements->end());
    }

    virtual void Print(std::ostream& os) const override {
        os << "Block: {";
        for (const auto stmt : Statements) {
            stmt->Print(os); os << ", ";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NMoleculeReaction : public NExpression {
public:
    MoleculeList Reactants;
    MoleculeList Products;
    bool bBiDirection;

    NMoleculeReaction() : bBiDirection(false) {}
    NMoleculeReaction(const MoleculeList& InReactants, const MoleculeList& InProducts, bool bInBiDirection)
        : Reactants(InReactants), Products(InProducts), bBiDirection(bInBiDirection) {}

    virtual void Print(std::ostream& os) const override {
        os << "MoleculeReaction: {";

        os << "Reactants: [";
        for (const auto& reactant : Reactants) {
            reactant->Print(os); os << ", ";
        }

        os << "], ";
        os << "Products: [";
        for (const auto& product : Products) {
            product->Print(os); os << ", ";
        }
        os << "], ";
        os << "BiDirection: " << bBiDirection;
        os << ", ";
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPropertyStatement : public NStatement {
public:
    const std::string Key;
    //const std::string Value;
    std::shared_ptr<NExpression> Value;

    NPropertyStatement() {}
    NPropertyStatement(const std::string& InKey)
        : Key(InKey) {}
    NPropertyStatement(const std::string& InKey, NExpression* InValue)
        : Key(InKey), Value(InValue) {}

    virtual void Print(std::ostream& os) const override {
        os << "Property: {";
        os << Key;
        if (Value) {
            os << ": ";
            Value->Print(os);
        }
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override;

};


class NReaction : public NStatement {
public:
    NIdentifier Id;
    //IdentifierList Reactants;
    //IdentifierList Products;

    SubstrateList Reactants;
    SubstrateList Products;

    bool bEffect; // positive: true, negative: false
//    bool bBiDirection;
	NIdentifier Location;
    PropertyList Property;

//    NReaction() : bBiDirection(false) {}
//    NReaction(const SubstrateList& InReactants, const SubstrateList& InProducts, bool bInBiDirection)
//    : Reactants(InReactants), Products(InProducts), bBiDirection(bInBiDirection) {}

    NReaction() {}
    NReaction(const SubstrateList& InReactants, const SubstrateList& InProducts, bool bInEffect)
    : Reactants(InReactants), Products(InProducts), bEffect(bInEffect) {}

    void SetID(const NIdentifier& InId) {
        Id = InId;
    }

	void SetLocation(const NIdentifier& InLocation) {
		Location = InLocation;
	}

    void SetProperty(const PropertyList& InProperty) {
        Property = InProperty;
    }

    void AddProperty(NPropertyStatement* InProperty) {
        Property.emplace_back(InProperty);
    }

    void AddProperty(const std::string& InName, NExpression* InValue) {
        Property.emplace_back(new NPropertyStatement(InName, InValue));
    }

    virtual void Print(std::ostream& os) const override {
        os << "Reaction: {";
        os << "Reactants: [";
        for (const auto& item: Reactants) {
             item->Print(os); os << ", ";
        }
        os << "], ";
        os << "Products: [";
        for (const auto& item: Products) {
            item->Print(os); os << ", ";
        }
        os << "], ";
        os << "Effect: ";
        if (bEffect)       { os << "positive" ; } 
        else if (!bEffect) { os << "negative" ; }
        os << ", ";

        for (const auto& property: Property) {
            property->Print(os); os << ", ";
        }
        os << "}" << std::endl;
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NReactionDeclaration : public NStatement {
public:
    const NIdentifier Id;
    const NReaction* OverallReaction;
//    std::shared_ptr<NBlock> Block;

//    NReactionDeclaration(const NIdentifier& InId, const NReaction& InOverallReaction, NBlock* InBlock)
//        : Id(InId), OverallReaction(InOverallReaction), Block(InBlock) {}

    NReactionDeclaration(const NIdentifier& InId, const NReaction* InOverallReaction)
            : Id(InId), OverallReaction(InOverallReaction) {}

    NReactionDeclaration(const NIdentifier& InId)
            : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "Reaction: {";
        Id.Print(os); os << ", ";
        if (OverallReaction) {
            OverallReaction->Print(os); os << ", ";
        }
//        if (Block) {
//            for (const auto& stmt: Block->Statements) {
//                stmt->Print(os); os << ", ";
//            }
//        }
        os << "}";
    }

    virtual void Visit(FTraversalContext &Context) const override;
};

class NProteinDeclaration : public NStatement {
public:
    const NIdentifier Id;
    const NReaction* OverallReaction;
    std::shared_ptr<NBlock> Block;

    NProteinDeclaration(const NIdentifier& InId, const NReaction* InOverallReaction, NBlock* InBlock)
        : Id(InId), OverallReaction(InOverallReaction), Block(InBlock) {}

    NProteinDeclaration(const NIdentifier& InId, const NReaction* InOverallReaction)
            : Id(InId), OverallReaction(InOverallReaction) {}

    NProteinDeclaration(const NIdentifier& InId, NBlock* InBlock)
        : Id(InId), Block(InBlock), OverallReaction(nullptr) {}

    NProteinDeclaration(const NIdentifier& InId)
            : Id(InId), OverallReaction(nullptr) {}

    virtual void Print(std::ostream& os) const override {
        os << "Protein: {";
        Id.Print(os); os << ", ";
        if (OverallReaction) {
            OverallReaction->Print(os); os << ", " << std::endl;
        }
        if (Block) {
            for (const auto& stmt: Block->Statements) {
                stmt->Print(os); os << ", " << std::endl;
            }
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext &Context) const override;
};

class NProteinComplexDeclaration : public NStatement {
public:
    const NIdentifier Id;
    IdentifierList Components;

    NProteinComplexDeclaration(const NIdentifier& InId, const IdentifierList& InComponents)
    : Id(InId), Components(InComponents) {}

    virtual void Print(std::ostream& os) const override {
        os << "ProteinComplex: {";
        os << Id.Name << ": [";
        for(const auto& item: Components) {
            os << item->Name << ", ";
        }
        os << "]";
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override {};
};

class NProcessDeclaration : public NStatement {
public:
    const NIdentifier Id;
    const NReaction* OverallReaction;
    std::shared_ptr<NBlock> Block;

    NProcessDeclaration(const NIdentifier& InId, const NReaction* InOverallReaction, NBlock* InBlock)
        : Id(InId), OverallReaction(InOverallReaction), Block(InBlock) {}

    NProcessDeclaration(const NIdentifier& InId, const NReaction* InOverallReaction)
        : Id(InId), OverallReaction(InOverallReaction) {}

    NProcessDeclaration(const NIdentifier& InId)
        : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "Process: {";
        Id.Print(os);
        os << ", ";
        if (OverallReaction) {
            OverallReaction->Print(os);
        }
        if (Block) {
            os << ", ";
            Block->Print(os);
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;

};

class NChainReactionExpression : public NExpression {
public:
    std::vector<NIdentifier> Identifiers;
    std::vector<int> Operations;

    void Add(const NIdentifier& InIdentifier, int Operation = 0) {
        Identifiers.emplace_back(InIdentifier);
        Operations.emplace_back(Operation);

        assert(Identifiers.size() == Operations.size());
    }

    virtual void Print(std::ostream& os) const override {
        assert(Identifiers.size() == Operations.size());
        os << "ChainReactionExpr: {";
        os << "Components: [";
        for(int i = 0; i < Operations.size(); i++) {
            os << Identifiers[i].Name << ", ";
        }
        os << "], ";
        os << "Operations: [";
        for(int i = 0; i < Operations.size(); i++) {
            os << Operations[i] << ", ";
        }
        os << "], ";
        os << "}";
    }
};

class NChainReaction : public NExpression {
public:
    std::vector<std::shared_ptr<NChainReactionExpression>> Exprs;
    std::vector<int> Operators;

    void Add(NChainReactionExpression *InExpr, int InOperator = 0) {
        Exprs.emplace_back(InExpr);
        Operators.emplace_back(InOperator);
        assert(Exprs.size() == Operators.size());
    }

    virtual void Print(std::ostream& os) const override {
        assert(Exprs.size() == Operators.size());

        os << "ChainReaction: {";
        os << "Expr: [";
        for(int i = 0; i < Exprs.size(); i++) {
            Exprs[i]->Print(os); os << ", ";
        }
        os << "], ";
        os << "Operators: [";
        for(int i = 0; i < Exprs.size(); i++) {
            os << Operators[i] << ", ";
        }
        os << "]";
        os << "}";
    }

};

class NPathwayExpression : public NExpression {
public:
    int Type;
    std::shared_ptr<NExpression> Lhs;
    std::shared_ptr<NExpression> Rhs;

    NPathwayExpression(NPathwayExpression* InLhs, NPathwayExpression* InRhs, int InType)
        : Lhs(InLhs), Rhs(InRhs), Type(InType) {}

    virtual void Print(std::ostream& os) const override {
        os << "PathwayExpression: {";
        os << "Type: " << Type << ", ";
        Lhs->Print(os); os << ", ";
        Rhs->Print(os); os << ", ";
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPathwayDeclaration : public NStatement {
public:
    const NIdentifier Id;

    std::shared_ptr<NPathwayExpression> PathwayExpression;
    std::shared_ptr<NBlock> Block;
    std::shared_ptr<NChainReaction> PathwayChainReaction;

    NPathwayDeclaration(const NIdentifier& InId, NPathwayExpression* InPathwayExpression)
        : Id(InId), PathwayExpression(InPathwayExpression) {}

    /* New */
    NPathwayDeclaration(const NIdentifier& InId, NBlock* InBlock)
        : Id(InId), Block(InBlock) {};

    NPathwayDeclaration(const NIdentifier& InId, NChainReaction* InChainReaction)
    : Id(InId), PathwayChainReaction(InChainReaction) {};

    virtual void Print(std::ostream& os) const override {
        os << "PathwayDeclaration: {";
        Id.Print(os);
        if (PathwayExpression) {
            os << ", ";
            PathwayExpression->Print(os);
        }
        if (Block) {
            os << ", ";
            Block->Print(os);
        }
        if (PathwayChainReaction) {
            os << ", ";
            PathwayChainReaction->Print(os);
        }
        os << "}";
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
        os << "PathwayReactionID: {";
        Id.Print(os);
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPathwayReactionStatement : public NStatement {
public:
    NChainReaction PathwayExpression;

    NPathwayReactionStatement(NChainReaction& InPathwayExpression)
        : PathwayExpression(InPathwayExpression) {}


    virtual void Print(std::ostream& os) const override {
        os << "NPathwayExpression: {";
        PathwayExpression.Print(os);
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NOrganismDeclaration : public NStatement {
public:
    NIdentifier Id;
    std::string Description;
    std::shared_ptr<NBlock> Block;


    NOrganismDeclaration(const NIdentifier& InId, const std::string& InDescription)
        : Id(InId), Description(InDescription) {}
    NOrganismDeclaration(const NIdentifier& InId, NBlock* InBlock)
        : Id(InId), Block(InBlock) {}
    NOrganismDeclaration(const NIdentifier& InId)
        : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "Organism: {";
        Id.Print(os);
        if (Description.size() > 0) {
            os << ", ";
            os << "Description: " << Description;
        }
        if (Block) {
            os << ", ";
            Block->Print(os);
        }
        os << "}";
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
        os << "Experiment: {";
        Id.Print(os); os << ", ";
        os << "Description: " << Description;
        if (Block) {
            os << ", ";
            Block->Print(os);
        }
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NProteinCofactorStatement : public NStatement {
public:
    const NIdentifier Id;
    IdentifierList CofactorList;

    NProteinCofactorStatement(const NIdentifier& InId, const IdentifierList& InCofactorList) : Id(InId), CofactorList(InCofactorList) {};

    virtual void Print(std::ostream& os) const override {
        os << "Cofactor: {";
        os << Id.Name << ": [";
        for (const auto& arg: CofactorList) {
            os << arg->Name << ", ";
        }
        os << "]";
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NProteinDomainStatement : public NStatement {
public:
    const NIdentifier Id;
    IdentifierList DomainList;

    NProteinDomainStatement(const NIdentifier& InId, const IdentifierList& InDomainList) : Id(InId), DomainList(InDomainList) {};

    virtual void Print(std::ostream& os) const override {
        os << "Domain: {";
        os << Id.Name << ": [";
        for (const auto& domain: DomainList) {
            os << domain->Name << ", ";
        }
        os << "]";
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NStepStatement : public NStatement {
public:
    const NIdentifier Id;
    const NReaction Reaction;

    NStepStatement(const NIdentifier& InId, const NReaction& InReaction)
    : Id(InId), Reaction(InReaction) {}

    virtual void Print(std::ostream& os) const override {
        os << "Step: {";
        os << Id.Name << ": {";
        Reaction.Print(os);
        os << "}";
        os << "}";
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
        os << "Sequence: {";
        os << Id.Name << ": [";
        for (const auto& step: StepList) {
            os << step->Name << ", ";
        }
        os << "]";
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NInitiationStatement : public NStatement {
public:
    const NIdentifier Id;

    NInitiationStatement(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "Initiation: " << Id.Name;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NElongationStatement : public NStatement {
public:
    NReaction Reaction;

    NElongationStatement(const NReaction& InReaction) : Reaction(InReaction) {}

    virtual void Print(std::ostream& os) const override {
        os << "Elongation: {";
        Reaction.Print(os);
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NTerminationStatement : public NStatement {
public:
    const NIdentifier Id;

    NTerminationStatement(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "Termination: " << Id.Name;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NRibosomeDeclaration : public NStatement {
public:
    const NIdentifier Id;

    StatementList Statements;


    NRibosomeDeclaration(const NIdentifier& InId, const StatementList* InStatements) 
        : Id(InId) {
            Statements.insert(Statements.end(), InStatements->begin(), InStatements->end());
        }

    virtual void Print(std::ostream& os) const override {
        os << "Ribosome: {";
        Id.Print(os); os << ", ";
        for (const auto stmt : Statements) {
            stmt->Print(os); os << ", ";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NPolymeraseDeclaration : public NStatement {
public:
    const NIdentifier Id;

    StatementList Statements;

    NPolymeraseDeclaration(const NIdentifier& InId, const StatementList* InStatements)
        : Id(InId) 
		{Statements.insert(Statements.end(), InStatements->begin(), InStatements->end());
        }

    virtual void Print(std::ostream& os) const override {
        os << "Polymerase: {";
        Id.Print(os); os << ", ";
        for (const auto stmt : Statements) {
            stmt->Print(os); os << ", ";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};


class NRibosomeBindingSite : public NStatement {
public:
    const NIdentifier Id;

    NRibosomeBindingSite(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "RibosomeBindingSite: " << Id.Name;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NTranslationTerminator : public NStatement {
public:
    const NIdentifier Id;

    NTranslationTerminator(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "TranslationTerminator: " << Id.Name;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NReplicationOrigin : public NStatement {
public:
    const NIdentifier Id;

    NReplicationOrigin(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "ReplicationOrigin: " << Id.Name;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NReplicationTerminus : public NStatement {
public:
    const NIdentifier Id;

    NReplicationTerminus(const NIdentifier& InId) : Id(InId) {}

    virtual void Print(std::ostream& os) const override {
        os << "ReplicationTerminus: " << Id.Name;
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NUsingStatement : public NStatement {
public:
    int Type;
    const NIdentifier Id;

    NUsingStatement(int InType, const NIdentifier& InId) : Type(InType), Id(InId) {};

    virtual void Print(std::ostream& os) const override {
        os << "Using: {";
        os << "Type: " << Type; os << ", ";
        Id.Print(os); os << ", ";
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NDummyDeclaration : public NStatement {
public:
    const std::string StringLiteral;

    NDummyDeclaration(const std::string& InStringLiteral) : StringLiteral(InStringLiteral) {}

    virtual void Print(std::ostream& os) const override {
        os << "DummyDeclaration: " << StringLiteral;
    }
    virtual void Visit(FTraversalContext& Context) const override;

};

class NContainerStatement : public NStatement {
public:
    const NIdentifier Id;
    std::shared_ptr<NBlock> Body;

    NContainerStatement(const NIdentifier& InId, NBlock* InBody)
        : Id(InId), Body(InBody) {};
    NContainerStatement(const NIdentifier& InId)
        : Id(InId) {};

    virtual void Print(std::ostream& os) const override {
        os << "Container: {";
        Id.Print(os); os << ", ";
        if (Body) {
            for (const auto& stmt: Body->Statements) {
                stmt->Print(os); os << ", ";
            }
        }
        os << "}";
    }

};

class NPetridishStatement : public NContainerStatement {
public:
    NPetridishStatement(const NIdentifier& InId, NBlock* InBody)
        : NContainerStatement(InId, InBody) {};
    NPetridishStatement(const NIdentifier& InId)
        : NContainerStatement(InId) {};

    virtual void Print(std::ostream& os) const override {
        os << "Petridish: {";
        Id.Print(os); os << ", ";
        if (Body) {
            for (const auto& stmt: Body->Statements) {
                stmt->Print(os); os << ", ";
            }
        }
        os << "}";
    }
};

class NLoopStatement : public NStatement {
public:
    StatementList InitStatements;
    //std::shared_ptr<NStatement> InitStatement;

    std::shared_ptr<NExpression> CondExpression;
    std::shared_ptr<NExpression> LoopExpression;

    std::shared_ptr<NBlock> Body;

    // for-statement
    NLoopStatement(StatementList& InInitStatements, NExpression* InCondExpression, NExpression* InLoopExpression, NBlock* InBody)
        : InitStatements(InInitStatements), CondExpression(InCondExpression), LoopExpression(InLoopExpression), Body(InBody) {}

    // while-statement
    NLoopStatement(NExpression* InCondExpression, NBlock* InBody)
        : CondExpression(InCondExpression), Body(InBody) {}


    virtual void Print(std::ostream& os) const override {
        os << "Loop: {";
        if (InitStatements.size() > 0) {
            os << "Init: [";
            for (const auto& i: InitStatements) {
                i->Print(os); os << ", ";
            }
            //InitStatement->Print(os);
            os << "], ";
        }
        if (CondExpression) {
            os << "Cond: {";
            CondExpression->Print(os);
            os << "}, ";
        }
        if (LoopExpression) {
            os << "Loop: {";
            LoopExpression->Print(os);
            os << "}, ";
        }
        if (Body) {
            os << "Body: {";
            Body->Print(os);
        }
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override;
};

class NIfStatement : public NStatement {
public:
    std::shared_ptr<NExpression> CondExpression;
    std::shared_ptr<NBlock> Body;
    std::shared_ptr<NBlock> ElseBody;

    NIfStatement(NExpression* InCondExpression, NBlock* InBody)
        : CondExpression(InCondExpression)
        , Body(InBody) {}

    NIfStatement(NExpression* InCondExpression, NBlock* InBody, NBlock* InElseBody)
        : CondExpression(InCondExpression)
        , Body(InBody), ElseBody(InElseBody) {}


    virtual void Print(std::ostream& os) const override {
        os << "If: {";
        os << "Cond: {"; CondExpression->Print(os); os << "}";
        os << ", ";
        os << "Body: {"; Body->Print(os); os << "}";
        if (ElseBody) {
            os << ", ";
            os << "Else: {"; ElseBody->Print(os); os << "}";
        }
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override;
};

class NAExpression : public NExpression {
public:
    std::shared_ptr<NExpression> OpA;
    std::shared_ptr<NExpression> OpB;
    int Oper;

    NAExpression(int InOper, NExpression* InOpA)
        : Oper(InOper), OpA(InOpA) {}

    NAExpression(int InOper, NExpression* InOpA, NExpression* InOpB)
        : Oper(InOper), OpA(InOpA), OpB(InOpB) {}

    virtual void Print(std::ostream& os) const override {
        os << "Arithmetic: {";
        os << "OpCode: " << Oper;
        if (OpA) {
            os << ", ";
            os << "OpA: {";
            OpA->Print(os);
            os << "}";
        }
        if (OpB) {
            os << ", ";
            os << "OpB: {";
            OpB->Print(os);
            os << "}";
        }
        os << "}";
    }
    virtual void Visit(FTraversalContext& Context) const override {};
};

class NConstantExpression : public NExpression {
public:
    std::string Value;
    std::string Unit;

    NConstantExpression(const std::string& InValue)
        : Value(InValue) {}
    NConstantExpression(const std::string& InValue, const std::string& InUnit)
        : Value(InValue), Unit(InUnit) {}

    virtual void Print(std::ostream& os) const override {
        os << "Constant: " << Value;
        if (!Unit.empty()) {
            os << "(" << Unit << ")";
        }
    }
    virtual void Visit(FTraversalContext& Context) const override {};

    virtual std::string Evaluate() const override {
        return Value;
    }
};

class NVariableExpression : public NExpression {
public:
    std::shared_ptr<NExpression> Variable;
    std::shared_ptr<NExpression> Index;

    NVariableExpression(NExpression* InVariableExpr)
        : Variable(InVariableExpr) {}
    NVariableExpression(NExpression* InVariableExpr, NExpression* InIndexExpr)
        : Variable(InVariableExpr), Index(InIndexExpr) {}

    virtual void Print(std::ostream& os) const override {
        os << "Variable: {";
        os << "Name: ";
        Variable->Print(os);
        if (Index) {
            os << ", ";
            os << "Index: {";
            Index->Print(os);
            os << "}";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override {};

    virtual std::string Evaluate() const override {
        return Variable->Evaluate();
    }
};

class NRangeExpression : public NExpression {
public:
    std::shared_ptr<NExpression> Begin; /* Inclusive */
    std::shared_ptr<NExpression> End;   /* Exclusive */
    std::shared_ptr<NExpression> Step;

    NRangeExpression()
        : Begin(nullptr), End(nullptr), Step(nullptr) {};

    NRangeExpression(NExpression* InBegin, NExpression* InEnd, NExpression* InStep)
        : Begin(InBegin), End(InEnd), Step(InStep) {};


    virtual void Print(std::ostream& os) const override {
        os << "Range: {";
        if (Begin) {
            os << "Begin: {";
            Begin->Print(os);
            os << "}, ";
        }
        if (End) {
            os << "End: {";
            End->Print(os);
            os << "}, ";
        }
        if (Step) {
            os << "Step: {";
            Step->Print(os);
            os << "}, ";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NFunctionCallExpression : public NExpression {
public:
    std::shared_ptr<NExpression> Name;
    ExpressionList *Args;

    NFunctionCallExpression(NExpression* InName)
        : Name(InName), Args(nullptr) {};
    NFunctionCallExpression(NExpression* InName, ExpressionList* InArgs)
        : Name(InName), Args(InArgs) {};


    virtual void Print(std::ostream& os) const override {
        os << "FunctionCall: {";
        os << "Name: {";
        Name->Print(os);
        os << "}, ";
        if (Args) {
            os << "Args: [";
            for (const auto& arg: *Args) {
                arg->Print(os); os << ", ";
            }
            os << "]";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NExpressionStatement : public NStatement {
public:
    std::shared_ptr<NExpression> Expression;

    NExpressionStatement() {};
    NExpressionStatement(NExpression* InExpression)
        : Expression(InExpression) {};


    virtual void Print(std::ostream& os) const override {
        os << "Expression: {";
        Expression->Print(os);
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override;
};

class NInitializerExpression : public NExpression {
public:
    std::vector<std::shared_ptr<NExpression>> ExpressionList;

    NInitializerExpression() {};
    NInitializerExpression(NExpression* InExpr) {
        ExpressionList.emplace_back(InExpr);
    }

    void Append(NExpression* InExpr) {
        ExpressionList.emplace_back(InExpr);
    }


    virtual void Print(std::ostream& os) const override {
        os << "Initializer: [";
        for (const auto& item: ExpressionList) {
            item->Print(os);
            os << ", ";
        }
        os << "]";
    }

    virtual void Visit(FTraversalContext& Context) const override {};
};

class NDeclaraionStatement : public NStatement {
public:
    std::string Type;
    NIdentifier Id;
    std::shared_ptr<NInitializerExpression> Initializer;

    NDeclaraionStatement() {};
    NDeclaraionStatement(const std::string& InType, const NIdentifier& InId)
        : Type(InType), Id(InId) {};
    NDeclaraionStatement(const std::string& InType, const NIdentifier& InId, NExpression* InInit)
        : Type(InType), Id(InId), Initializer(std::make_shared<NInitializerExpression>(InInit)) {};
    NDeclaraionStatement(const std::string& InType, const NIdentifier& InId, NInitializerExpression* InInitList)
        : Type(InType), Id(InId), Initializer(InInitList) {};


    void SetType(const std::string& InType)
    {
        Type = InType;
    }
    void SetIdentifier(const NIdentifier& InId)
    {
        Id = InId;
    }
    void SetInitializer(NInitializerExpression* InInit)
    {
        Initializer = std::shared_ptr<NInitializerExpression>(InInit);
    }

    virtual void Print(std::ostream& os) const override {
        os << "Declaration: {";
        os << "Type: " << Type << ", ";
        os << "Name: "; Id.Print(os); os << ", ";
        if (Initializer) {
            os << "Value: [";
            Initializer->Print(os);
            os << "]";
        }
        os << "}";
    }

    virtual void Visit(FTraversalContext& Context) const override {};


    static void UpdateType(StatementList *Statements, const std::string& InType)
    {
        for (auto& stmt: *Statements) {
            NDeclaraionStatement* DeclStmt = dynamic_cast<NDeclaraionStatement *>(stmt.get());
            if (DeclStmt) {
                DeclStmt->SetType(InType);
            }
        }
    }
};


#endif /* LCC_NODE_H */
