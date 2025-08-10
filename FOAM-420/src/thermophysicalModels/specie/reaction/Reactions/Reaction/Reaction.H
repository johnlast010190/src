/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2021-2022 Esi Ltd.
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::Reaction

Description
    Simple extension of ReactionThermo to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    ReactionI.H
    Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Reaction_H
#define Reaction_H

#include "speciesTable/speciesTable.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "include/dummyThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ReactionThermo>
class Reaction;

template<class ReactionThermo>
inline Ostream& operator<<(Ostream&, const Reaction<ReactionThermo>&);


/*---------------------------------------------------------------------------*\
                           Class Reaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class Reaction
:
    public ReactionThermo::thermoType
{
protected:

    // Protected member functions

        //- Return string representation of the left of the reaction
        void reactionStrLeft(OStringStream& reaction) const;

        //- Return string representation of the right of the reaction
        void reactionStrRight(OStringStream& reaction) const;


public:

    // Static data

        //- Number of un-named reactions
        static label nUnNamedReactions;


    // Public data types

        //- Class to hold the specie index and its coefficients in the
        //  reaction rate expression
        struct specieCoeffs
        {
            label index;
            scalar stoichCoeff;
            scalar exponent;

            specieCoeffs()
            :
                index(-1),
                stoichCoeff(0),
                exponent(1)
            {}

            specieCoeffs(const speciesTable& species, Istream& is);

            bool operator==(const specieCoeffs& sc) const
            {
                return index == sc.index;
            }

            bool operator!=(const specieCoeffs& sc) const
            {
                return index != sc.index;
            }

            friend Ostream& operator<<(Ostream& os, const specieCoeffs& sc)
            {
                os  << sc.index << token::SPACE
                    << sc.stoichCoeff << token::SPACE
                    << sc.exponent;
                return os;
            }
        };


private:

    // Private data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const speciesTable& species_;

        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;

        //- Material tables pointer
        autoPtr<materialTables> mat_;

        //- Is it using new material library?
        mutable bool isMaterials_;

        //- Are the material model initialized?
        mutable bool isMatInit_;

        //- Temperature pointer
        mutable scalarField* T_;

        //- Pressure pointer
        mutable scalarField* p_;

        //- Mollar weights for the eqution left side
        mutable scalarField lhsWs_;

        //- Mollar weights for the eqution right side
        mutable scalarField rhsWs_;

        //- Chemical equation coefficents
        mutable scalar coeffs_;

        //- Pointer to the gibs energy models
        // for chem. equation left side
        mutable matScalarUTable lhsGModels_;

        //- Pointer to the gibs energy models
        // for chem. equation right side
        mutable matScalarUTable rhsGModels_;

        //- Phase name
        mutable word phaseName_;


    // Private Member Functions

        //- Return string representation of reaction
        string reactionStr(OStringStream& reaction) const;

        //- Construct reaction thermo
        void setThermo(const HashPtrTable<ReactionThermo>& thermoDatabase);

        //- Disallow default bitwise assignment
        void operator=(const Reaction<ReactionThermo>&) = delete;

        //- Return new reaction ID for un-named reactions
        label getNewReactionID();

        //- Initialize material models
        void initMaterialModels() const;


public:

    //- Runtime type information
    TypeName("Reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            dictionary,
            (
                const speciesTable& species,
                const HashPtrTable<ReactionThermo>& thermoDatabase,
                const dictionary& dict
            ),
            (species, thermoDatabase, dict)
        );


    // Constructors

        //- Construct from components
        Reaction
        (
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            bool initReactionThermo = true
        );

        //- Construct as copy given new speciesTable
        Reaction
        (
            const Reaction<ReactionThermo>&,
            const speciesTable& species
        );

        //- Construct from dictionary
        // NOTE: initReactionThermo is used by solidReaction where there is no
        // need of setting a lhs - rhs thermo type for each reaction. This is
        // needed for mechanism with reversible reactions
        Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict,
            bool initReactionThermo = true
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ReactionThermo>> clone() const //
        {
            return autoPtr<Reaction<ReactionThermo>>
            (
                new Reaction<ReactionThermo>(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction<ReactionThermo>>
            (
                new Reaction<ReactionThermo>(*this, species)
            );
        }


    // Selectors

        //- Return a pointer to new patchField created on freestore from dict
        static autoPtr<Reaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Reaction()
    {}


    // Member Functions

        // Access

            inline word& name();
            inline const word& name() const;

            // - Access to basic components of the reaction
            inline const List<specieCoeffs>& lhs() const;
            inline const List<specieCoeffs>& rhs() const;

            // - Access to gas components of the reaction
            virtual const List<specieCoeffs>& grhs() const;
            virtual const List<specieCoeffs>& glhs() const;

            //- Access to specie list
            const speciesTable& species() const;

            //- Access to gas specie list
            virtual const speciesTable& gasSpecies() const;

            //- Construct the left- and right-hand-side reaction coefficients
            void setLRhs
            (
                Istream&,
                const speciesTable&,
                List<specieCoeffs>& lhs,
                List<specieCoeffs>& rhs
            );


        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Test
            virtual scalar Kc
            (
                const scalar p,
                const scalar T
            ) const;


        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionThermo>
        (
            Ostream&,
            const Reaction<ReactionThermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reaction/Reactions/Reaction/ReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reaction/Reactions/Reaction/Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
