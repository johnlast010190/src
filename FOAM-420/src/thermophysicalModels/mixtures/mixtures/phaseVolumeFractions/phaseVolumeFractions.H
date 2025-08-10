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
    (c) 2023 Esi Ltd

Class
    Foam::phaseVolumeFractions

Group
    grpReactionThermophysicalMixtures

Description
    Multi-component mixture.

    Provides a list of volume fraction fields and helper functions to
    query mixture composition.

SourceFiles
    phaseVolumeFractions.C

\*---------------------------------------------------------------------------*/

#ifndef phaseVolumeFractions_H
#define phaseVolumeFractions_H

#include "fields/volFields/volFields.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/typeInfo/typeInfo.H"
#include "phasesTable/phasesTable.H"
#include "mixtures/volumeMassFractions/volumeMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class phaseVolumeFractions Declaration
\*---------------------------------------------------------------------------*/

class phaseVolumeFractions
:
    public volumeMassFractions
{

protected:

    // Protected data

        //- Table of phases names
        phasesTable phases_;

        //- List of phases active flags
        List<bool> active_;

        //- 'Passive' phase where
        //  phase fraction = 1-(sum of other phase fractions)
        label passivePhaseIndex_;

        //- Phase volume fractions
        PtrList<volScalarField> alphas_;

        //- Return phase divergence
        PtrList<volScalarField> dgdts_;

        //- Combined volume fractions for visualisation
        //  (only created if nPhases > 2)
        autoPtr<volScalarField> combinedAlphas_;

        //- Get mesh through object/solution registry
        const fvMesh& getMesh(const objectRegistry& obr) const;


public:

    //- Run time type information
    TypeName("phaseVolumeFractions");


    //- The base class of the mixture
    typedef phaseVolumeFractions basicMixtureType;


    // Constructors

        //- Read fraction fields
        void createFields
        (
            const dictionary& materialPropertiesDict,
            const objectRegistry& obr
        );

        //- Construct from dictionary, phases names, mesh
        phaseVolumeFractions
        (
            const dictionary& materialPropertiesDict,
            const wordList& phasesNames,
            const objectRegistry& obr,
            const word& dummy = word::null
        );

        //- Construct from dict, object registry
        phaseVolumeFractions
        (
            const dictionary& materialPropertiesDict,
            const objectRegistry& obr,
            const word& dummy = word::null
        );


    //- Destructor
    virtual ~phaseVolumeFractions()
    {}


    // Member functions

        //- Return the table of phases
        inline const phasesTable& phases() const;

        //- Does the mixture include this phase?
        inline bool contains(const word& phaseName) const;

        //- Return true for active phases
        inline bool active(label phasei) const override;

        //- Return the bool list of active phases
        inline const List<bool>& active() const override;

        //- Set phasei active
        inline void setActive(label phasei);

        //- Set phasei inactive
        inline void setInactive(label phasei);

        //- Return passive phase index, or -1 if none
        virtual label passiveIndex() const override
        {
            return passivePhaseIndex_;
        }

        //- Return the volume-fraction fields
        inline PtrList<volScalarField>& alphas();

        //- Return the const volume-fraction fields
        inline const PtrList<volScalarField>& alphas() const;

        //- Return the volume fraction fields
        inline PtrList<volScalarField>& fractions() override;

         //- Return the const volume fraction fields
        inline const PtrList<volScalarField>& fractions() const override;

        //- Return the volume-fraction field for a phase given by index
        inline volScalarField& alphas(const label i);

        //- Return the const volume-fraction field for a phase given by index
        inline const volScalarField& alphas(const label i) const;

        //- Return the volume-fraction field for a phase given by name
        inline volScalarField& alphas(const word& phaseName);

        //- Return the const volume-fraction field for a phase given by name
        inline const volScalarField& alphas(const word& phaseName) const;

        //- Return phase divergence field
        const PtrList<volScalarField>& dgdts() const
        {
            return dgdts_;
        }

        //- Recreate combined alphas
        void recomputeCombinedAlphas();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseVolumeFractionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
