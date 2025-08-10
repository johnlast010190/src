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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::speciesMassFractions

Group
    grpReactionThermophysicalMixtures

Description
    Multi-component mixture.

    Provides a list of mass fraction fields and helper functions to
    query mixture composition.

SourceFiles
    speciesMassFractions.C

\*---------------------------------------------------------------------------*/

#ifndef speciesMassFractions_H
#define speciesMassFractions_H

#include "fields/volFields/volFields.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/typeInfo/typeInfo.H"
#include "speciesTable/speciesTable.H"
#include "mixtures/volumeMassFractions/volumeMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class speciesMassFractions Declaration
\*---------------------------------------------------------------------------*/

class speciesMassFractions
:
    public volumeMassFractions
{

protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- List of specie active flags
        List<bool> active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;

        //- Get mesh through object/solution registry
        const fvMesh& getMesh(const objectRegistry& obr) const;


public:

    //- Run time type information
    TypeName("speciesMassFractions");


    //- The base class of the mixture
    typedef speciesMassFractions basicMixtureType;


    // Constructors

        //- Read fraction fields
        void createFields
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const objectRegistry& obr,
            const word& phaseName
        );

        //- Construct from dictionary, species names, mesh and phase name
        speciesMassFractions
        (
            const dictionary&,
            const wordList& specieNames,
            const objectRegistry&,
            const word&
        );

        //- Construct from dict, object registry and phase name
        speciesMassFractions
        (
            const dictionary& thermoDict,
            const objectRegistry& obr,
            const word& phaseName
        );


    //- Destructor
    virtual ~speciesMassFractions()
    {}


    // Member functions

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return true for active species
        inline bool active(label speciei) const override;

        //- Return the bool list of active species
        inline const List<bool>& active() const override;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass/volume fraction fields
        inline PtrList<volScalarField>& fractions() override;

         //- Return the const mass/volume fraction fields
        inline const PtrList<volScalarField>& fractions() const override;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "speciesMassFractionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
