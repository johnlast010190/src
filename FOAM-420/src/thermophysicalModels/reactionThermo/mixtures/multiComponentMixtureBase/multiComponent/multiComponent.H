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
    Foam::multiComponent

Group
    grpReactionThermophysicalMixtures

Description
    Foam::multiComponent

SourceFiles
    multiComponent.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponent_H
#define multiComponent_H

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponent Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponent
:
public multiComponentMixtureBase<ThermoType>
{
    // Private data


    // Private Member Functions

        //- Construct as copy (not implemented)
        multiComponent(const multiComponent<ThermoType>&);


public:


    //- Runtime type information
    TypeName("multiComponent");


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  obr and phase name
        multiComponent
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const objectRegistry&,
            const speciesTable& species,
            PtrList<volScalarField>& Y,
            const word&
        );

        //- Construct from dictionary, obr and phase name
        multiComponent
        (
            const dictionary&,
            const objectRegistry&,
            const speciesTable& species,
            PtrList<volScalarField>& Y,
            const word&
        );


    //- Destructor
    virtual ~multiComponent()
    {}


    // Member functions

        virtual const ThermoType& cellMixture(const label celli) const;

        virtual const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        virtual const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        virtual const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtures/multiComponentMixtureBase/multiComponent/multiComponent.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
