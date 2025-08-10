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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::solidReactionThermo

Description
    Foam::solidReactionThermo

SourceFiles
    solidReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidReactionThermo_H
#define solidReactionThermo_H

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "solidThermo/solidThermo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidReactionThermo Declaration
\*---------------------------------------------------------------------------*/


class solidReactionThermo
:
    public solidThermo
{

public:

    //- Runtime type information
    TypeName("solidReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidReactionThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );

     // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidReactionThermo,
        dictionary,
        (const objectRegistry& obr, const dictionary& dict, const word& phaseName),
        (obr, dict, phaseName)
    );

    // Constructors

        //- Construct from mesh and phase name
        solidReactionThermo
        (
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        solidReactionThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based onobjectRegistry
        static autoPtr<solidReactionThermo> New
        (
            const objectRegistry&,
            const word& phaseName=word::null
        );

         //- Standard selection based on objectRegistry amd dictionary
        static autoPtr<solidReactionThermo> New
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~solidReactionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
