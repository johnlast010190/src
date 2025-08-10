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
    (c) 2022 Esi Ltd.

Class
    Foam::fieldInits::equationOfState

Description
    equationOfState initialization method, calculating the density from
    the equation of state.

SourceFiles
    equationOfState.C

\*--------------------------------------------------------------------*/

#ifndef equationOfState_H
#define equationOfState_H

#include "fieldInitialisation/fieldInit/fieldInit.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{

/*--------------------------------------------------------------------*\
                              Class equationOfState
\*--------------------------------------------------------------------*/

class equationOfState
:
    public fieldInit
{
public:

    //- Runtime type information
    TypeName("equationOfState");


    // Constructors

        //- construct from initMethod
        equationOfState
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );


    //- Destructor
    virtual ~equationOfState()
    {}


    // Member Functions

        //- Correct the fieldInit model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldInits
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
