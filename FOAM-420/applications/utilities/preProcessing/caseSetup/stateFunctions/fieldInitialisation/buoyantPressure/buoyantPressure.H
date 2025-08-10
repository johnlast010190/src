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
    (c) 2020 Esi Ltd.

Class
    Foam::fieldInits::buoyantPressure

Description
    buoyantPressure initialization method, solving the stationary
    pressure Laplace equation to initialise non-linear pressure fields
    for buoyant solvers (e.g. large air stacks in atmospheric flows).

SourceFiles
    buoyantPressure.C

\*--------------------------------------------------------------------*/

#ifndef buoyantPressure_H
#define buoyantPressure_H

#include "fieldInitialisation/fieldInit/fieldInit.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{

/*--------------------------------------------------------------------*\
                              Class buoyantPressure
\*--------------------------------------------------------------------*/

class buoyantPressure
:
    public fieldInit
{
public:

    //- Runtime type information
    TypeName("buoyantPressure");


    // Constructors

        //- construct from initMethod
        buoyantPressure
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );


    //- Destructor
    virtual ~buoyantPressure()
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
