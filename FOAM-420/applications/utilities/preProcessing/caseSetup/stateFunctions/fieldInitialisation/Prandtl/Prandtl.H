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
    (c) 2016 Esi Ltd.

Class
    Foam::fieldInits::Prandtl

Description
    Prandtl initialization method.

SourceFiles
    Prandtl.C

\*--------------------------------------------------------------------*/

#ifndef Prandtl_H
#define Prandtl_H

#include "fieldInitialisation/fieldInit/fieldInit.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{

/*--------------------------------------------------------------------*\
                              Class Prandtl
\*--------------------------------------------------------------------*/

class Prandtl
:
    public fieldInit
{
    // Private data

    // Private Functions

public:

    //- Runtime type information
    TypeName("Prandtl");


    // Constructors

        //- construct from initMethod
        Prandtl
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );


    //- Destructor
    virtual ~Prandtl()
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
