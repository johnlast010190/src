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
    (c) 2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::Function1Types::quarterCosineRamp

Description
    Quarter-cosine ramp function starting from 0 and increasing to 1 from \c
    start over the \c duration and remaining at 1 thereafter.

See also
    Foam::Function1Types::ramp

SourceFiles
    quarterCosineRamp.C

\*---------------------------------------------------------------------------*/

#ifndef quarterCosineRamp_H
#define quarterCosineRamp_H

#include "primitives/functions/Function1/ramp/ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class quarterCosineRamp Declaration
\*---------------------------------------------------------------------------*/

class quarterCosineRamp
:
    public ramp
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const quarterCosineRamp&);


public:

    // Runtime type information
    TypeName("quarterCosineRamp");


    // Constructors

        //- Construct from entry name and dictionary
        quarterCosineRamp
        (
            const word& entryName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~quarterCosineRamp();


    // Member Functions

        //- Return value for time t
        scalar value(const scalar t) const;

        //- Derivative
        scalar derivative(const scalar t) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
