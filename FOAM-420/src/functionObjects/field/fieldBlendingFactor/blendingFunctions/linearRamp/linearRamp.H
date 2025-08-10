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
    (c) 2011 Esi Ltd.

Class
    Foam::linearRamp

Description
    linearRamp blending functions

SourceFiles
    linearRamp.C

\*---------------------------------------------------------------------------*/

#ifndef linearRamp_H
#define linearRamp_H

#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"
#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingFunctions
{

/*---------------------------------------------------------------------------*\
                          Class linearRamp Declaration
\*---------------------------------------------------------------------------*/

class linearRamp
:
    public blendingFunction
{
private:

    // Private Data

        //- minimum value for blending
        scalar minValue_;

        //- maximum value for blending
        scalar maxValue_;

        //- blending coeff at minimum value
        scalar minCoeff_;

        //- blending coeff at maximum value
        scalar maxCoeff_;

        //- blending polynomial
        Polynomial<2> linearF_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearRamp(const linearRamp&);

        //- Disallow default bitwise assignment
        void operator=(const linearRamp&);


public:

    //- Runtime type information
        TypeName("linear");

    // Constructors

        //- Construct from components
        linearRamp(const dictionary& data, bool dum = true);

    //- Destructor
    virtual ~linearRamp(){}


    // Member Functions

        //- Update blending function
        virtual void update(scalar);


    // Member Operators

        //- convert value into blending factor
        virtual scalar operator()(scalar) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
