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
    (c) 2013 Esi Ltd.

Class
    Foam::tanhRamp

Description
    tanhRamp blending functions
    min/maxValue is mapped to x = +/-3
    Output is mapped from -1/+1 to user supplied coefficients

SourceFiles
    tanhRamp.C

\*---------------------------------------------------------------------------*/

#ifndef tanhRamp_H
#define tanhRamp_H

#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tanhRamp Declaration
\*---------------------------------------------------------------------------*/

class tanhRamp
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

        //- linear input/output scaling coefficients
        scalar inOffset_;
        scalar inScale_;
        scalar outOffset_;
        scalar outScale_;

        //- switch to Pade approximation of tanh
        Switch approxTanh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        tanhRamp(const tanhRamp&);

        //- Disallow default bitwise assignment
        void operator=(const tanhRamp&);


public:

    //- Runtime type information
        TypeName("tanh");

    // Constructors

        //- Construct from components
        tanhRamp(const dictionary& data, bool dum = true);

    //- Destructor
    virtual ~tanhRamp(){}


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
