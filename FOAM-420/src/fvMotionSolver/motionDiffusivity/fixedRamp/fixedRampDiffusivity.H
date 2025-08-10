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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::fixedRampDiffusivity

Description
    Ramp based diffusivity. Two additional coefficent specify:
    Lbuff - thickness of near-wall buffer zone
    DbuffRatio - ration of buffer zone diffusivity to outer diffusivity

    Uses cubic ramp of extent Lbuff between buffer and outer zones

SourceFiles
    fixedRampDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRampDiffusivity_H
#define fixedRampDiffusivity_H

#include "motionDiffusivity/uniform/uniformDiffusivity.H"
#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedRampDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class fixedRampDiffusivity
:
    public uniformDiffusivity
{
    // Private data

        //- Patches selected to base the distance on
        wordList patchNames_;

        //- Ramp coefficients

            //- buffer length
            scalar Lbuff_;

            //- buffer diffusivity
            scalar DbuffRatio_;

            //- blending polynomial
            Polynomial<4> cubicF_;

    // Private Member Functions

        //- Ramp function
        scalar fDiffusivity(scalar y);

        //- Disallow default bitwise copy construct
        fixedRampDiffusivity(const fixedRampDiffusivity&);

        //- Disallow default bitwise assignment
        void operator=(const fixedRampDiffusivity&);


public:

    //- Runtime type information
    TypeName("fixedRamp");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        fixedRampDiffusivity
        (
            const fvMesh& mSolver,
            Istream& mdData
        );


    //- Destructor
    virtual ~fixedRampDiffusivity();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
