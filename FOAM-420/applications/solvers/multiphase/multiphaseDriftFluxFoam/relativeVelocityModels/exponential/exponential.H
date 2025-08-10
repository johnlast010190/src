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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

Class
    Foam::exponential

Description
    exponential relative velocity model

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "relativeVelocityModel/relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public relativeVelocityModel
{
    // Private data

        //- A coefficient
        dimensionedScalar a_;

        //- Drift velocity
        dimensionedVector V0_;

        //- Stokes velocity (=max)
        dimensionedVector Vstokes_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("exponential");


    // Constructors

        //- Construct from components
        exponential
        (
            const word& phaseName,
            const dictionary& dict,
            const volScalarField& alpha,
            const scalar& alphaMax
        );


    //- Destructor
    ~exponential();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();

        // drift velocity monotonous so max at alpha=0
        virtual scalar alpha23() const
        {
            return 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
