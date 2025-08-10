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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::algebraicSlipModel

Description
    algebraicSlipModel relative velocity model

SourceFiles
    algebraicSlipModel.C

\*---------------------------------------------------------------------------*/

#ifndef algebraicSlipModel_H
#define algebraicSlipModel_H

#include "relativeVelocityModel/relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class algebraicSlipModel Declaration
\*---------------------------------------------------------------------------*/

class algebraicSlipModel
:
    public relativeVelocityModel
{
    // Private data

        //- A coefficient
        const dimensionedScalar dd_;

        //- Drift velocity
        const dimensionedScalar Cd_;

        //- Mixture velocity
        const volVectorField& U_;

        //- momentum flux
        const surfaceScalarField& rhoPhi_;

        //- read gravitational constant
        const uniformDimensionedVectorField g_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("algebraicSlipModel");


    // Constructors

        //- Construct from components
        algebraicSlipModel
        (
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        );


    //- Destructor
    ~algebraicSlipModel();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
