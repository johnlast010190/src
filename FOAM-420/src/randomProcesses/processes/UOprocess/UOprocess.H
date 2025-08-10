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

Class
    Foam::UOprocess

Description
    Random UO process.

SourceFiles
    UOprocess.C

\*---------------------------------------------------------------------------*/

#ifndef UOprocess_H
#define UOprocess_H

#include "fields/Fields/complexFields/complexFields.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Kmesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class UOprocess Declaration
\*---------------------------------------------------------------------------*/

class UOprocess
{
    // Private data

        Random GaussGen;

        const Kmesh& Mesh;
        const scalar DeltaT, RootDeltaT;
        complexVectorField UOfield;

        scalar Alpha;
        scalar Sigma;
        scalar Kupper;
        scalar Klower;
        scalar Scale;


    // Private Member Functions

        complexVector WeinerProcess();


public:

    // Constructors

        //- Construct from wavenumber mesh and timestep
        UOprocess
        (
            const Kmesh& kmesh,
            const scalar deltaT,
            const dictionary&
        );


    // Member Functions

        const complexVectorField& newField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
