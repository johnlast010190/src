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
    (c) 2020 Esi Ltd.

Class
    Foam::volWeighted

Description
    Interpolation scheme class that generates weightes from volume field

    Weights are calculated inverse:

    w = A[neightbour]/(A[owner] + A[neighbour])

SourceFiles
    volWeighted.C

\*---------------------------------------------------------------------------*/

#ifndef volWeighted_H
#define volWeighted_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class volWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const volWeighted&);


public:

    //- Runtime type information
    TypeName("volWeighted");


    // Constructors

        //- Construct from Istream
        volWeighted
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        volWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            tmp<surfaceScalarField> volFieldWeights
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "volWeights",
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    this->mesh(),
                    dimless
                )
            );

            const labelUList& owner = this->mesh().owner();
            const labelUList& neighbour = this->mesh().neighbour();

            const scalarField& vol = this->mesh().V();

            // ... and reference to the internal field of the weighting factors
            scalarField& w = volFieldWeights.ref().primitiveFieldRef();

            forAll(owner, facei)
            {
                w[facei]
                    = vol[neighbour[facei]]
                    /(vol[owner[facei]] + vol[neighbour[facei]]);
            }

            forAll(this->mesh().boundary(), patchI)
            {
                volFieldWeights.ref().boundaryFieldRef()[patchI] = 1.0;
            }

            return volFieldWeights;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
