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
    Foam::weighted

Group
    grpFvSurfaceInterpolationSchemes

Description
    Interpolation scheme class using weights looked-up from the objectRegistry.

SourceFiles
    weighted.C

\*---------------------------------------------------------------------------*/

#ifndef weighted_H
#define weighted_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class weighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class weighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const surfaceScalarField& weights_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const weighted&);


public:

    //- Runtime type information
    TypeName("weighted");


    // Constructors

        //- Construct from weights
        weighted(const surfaceScalarField& weights)
        :
            surfaceInterpolationScheme<Type>(weights.mesh()),
            weights_(weights)
        {}

        //- Construct from Istream
        weighted
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            weights_
            (
                db.objectRegistry::template
                lookupObject<const surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        weighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            weights_
            (
                this->mesh().objectRegistry::template
                lookupObject<const surfaceScalarField>
                (
                    word(is)
                )
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return weights_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
