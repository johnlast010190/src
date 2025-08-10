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
    Foam::weightBlended

Description
    Interpolation scheme that blends two sets of weights.

SourceFiles
    weightBlended.C

\*---------------------------------------------------------------------------*/

#ifndef weightBlended_H
#define weightBlended_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class weightBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class weightBlended
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        tmp<surfaceInterpolationScheme<Type>> tScheme2_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        weightBlended(const weightBlended&);

        //- Disallow default bitwise assignment
        void operator=(const weightBlended&);


public:

    //- Runtime type information
    TypeName("weightBlended");


    // Constructors

        //- Construct from mesh and Istream
        weightBlended
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, db, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, db, is)
            )
        {}


        //- Construct from mesh, faceFlux and Istream
        weightBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            tmp<surfaceScalarField> w1(tScheme1_().weights(vf));
            tmp<surfaceScalarField> w2(tScheme1_().weights(vf));
            tmp<surfaceScalarField> w12(w1*w2);
            dimensionedScalar one("one", dimless, 1.0);

            return ((w12)/(2.0*w12 - w1 - w2 + one));
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
