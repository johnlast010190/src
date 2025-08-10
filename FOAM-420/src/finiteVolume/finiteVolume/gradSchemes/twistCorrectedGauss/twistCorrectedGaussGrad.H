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
    (c) 2016 Esi Ltd.

Class
    Foam::fv::twistCorrectedGaussGrad

Description
    Gradient scheme similar to Gauss but with skewness correction and
    correction for face twist by performing sum of skew corrections for face
    triangle decomposition.

    Usage:
    gradSchemes
    {
       grad(...) twistCorrectedGauss linear <initialCorrScheme> <nIterations>;
    }
    where
    <initialCorrScheme> is the scheme data for the gradient used to perform the
    initial skew correction. Defaults to "Gauss linear" if unspecified.
    <nIterations> is the number of iterations performed, where the calculated
    gradient is used as the skew correction gradient for the following
    iteration.
    Defaults to 2 if <initialCorrScheme> is not set, otherwise defaults to 1.


SourceFiles
    twistCorrectedGaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef twistCorrectedGaussGrad_H
#define twistCorrectedGaussGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/gradSchemes/taylorGauss/taylorGaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class twistCorrectedGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class twistCorrectedGaussGrad
:
    public fv::gradScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;

        tmp<gradScheme<Type>> tgradCorrectionScheme_;

        label nIterations_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        twistCorrectedGaussGrad(const twistCorrectedGaussGrad&);

        //- Disallow default bitwise assignment
        void operator=(const twistCorrectedGaussGrad&);


public:

    //- Runtime type information
    TypeName("twistCorrectedGauss");


    // Constructors

        //- Construct from mesh
        twistCorrectedGaussGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linear<Type>(mesh)),
            tgradCorrectionScheme_(new taylorGaussGrad<Type>(mesh)),
            nIterations_(1)
        {}

        //- Construct from Istream
        twistCorrectedGaussGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(nullptr),
            tgradCorrectionScheme_(nullptr),
            nIterations_(1)
        {
            tinterpScheme_ =
                tmp<surfaceInterpolationScheme<Type>>
                (
                    surfaceInterpolationScheme<Type>::New(mesh, db, is)
                );
            if (is.eof())
            {
                tgradCorrectionScheme_ =
                    tmp<gradScheme<Type>>
                    (
                        new taylorGaussGrad<Type>(mesh)
                    );
                nIterations_ = 1;
            }
            else
            {
                tgradCorrectionScheme_ =
                    tmp<gradScheme<Type>>
                    (
                        gradScheme<Type>::New(mesh, db, is)
                    );
                if (!is.eof())
                {
                    is >> nIterations_;
                }
                else
                {
                    nIterations_ = 1;
                }
            }
        }


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteVolume/gradSchemes/twistCorrectedGauss/twistCorrectedGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
