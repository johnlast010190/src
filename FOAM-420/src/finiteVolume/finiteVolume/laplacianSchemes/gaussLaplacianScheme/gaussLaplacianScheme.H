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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2021 Esi Ltd.

Class
    Foam::fv::gaussLaplacianScheme

Group
    grpFvLaplacianSchemes

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    gaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussLaplacianScheme_H
#define gaussLaplacianScheme_H

#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class gaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class gaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Disallow default bitwise copy construct
        gaussLaplacianScheme(const gaussLaplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gaussLaplacianScheme&);


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacianScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            laplacianScheme<Type, GType>(mesh, db, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        gaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    //- Destructor
    virtual ~gaussLaplacianScheme()
    {}


    // Member Functions

        static tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp
        <
            fvBlockMatrix<Type>
        > fvmBLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};




// Use macros to emulate partial-specialisation of the the Laplacian functions
// for scalar diffusivity gamma

#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
inline tmp<fvMatrix<Type>> gaussLaplacianScheme<Type, scalar>::fvmLaplacian    \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<GeometricField<Type, fvPatchField, volMesh>>                               \
inline gaussLaplacianScheme<Type, scalar>::fvcLaplacian                        \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);


defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);



#define defineFvmScalarBLaplacianGamma(Type)                                   \
                                                                               \
template<>                                                                     \
inline tmp<fvBlockMatrix<scalar>> gaussLaplacianScheme<scalar, Type>::fvmBLaplacian   \
(                                                                              \
    const GeometricField<Type, fvsPatchField, surfaceMesh>&,                   \
    const GeometricField<scalar, fvPatchField, volMesh>&                       \
);                                                                             \


defineFvmScalarBLaplacianGamma(sphericalTensor);
defineFvmScalarBLaplacianGamma(symmTensor);
defineFvmScalarBLaplacianGamma(tensor);


#define defineFvmVectorBLaplacianGamma(Type)                                   \
                                                                               \
template<>                                                                     \
inline tmp<fvBlockMatrix<vector>> gaussLaplacianScheme<vector, Type>::fvmBLaplacian   \
(                                                                              \
    const GeometricField<Type, fvsPatchField, surfaceMesh>&,                   \
    const GeometricField<vector, fvPatchField, volMesh>& vf                    \
)                                                                              \
{                                                                              \
    tmp<fvBlockMatrix<vector>> tfvbm                                           \
    (                                                                          \
        new fvBlockMatrix<vector>                                              \
        (                                                                      \
            const_cast<GeometricField<vector, fvPatchField, volMesh>&>(vf)     \
        )                                                                      \
    );                                                                         \
    FatalErrorInFunction                                                       \
        << " BLaplacianscalar(scalar, scalar-vector) suppoted"                 \
        << abort(FatalError);                                                  \
    return tfvbm;                                                              \
}                                                                              \


defineFvmVectorBLaplacianGamma(sphericalTensor);
defineFvmVectorBLaplacianGamma(symmTensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.C"
    #include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussBlockLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
