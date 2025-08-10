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
    (c) 2010-2016 Esi Ltd.

Class
    Foam::omegaGaussLaplacianScheme

Description
    omegaGaussLaplacianScheme interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction.

SourceFiles
    omegaGaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef omegaGaussLaplacianScheme_H
#define omegaGaussLaplacianScheme_H

#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"
#include "turbulenceModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{
/*---------------------------------------------------------------------------*\
                           Class omegaGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class omegaGaussLaplacianScheme
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
        omegaGaussLaplacianScheme(const omegaGaussLaplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const omegaGaussLaplacianScheme&);

        scalar getDeltaStarCoefficient
        (
            const label& patchID,
            const label& patchFaceLabel,
            const label& faceLabel,
            const scalar& yPlus
        );
public:

    //- Runtime type information
    TypeName("omegaGauss");


    // Constructors
        //- Construct null
        omegaGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from faceFlux
        omegaGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        omegaGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            laplacianScheme<Type, GType>(mesh)
        {
        }

        //- Construct from mesh, interpolation and snGradScheme schemes
        omegaGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}

        virtual ~omegaGaussLaplacianScheme(){}
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

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        )
        {
            NotImplemented;
        }

        tmp
        <
            fvBlockMatrix<Type>
        > fvmBLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        )
        {
            FatalErrorInFunction << "Not Implemented " << abort(FatalError);
            return tmp<fvBlockMatrix<Type>>();
        }


        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "omegaGaussLaplacianScheme.C"
#endif

#endif

// ************************************************************************* //
