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
    Foam::fv::compactGaussLaplacianScheme

Group
    grpFvLaplacianSchemes

Description
    (An)isotropic Laplacian using compact stencil formulation with over-
    relaxation for face non-orthogonality (if 'corrected' snGrad scheme is
    used).

SourceFiles
    compactGaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef compactGaussLaplacianScheme_H
#define compactGaussLaplacianScheme_H

#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                Class compactGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class compactGaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        compactGaussLaplacianScheme(const compactGaussLaplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const compactGaussLaplacianScheme&);

public:

    //- Runtime type information
    TypeName("compactGauss");


    // Constructors

        //- Construct null
        compactGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        compactGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            laplacianScheme<Type, GType>(mesh, db, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        compactGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    //- Destructor
    virtual ~compactGaussLaplacianScheme()
    {}


    // Member Functions

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
    #include "compactGaussLaplacianScheme/compactGaussLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
