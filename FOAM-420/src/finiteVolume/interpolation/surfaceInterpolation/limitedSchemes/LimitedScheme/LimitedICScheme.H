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
    (c) 2019 Esi Ltd.

Class
    Foam::LimitedICScheme

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    Class to create NVD/TVD limited weighting-factors.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches). The weight function is
    supplied the central-differencing weighting factor, the face-flux, the
    cell and face gradients (from which the normalised variable
    distribution may be created) and the cell centre distance.

    This code organisation is both neat and efficient, allowing for
    convenient implementation of new schemes to run on parallelised cases.

SourceFiles
    LimitedICScheme.C

\*---------------------------------------------------------------------------*/

#ifndef LimitedICScheme_H
#define LimitedICScheme_H

#include "interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/LimitFuncs.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/NVDTVD.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/NVDVTVDV.H"

#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/blendingFactors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class LimitedICScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Limiter, template<class> class LimitFunc, class BlendingFactor>
class LimitedICScheme
:
    public limitedSurfaceInterpolationScheme<Type>,
    public Limiter
{
    // Private Data Member
        word fieldName_;


    // Private Member Functions

        //- Calculate the limiter
        void calcLimiter
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            surfaceScalarField& limiterField
        ) const;

        //- Disallow default bitwise copy construct
        LimitedICScheme(const LimitedICScheme&);

        //- Disallow default bitwise assignment
        void operator=(const LimitedICScheme&);


public:

    //- Runtime type information
    TypeName("LimitedICScheme");

    typedef Limiter LimiterType;

    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        LimitedICScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const Limiter& weight
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(weight),
            fieldName_("alpha1")
        {}

        //- Construct from mesh, database and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the db objectRegistry
        LimitedICScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, db, is),
            Limiter(is),
            fieldName_(word(is))
        {}

        //- Construct from mesh, faceFlux and Istream
        LimitedICScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(is),
            fieldName_(word(is))
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeLimitedICSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    BLENFAC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef LimitedICScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>              \
    LimitedICScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##BLENFAC##_;                          \
defineTemplateTypeNameAndDebugWithName                                         \
    (LimitedICScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##BLENFAC##_, #SS, 0);                \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<LimitedICScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>                    \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<LimitedICScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>                    \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                       \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable             \
<LimitedICScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>                    \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                    \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable         \
<LimitedICScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>                    \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


#define makeLimitedICSurfaceInterpolationScheme(SS, LIMITER)                     \
                                                                               \
makeLimitedICSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,null,scalar)     \
makeLimitedICSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,null,vector)     \
makeLimitedICSurfaceInterpolationTypeScheme                                      \
(                                                                              \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    magSqr,                                                                    \
    null,                                                                      \
    sphericalTensor                                                            \
)                                                                              \
makeLimitedICSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,null,symmTensor)\
makeLimitedICSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,null,tensor)


#define makeLimitedICVSurfaceInterpolationScheme(SS, LIMITER)                    \
makeLimitedICSurfaceInterpolationTypeScheme(SS,LIMITER,NVDVTVDV,null,null,vector)


#define makeLLimitedICSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    LLIMITER,                                                                  \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    BLENFAC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef LimitedICScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>    \
    LimitedICScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##BLENFAC##_;                \
defineTemplateTypeNameAndDebugWithName                                         \
    (LimitedICScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##BLENFAC##_, #SS, 0);      \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<LimitedICScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>          \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<LimitedICScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>          \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                       \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable             \
<LimitedICScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>          \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                    \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable         \
<LimitedICScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC, blendingFactors::BLENFAC>>          \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/LimitedICScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
