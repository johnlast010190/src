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

Class
    Foam::limitedSurfaceInterpolationScheme

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    Abstract base class for limited surface interpolation schemes.

SourceFiles
    limitedSurfaceInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef limitedSurfaceInterpolationScheme_H
#define limitedSurfaceInterpolationScheme_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class limitedSurfaceInterpolationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedSurfaceInterpolationScheme
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow copy construct
        limitedSurfaceInterpolationScheme
        (
            const limitedSurfaceInterpolationScheme&
        );

        //- Disallow default bitwise assignment
        void operator=(const limitedSurfaceInterpolationScheme&);


protected:

    // Protected data

        const surfaceScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("limitedSurfaceInterpolationScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            limitedSurfaceInterpolationScheme,
            Mesh,
            (
                const fvMesh& mesh,
                const objectRegistry& db,
                Istream& schemeData
            ),
            (mesh, db, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            limitedSurfaceInterpolationScheme,
            MeshFlux,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh and faceFlux
        limitedSurfaceInterpolationScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        limitedSurfaceInterpolationScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                db.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<limitedSurfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<limitedSurfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~limitedSurfaceInterpolationScheme();


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        //- Return the interpolation weighting factors for the given field,
        //  by limiting the given weights with the given limiter
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const surfaceScalarField& CDweights,
            tmp<surfaceScalarField> tLimiter
        ) const;

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return the interpolation weighting factors
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        flux
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makelimitedSurfaceInterpolationTypeScheme(SS, Type)                    \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>      \
    add##SS##Type##MeshFluxConstructorToTable_;                                \
                                                                               \
limitedSurfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>   \
    add##SS##Type##MeshConstructorToLimitedTable_;                             \
                                                                               \
limitedSurfaceInterpolationScheme<Type>::                                      \
    addMeshFluxConstructorToTable<SS<Type>>                                    \
    add##SS##Type##MeshFluxConstructorToLimitedTable_;

#define makelimitedSurfaceInterpolationScheme(SS)                              \
                                                                               \
makelimitedSurfaceInterpolationTypeScheme(SS, scalar)                          \
makelimitedSurfaceInterpolationTypeScheme(SS, vector)                          \
makelimitedSurfaceInterpolationTypeScheme(SS, sphericalTensor)                 \
makelimitedSurfaceInterpolationTypeScheme(SS, symmTensor)                      \
makelimitedSurfaceInterpolationTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
