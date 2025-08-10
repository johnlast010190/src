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
    (c) 2019-2022 Esi Ltd.

Class
    Foam::fv::convectionScheme

Group
    grpFvConvectionSchemes

Description
    Abstract base class for convection schemes.

SourceFiles
    convectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef convectionScheme_H
#define convectionScheme_H

#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class convectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class convectionScheme
:
    public tmp<convectionScheme<Type>>::refCount
{
    // Private data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            convectionScheme,
            Istream,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            convectionScheme,
            Multivariate,
            (
                const fvMesh& mesh,
                const typename multivariateSurfaceInterpolationScheme<Type>::
                    fieldTable& fields,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, fields, faceFlux, schemeData)
        );


    // Constructors

        //- Copy construct
        convectionScheme(const convectionScheme&);

        //- Construct from mesh, flux and Istream
        convectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField&
        )
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new convectionScheme created on freestore
        static tmp<convectionScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


        //- Return a pointer to a new multivariate convectionScheme
        //  created on freestore
        static tmp<convectionScheme<Type>> New
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~convectionScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        virtual tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        virtual tmp<fvBlockMatrix<Type>> fvmBDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;


    // Member operators

        void operator=(const convectionScheme<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables
// and instantiate explicitly

#define makeFvConvectionTypeScheme(SS, Type)                                   \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            convectionScheme<Type>::addIstreamConstructorToTable<SS<Type>>     \
                add##SS##Type##IstreamConstructorToTable_;                     \
            template class SS<Type>;                                           \
        }                                                                      \
    }

#define makeFvConvectionScheme(SS)                                             \
                                                                               \
makeFvConvectionTypeScheme(SS, scalar)                                         \
makeFvConvectionTypeScheme(SS, vector)                                         \
makeFvConvectionTypeScheme(SS, sphericalTensor)                                \
makeFvConvectionTypeScheme(SS, symmTensor)                                     \
makeFvConvectionTypeScheme(SS, tensor)


#define makeMultivariateFvConvectionTypeScheme(SS, Type)                       \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            convectionScheme<Type>::                                           \
                addMultivariateConstructorToTable<SS<Type>>                    \
                add##SS##Type##MultivariateConstructorToTable_;                \
            template class SS<Type>;                                           \
        }                                                                      \
    }


#define makeMultivariateFvConvectionScheme(SS)                                 \
                                                                               \
makeMultivariateFvConvectionTypeScheme(SS, scalar)                             \
makeMultivariateFvConvectionTypeScheme(SS, vector)                             \
makeMultivariateFvConvectionTypeScheme(SS, sphericalTensor)                    \
makeMultivariateFvConvectionTypeScheme(SS, symmTensor)                         \
makeMultivariateFvConvectionTypeScheme(SS, tensor)


#define makeFvConvectionAndMultivariateTypeScheme(SS, Type)          \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            convectionScheme<Type>::addIstreamConstructorToTable<SS<Type>>     \
                add##SS##Type##IstreamConstructorToTable_;                     \
            convectionScheme<Type>::                                           \
                addMultivariateConstructorToTable<SS<Type>>                    \
                add##SS##Type##MultivariateConstructorToTable_;                \
            template class SS<Type>;                                           \
        }                                                                      \
    }


#define makeFvConvectionAndMultivariateScheme(SS)                              \
                                                                               \
makeFvConvectionAndMultivariateTypeScheme(SS, scalar)                          \
makeFvConvectionAndMultivariateTypeScheme(SS, vector)                          \
makeFvConvectionAndMultivariateTypeScheme(SS, sphericalTensor)                 \
makeFvConvectionAndMultivariateTypeScheme(SS, symmTensor)                      \
makeFvConvectionAndMultivariateTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/convectionSchemes/convectionScheme/convectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
