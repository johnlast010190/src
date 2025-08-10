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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Class
    faLaplacianScheme

Description
    Abstract base class for laplacian schemes.

SourceFiles
    faLaplacianScheme.C
    faLaplacianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef faLaplacianScheme_H
#define faLaplacianScheme_H

#include "memory/tmp/tmp.H"
#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"
#include "interpolation/edgeInterpolation/schemes/linear/linearEdgeInterpolation.H"
#include "finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrad.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class faMatrix;

class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class laplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class laplacianScheme
:
    public refCount
{

protected:

    // Protected data

        const faMesh& mesh_;
        tmp<edgeInterpolationScheme<GType>> tinterpGammaScheme_;
        tmp<lnGradScheme<Type>> tlnGradScheme_;


    // Private Member Functions

        //- Disallow copy construct
        laplacianScheme(const laplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const laplacianScheme&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            laplacianScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        laplacianScheme(const faMesh& mesh)
        :
            mesh_(mesh),
            tinterpGammaScheme_(new linearEdgeInterpolation<GType>(mesh)),
            tlnGradScheme_(new correctedLnGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        laplacianScheme(const faMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpGammaScheme_(nullptr),
            tlnGradScheme_(nullptr)
        {
            tinterpGammaScheme_ = tmp<edgeInterpolationScheme<GType>>
            (
                edgeInterpolationScheme<GType>::New(mesh, is)
            );

            tlnGradScheme_ = tmp<lnGradScheme<Type>>
            (
                lnGradScheme<Type>::New(mesh, is)
            );
        }


        //- Construct from mesh, interpolation and snGradScheme scheme
        laplacianScheme
        (
            const faMesh& mesh,
            const tmp<edgeInterpolationScheme<GType>>& igs,
            const tmp<lnGradScheme<Type>>& sngs
        )
        :
            mesh_(mesh),
            tinterpGammaScheme_(igs),
            tlnGradScheme_(sngs)
        {}

    // Selectors

        //- Return a pointer to a new laplacianScheme created on freestore
        static tmp<laplacianScheme<Type, GType>> New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~laplacianScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<faMatrix<Type>> famLaplacian
        (
//            const edgeScalarField&,
            const GeometricField<GType, faePatchField, faEdgeMesh>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famLaplacian
        (
//            const areaScalarField&,
            const GeometricField<GType, faPatchField, areaMesh>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
//            const edgeScalarField&,
            const GeometricField<GType, faePatchField, faEdgeMesh>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
//            const areaScalarField&,
            const GeometricField<GType, faPatchField, areaMesh>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables
#define makeFaLaplacianTypeScheme(SS, GType, Type)                             \
    typedef Foam::fa::SS<Foam::Type, Foam::GType> SS##Type##GType;             \
    defineNamedTemplateTypeNameAndDebug(SS##Type##GType, 0);                   \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fa                                                           \
        {                                                                      \
            typedef SS<Type, GType> SS##Type##GType;                           \
                                                                               \
            laplacianScheme<Type, GType>::                                     \
                addIstreamConstructorToTable<SS<Type, GType>>                 \
                add##SS##Type##GType##IstreamConstructorToTable_;              \
        }                                                                      \
    }


#define makeFaLaplacianScheme(SS)                                              \
                                                                               \
makeFaLaplacianTypeScheme(SS, scalar, scalar)                                  \
makeFaLaplacianTypeScheme(SS, symmTensor, scalar)                              \
makeFaLaplacianTypeScheme(SS, tensor, scalar)                                  \
makeFaLaplacianTypeScheme(SS, scalar, vector)                                  \
makeFaLaplacianTypeScheme(SS, symmTensor, vector)                              \
makeFaLaplacianTypeScheme(SS, tensor, vector)                                  \
makeFaLaplacianTypeScheme(SS, scalar, sphericalTensor)                         \
makeFaLaplacianTypeScheme(SS, symmTensor, sphericalTensor)                     \
makeFaLaplacianTypeScheme(SS, tensor, sphericalTensor)                         \
makeFaLaplacianTypeScheme(SS, scalar, symmTensor)                              \
makeFaLaplacianTypeScheme(SS, symmTensor, symmTensor)                          \
makeFaLaplacianTypeScheme(SS, tensor, symmTensor)                              \
makeFaLaplacianTypeScheme(SS, scalar, tensor)                                  \
makeFaLaplacianTypeScheme(SS, symmTensor, tensor)                              \
makeFaLaplacianTypeScheme(SS, tensor, tensor)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/laplacianSchemes/faLaplacianScheme/faLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
