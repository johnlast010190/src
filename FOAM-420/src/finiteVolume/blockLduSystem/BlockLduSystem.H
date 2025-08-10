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
    (c) Vuko Vukcevic, FMENA Zagreb.
    (c) 2019 Esi Ltd.

Class
    Foam::BlockLduSystem

Description
    BlockLduSystem is a wrapper for BlockLduMatrix with source field

SourceFiles
    BlockLduSystem.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduSystem_H
#define BlockLduSystem_H

#include "matrices/blockLduMatrix/BlockLduMatrix/blockLduMatrices.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class blockType, class sourceType>
class BlockLduSystem;

template<class blockType, class sourceType>
Ostream& operator<<(Ostream&, const BlockLduSystem<blockType, sourceType>&);


/*---------------------------------------------------------------------------*\
                        Class BlockLduSystem Declaration
\*---------------------------------------------------------------------------*/

template<class blockType, class sourceType>
class BlockLduSystem
:
    public BlockLduMatrix<blockType>
{
    // Private data

        //- Source term
        Field<sourceType> source_;

protected:

        PtrList<CoeffField<blockType>> internalCoeffs_;
        PtrList<Field<sourceType>> boundaryCoeffs_;

        //- Dimension set for each solving component.
        //  It is not like the segregated matrices that we solve only one field
        //  Here we solve all the equations together.
        //  So they can have different dimensions
        PtrList<dimensionSet> dimensions_;

        //- Face flux field for non-orthogonal correction
        mutable GeometricField<sourceType, fvsPatchField, surfaceMesh>
            *faceFluxCorrectionPtr_;

public:

    // Constructors

        //- Construct given addressing
        explicit BlockLduSystem(const fvMesh&);

        //- Construct given mesh and dimensions
        explicit BlockLduSystem(const fvMesh&, const dimensionSet&);

        //- Construct given addressing and source field
        BlockLduSystem(const fvMesh&, const Field<sourceType>&);

        //- Construct from components
        BlockLduSystem
        (
            const BlockLduMatrix<blockType>&,
            const Field<sourceType>&
        );

        //- Construct as copy
        BlockLduSystem(const BlockLduSystem<blockType, sourceType>&);


    //- Destructor
    virtual ~BlockLduSystem();

    // Member functions

        //- Operators specific to NS-eq.
            //- divU - laplacian(p) (like total-pressure BC)
            void addContinuityCoupledBC
            (
                const GeometricField<blockType, fvPatchField, volMesh>&,
                const volScalarField&,
                const surfaceTensorField&
            );


        //- Access

            Field<sourceType>& source()
            {
                return source_;
            }

            const Field<sourceType>& source() const
            {
                return source_;
            }

            const PtrList<CoeffField<blockType>>& internalCoeffs() const
            {
                return internalCoeffs_;
            }

            const PtrList<Field<sourceType>>& boundaryCoeffs() const
            {
                return boundaryCoeffs_;
            }

            PtrList<CoeffField<blockType>>& internalCoeffs()
            {
                return internalCoeffs_;
            }

            PtrList<Field<sourceType>>& boundaryCoeffs()
            {
                return boundaryCoeffs_;
            }

            typedef GeometricField<sourceType, fvsPatchField, surfaceMesh>
                *surfaceTypeFieldPtr;
            surfaceTypeFieldPtr& faceFluxCorrectionPtr()
            {
                return faceFluxCorrectionPtr_;
            }

            const surfaceTypeFieldPtr& faceFluxCorrectionPtr() const
            {
                return faceFluxCorrectionPtr_;
            }

            const PtrList<dimensionSet>& dimensionSets() const
            {
                return dimensions_;
            }

            PtrList<dimensionSet>& dimensionSets()
            {
                return dimensions_;
            }


    // Member operators

        void operator=(const BlockLduSystem<blockType, sourceType>&);

        void negate();

        void operator+=(const BlockLduSystem<blockType, sourceType>&);
        void operator-=(const BlockLduSystem<blockType, sourceType>&);

        void operator*=(const scalarField&);
        void operator*=(const scalar);


    // Ostream operator

        friend Ostream& operator<< <blockType, sourceType>
        (
            Ostream&,
            const BlockLduSystem<blockType, sourceType>&
        );
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
