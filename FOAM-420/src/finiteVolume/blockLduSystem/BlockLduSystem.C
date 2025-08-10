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

\*---------------------------------------------------------------------------*/

#include "blockLduSystem/BlockLduSystem.H"
#include "db/IOstreams/IOstreams.H"
#include "volMesh/volMesh.H"
#include "fields/GeometricFields/GeometricField/GeometricFields.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const fvMesh& mesh
)
:
    BlockLduMatrix<blockType>(mesh),
    source_(mesh.lduAddr().size(), pTraits<sourceType>::zero),
    internalCoeffs_(mesh.boundary().size()),
    boundaryCoeffs_(mesh.boundary().size()),
    dimensions_(pTraits<sourceType>::nComponents),
    faceFluxCorrectionPtr_(nullptr)
{}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const fvMesh& mesh,
    const dimensionSet& dSet
)
:
    BlockLduMatrix<blockType>(mesh),
    source_(mesh.lduAddr().size(), pTraits<sourceType>::zero),
    internalCoeffs_(mesh.boundary().size()),
    boundaryCoeffs_(mesh.boundary().size()),
    dimensions_(pTraits<sourceType>::nComponents),
    faceFluxCorrectionPtr_(nullptr)
{
    forAll(dimensions_, dI)
    {
        dimensions_.set(dI, new dimensionSet(dSet));
    }
}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const fvMesh& mesh,
    const Field<sourceType>& s
)
:
    BlockLduMatrix<blockType>(mesh),
    source_(s),
    internalCoeffs_(mesh.boundary().size()),
    boundaryCoeffs_(mesh.boundary().size()),
    dimensions_(pTraits<sourceType>::nComponents),
    faceFluxCorrectionPtr_(nullptr)
{
    if (mesh.lduAddr().size() != s.size())
    {
        FatalErrorInFunction
            << "Sizes of ldu addressing and source field are not the same."
            << abort(FatalError);
    }
}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const BlockLduMatrix<blockType>& bm,
    const Field<sourceType>& s
)
:
    BlockLduMatrix<blockType>(bm),
    source_(s),
    internalCoeffs_(),
    boundaryCoeffs_(),
    dimensions_(pTraits<sourceType>::nComponents),
    faceFluxCorrectionPtr_(nullptr)
{
    if (this->lduAddr().size() != s.size())
    {
        FatalErrorInFunction
            << "Sizes of block matrix and source field are not the same."
            << abort(FatalError);
    }
}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::BlockLduSystem
(
    const BlockLduSystem<blockType, sourceType>& bs
)
:
    BlockLduMatrix<blockType>(bs),
    source_(bs.source()),
    internalCoeffs_(bs.internalCoeffs().size()),
    boundaryCoeffs_(bs.boundaryCoeffs().size()),
    dimensions_(bs.dimensions_.size()),
    faceFluxCorrectionPtr_(nullptr)
{
    const PtrList<CoeffField<blockType>>& iC = bs.internalCoeffs();
    const PtrList<Field<sourceType>>& bC = bs.boundaryCoeffs();
    const PtrList<dimensionSet>& dimSetList = bs.dimensionSets();

    forAll(dimensions_, dI)
    {
        if (dimSetList.set(dI))
        {
            dimensions_.set(dI, new dimensionSet(dimSetList[dI]));
        }
    }

    forAll(iC, pI)
    {
        if (iC.set(pI))
        {
            internalCoeffs_.set(pI, new CoeffField<blockType>(iC[pI]));
        }
    }
    forAll(bC, pI)
    {
        if (bC.set(pI))
        {
            boundaryCoeffs_.set
            (
                pI,
                new Field<sourceType>
                (
                    bC[pI]
                )
            );
        }
    }
    if (bs.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
            GeometricField<sourceType, fvsPatchField, surfaceMesh>
            (
                *(bs.faceFluxCorrectionPtr_)
            );
    }
}


template<class blockType, class sourceType>
Foam::BlockLduSystem<blockType, sourceType>::~BlockLduSystem()
{
    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::addContinuityCoupledBC
(
    const GeometricField<blockType, fvPatchField, volMesh>& U,
    const volScalarField& p,
    const surfaceTensorField& rDf
)
{
    NotImplemented;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class blockType, class sourceType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    os  << static_cast<const BlockLduMatrix<blockType>&>(bs) << nl
        << bs.source() << endl;

    return os;
}


// ************************************************************************* //
