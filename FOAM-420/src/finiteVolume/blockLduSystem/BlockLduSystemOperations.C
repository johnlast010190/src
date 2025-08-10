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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::negate()
{
    BlockLduMatrix<blockType>::negate();
    source_.negate();
    forAll(internalCoeffs_, pI)
    {
        if (internalCoeffs_.set(pI))
        {
            internalCoeffs_[pI].negate();
        }
    }
    forAll(boundaryCoeffs_, pI)
    {
        if (boundaryCoeffs_.set(pI))
        {
            boundaryCoeffs_[pI].negate();
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator=
(
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    if (this == &bs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    BlockLduMatrix<blockType>::operator=(bs);
    source_ = bs.source();
    const PtrList<CoeffField<blockType>>& iC = bs.internalCoeffs();
    const PtrList<Field<sourceType>>& bC = bs.boundaryCoeffs();
    forAll(iC, pI)
    {
        if (iC.set(pI))
        {
            internalCoeffs_[pI] = iC[pI];
        }
    }
    forAll(bC, pI)
    {
        if (bC.set(pI))
        {
            boundaryCoeffs_[pI] = bC[pI];
        }
    }
    if (faceFluxCorrectionPtr_ && bs.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *bs.faceFluxCorrectionPtr_;
    }
    else if (bs.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<sourceType, fvsPatchField, surfaceMesh>
        (*bs.faceFluxCorrectionPtr_);
    }


    const PtrList<dimensionSet>& dimSetList = bs.dimensionSets();
    forAll(dimSetList, dI)
    {
        if (dimSetList.set(dI))
        {
            if (!dimensions_.set(dI))
            {
                dimensions_.set(dI, new dimensionSet(dimSetList[dI]));
            }
            else
            {
                if (dimensions_[dI] != dimSetList[dI])
                {
                    FatalErrorInFunction
                        << "Inconsistent dimensions"
                        << dimSetList[dI] << " = "
                        << dimensions_[dI]
                        << abort(FatalError);
                }
            }
        }
    }
}


template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator+=
(
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    BlockLduMatrix<blockType>::operator+=(bs);
    source_ += bs.source();
    const PtrList<CoeffField<blockType>>& iC = bs.internalCoeffs();
    const PtrList<Field<sourceType>>& bC = bs.boundaryCoeffs();
    forAll(iC, pI)
    {
        if (iC.set(pI))
        {
            if (internalCoeffs_.set(pI))
            {
                internalCoeffs_[pI] += iC[pI];
            }
            else
            {
                internalCoeffs_.set
                (
                    pI, new CoeffField<blockType>(iC[pI].size())
                );
                internalCoeffs_[pI] = iC[pI];
            }
        }
    }
    forAll(bC, pI)
    {
        if (bC.set(pI))
        {
            if (boundaryCoeffs_.set(pI))
            {
                boundaryCoeffs_[pI] += bC[pI];
            }
            else
            {
                boundaryCoeffs_.set(pI, new Field<sourceType>(bC[pI]));
            }
        }
    }
    if (faceFluxCorrectionPtr_ && bs.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *bs.faceFluxCorrectionPtr_;
    }
    else if (bs.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<sourceType, fvsPatchField, surfaceMesh>
        (*bs.faceFluxCorrectionPtr_);
    }

    const PtrList<dimensionSet>& dimSetList = bs.dimensionSets();
    forAll(dimSetList, dI)
    {
        if (dimSetList.set(dI) && dimensions_.set(dI))
        {
            if (dimensions_[dI] != dimSetList[dI])
            {
                FatalErrorInFunction
                    << "Inconsistent dimensions"
                    << dimSetList[dI] << " += "
                    << dimensions_[dI]
                    << abort(FatalError);
            }
        }
    }
}

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator-=
(
    const BlockLduSystem<blockType, sourceType>& bs
)
{
    BlockLduMatrix<blockType>::operator-=(bs);
    source_ -= bs.source();
    const PtrList<CoeffField<blockType>>& iC = bs.internalCoeffs();
    const PtrList<Field<sourceType>>& bC = bs.boundaryCoeffs();
    forAll(iC, pI)
    {
        if (iC.set(pI))
        {
            internalCoeffs_[pI] -= iC[pI];
        }
    }
    forAll(bC, pI)
    {
        if (bC.set(pI))
        {
            boundaryCoeffs_[pI] -= bC[pI];
        }
    }
    if (faceFluxCorrectionPtr_ && bs.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *bs.faceFluxCorrectionPtr_;
    }
    else if (bs.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<sourceType, fvsPatchField, surfaceMesh>
        (-*bs.faceFluxCorrectionPtr_);
    }

    const PtrList<dimensionSet>& dimSetList = bs.dimensionSets();
    forAll(dimSetList, dI)
    {
        if (dimSetList.set(dI) && dimensions_.set(dI))
        {
            if (dimensions_[dI] != dimSetList[dI])
            {
                FatalErrorInFunction
                    << "Inconsistent dimensions"
                    << dimSetList[dI] << " -= "
                    << dimensions_[dI]
                    << abort(FatalError);
            }
        }
    }
}

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator*=
(
    const scalarField& sf
)
{
    BlockLduMatrix<blockType>::operator*=(sf);
    source_ *= sf;
}

template<class blockType, class sourceType>
void Foam::BlockLduSystem<blockType, sourceType>::operator*=
(
    const scalar s
)
{
    BlockLduMatrix<blockType>::operator*=(s);
    source_ *= s;
    forAll(internalCoeffs_, pI)
    {
        if (internalCoeffs_.set(pI))
        {
            internalCoeffs_[pI] *= s;
        }
    }
    forAll(boundaryCoeffs_, pI)
    {
        if (boundaryCoeffs_.set(pI))
        {
            boundaryCoeffs_[pI] *= s;
        }
    }
    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


// ************************************************************************* //
