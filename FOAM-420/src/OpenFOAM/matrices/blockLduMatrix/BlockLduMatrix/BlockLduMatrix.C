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
    (c) 2004-6 H. Jasak All rights reserved
    (c) 2019 Esi Ltd.

Description
    BlockLduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.  Addressing object must
    be supplied for the upper and lower triangles.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"
#include "db/IOstreams/IOstreams.H"
#include "include/demandDrivenData.H"
#include "db/IOobjects/IOField/IOField.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::label Foam::BlockLduMatrix<Type>::fixFillIn
(
    debug::optimisationSwitch("matrixConstraintFillIn", 4)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix(const lduMesh& ldu)
:
    lduMesh_(ldu),
    diagPtr_(nullptr),
    upperPtr_(nullptr),
    lowerPtr_(nullptr),
    interfaces_(ldu.interfaces().size()),
    coupleUpper_(ldu.interfaces().size()),
    coupleLower_(ldu.interfaces().size()),
    fixedEqns_(ldu.lduAddr().size()/fixFillIn)
{
    const lduAddressing& addr = ldu.lduAddr();

    forAll(coupleUpper_, i)
    {
        if (ldu.interfaces().set(i))
        {
            coupleUpper_.set(i, new CoeffField<Type>(addr.patchAddr(i).size()));
            coupleLower_.set(i, new CoeffField<Type>(addr.patchAddr(i).size()));
        }
        else
        {
            coupleUpper_.set(i, new CoeffField<Type>(0));
            coupleLower_.set(i, new CoeffField<Type>(0));
        }
    }
}


template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix(const BlockLduMatrix<Type>& A)
:
    refCount(),
    lduMesh_(A.lduMesh_),
    diagPtr_(nullptr),
    upperPtr_(nullptr),
    lowerPtr_(nullptr),
    interfaces_(A.interfaces()),
    coupleUpper_(A.interfaces().size()),
    coupleLower_(A.interfaces().size()),
    fixedEqns_(A.fixedEqns_)
{
    if (A.diagPtr_)
    {
        diagPtr_ = new TypeCoeffField(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = new TypeCoeffField(*(A.upperPtr_));
    }

    if (A.lowerPtr_)
    {
        lowerPtr_ = new TypeCoeffField(*(A.lowerPtr_));
    }

  forAll(A.lduMesh_.interfaces(), i)
    {
        if (A.lduMesh_.interfaces().set(i))
        {
            coupleUpper_.set(i, A.coupleUpper_[i]);
            coupleLower_.set(i, A.coupleUpper_[i]);
        }
        else
        {
            coupleUpper_.set(i, new CoeffField<Type>(0));
            coupleLower_.set(i, new CoeffField<Type>(0));
        }
    }
}


//HJ, problematic: memmory management.
// Reconsider.  HJ, 7/Nov/2010
template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix(BlockLduMatrix<Type>& A, bool reUse)
:
    refCount(),
    lduMesh_(A.lduMesh_),
    diagPtr_(nullptr),
    upperPtr_(nullptr),
    lowerPtr_(nullptr),
    interfaces_(A.interfaces_, reUse),
    coupleUpper_(A.coupleUpper_, reUse),
    coupleLower_(A.coupleLower_, reUse),
    fixedEqns_(A.fixedEqns_)
{
    if (reUse)
    {
        if (A.lowerPtr_)
        {
            lowerPtr_ = A.lowerPtr_;
            A.lowerPtr_ = nullptr;
        }

        if (A.diagPtr_)
        {
            diagPtr_ = A.diagPtr_;
            A.diagPtr_ = nullptr;
        }

        if (A.upperPtr_)
        {
            upperPtr_ = A.upperPtr_;
            A.upperPtr_ = nullptr;
        }
    }
    else
    {
        if (A.diagPtr_)
        {
            diagPtr_ = new TypeCoeffField(*(A.diagPtr_));
        }

        if (A.upperPtr_)
        {
            upperPtr_ = new TypeCoeffField(*(A.upperPtr_));
        }

        if (A.lowerPtr_)
        {
            lowerPtr_ = new TypeCoeffField(*(A.lowerPtr_));
        }
    }
}


template<class Type>
Foam::BlockLduMatrix<Type>::BlockLduMatrix
(
    const lduMesh& mesh,
    Istream& is
)
:
    lduMesh_(mesh),
    diagPtr_(nullptr),
    upperPtr_(nullptr),
    lowerPtr_(nullptr),
    interfaces_(mesh.interfaces().size()),
    coupleUpper_(mesh.interfaces().size()),
    coupleLower_(mesh.interfaces().size()),
    fixedEqns_(mesh.lduAddr().size()/fixFillIn)
{
//    Switch hasDiag(is);
//    Switch hasUp(is);
//    Switch hasLow(is);

//    if (hasDiag)
    {
        diagPtr_ = new TypeCoeffField(is);
    }
//    if (hasUp)
    {
        upperPtr_ = new TypeCoeffField(is);
    }
//    if (hasLow)
    {
        lowerPtr_ = new TypeCoeffField(is);
    }
    forAll(mesh.interfaces(), i)
    {
        if (mesh.interfaces().set(i))
        {
            coupleUpper_.set(i, new TypeCoeffField(is));
            coupleLower_.set(i, new TypeCoeffField(is));
        }
        else
        {
            coupleUpper_.set(i, new CoeffField<Type>(0));
            coupleLower_.set(i, new CoeffField<Type>(0));
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockLduMatrix<Type>::~BlockLduMatrix()
{
    deleteDemandDrivenData(diagPtr_);
    deleteDemandDrivenData(upperPtr_);
    deleteDemandDrivenData(lowerPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ = new TypeCoeffField(lduAddr().size());
    }

    return *diagPtr_;
}


template<class Type>
const typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorInFunction
            << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


template<class Type>
typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::upper()
{
    if (!upperPtr_)
    {
        upperPtr_ = new TypeCoeffField(lduAddr().lowerAddr().size());
    }

    return *upperPtr_;
}


template<class Type>
const typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::upper() const
{
    if (!upperPtr_)
    {
        FatalErrorInFunction
            << "upperPtr_ unallocated"
            << abort(FatalError);
    }

    return *upperPtr_;
}


template<class Type>
typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::lower()
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
            lowerPtr_ = new TypeCoeffField(upperPtr_->transpose());
        }
        else
        {
            lowerPtr_ = new TypeCoeffField(lduAddr().lowerAddr().size());
        }
    }

    return *lowerPtr_;
}


template<class Type>
const typename Foam::BlockLduMatrix<Type>::TypeCoeffField&
Foam::BlockLduMatrix<Type>::lower() const
{
    if (!lowerPtr_)
    {
        FatalErrorInFunction
            << "lowerPtr_ unallocated"
            << abort(FatalError);
    }

    return *lowerPtr_;
}


template<class Type>
void Foam::BlockLduMatrix<Type>::clearInterfaces()
{
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            delete interfaces_(i);
        }
    }
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::empty() const
{
    return (!diagPtr_ && !lowerPtr_ && !upperPtr_);
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::diagonal() const
{
    return (diagPtr_ && !lowerPtr_ && !upperPtr_);
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::symmetric() const
{
    if (lowerPtr_ && !upperPtr_)
    {
        FatalErrorInFunction
            << "Matrix assembly error: symmetric matrix but only lower "
            << "triangle is allocated.  This is not allowed."
            << abort(FatalError);
    }

    return (diagPtr_ && (!lowerPtr_ && upperPtr_));
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::asymmetric() const
{
    return (diagPtr_ && lowerPtr_ && upperPtr_);
}


template<class Type>
bool Foam::BlockLduMatrix<Type>::componentCoupled() const
{
    // Return true if the matrix coefficient couple the components
    if (thereIsDiag())
    {
        if (diag().activeType() == blockCoeffBase::SQUARE)
        {
            return true;
        }
    }

    if (thereIsUpper())
    {
        if (upper().activeType() == blockCoeffBase::SQUARE)
        {
            return true;
        }
    }

    if (thereIsLower())
    {
        if (lower().activeType() == blockCoeffBase::SQUARE)
        {
            return true;
        }
    }

    return false;
}


template<class Type>
void Foam::BlockLduMatrix<Type>::setResidualField
(
    const Field<Type>& residual,
    const word& fieldName
) const
{
    if (!lduMesh_.hasDb())
    {
        return;
    }

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        const word resultName = word("initialResidual_")
            + fieldName + std::to_string(cmpt);
        IOField<scalar>* residualPtr =
            this->mesh().thisDb().template lookupObjectRefPtr<IOField<scalar>>
            (
                resultName
            );

        if (residualPtr)
        {
            const IOdictionary* dataPtr =
                this->mesh().thisDb().template lookupObjectPtr<IOdictionary>
                (
                    "data"
                );

            if (dataPtr)
            {
                if (dataPtr->found("firstIteration"))
                {
                    *residualPtr = residual.component(cmpt);
                    DebugInformation
                        << "Setting residual field for first solver iteration "
                        << "for solver field: " << fieldName << endl;
                }
            }
            else
            {
                *residualPtr = residual.component(cmpt);
                DebugInformation
                    << "Setting residual field for solver field "
                    << fieldName << endl;
            }
        }

    }
}




// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const BlockLduMatrix<Type>& ldum)
{

    Switch hasDiag = ldum.thereIsDiag();
    Switch hasUp = ldum.thereIsUpper();
    Switch hasLow = ldum.thereIsLower();

    if (hasDiag)
    {
        os  << *ldum.diagPtr_;
    }

    if (hasUp)
    {
        os  << *ldum.upperPtr_;
    }

    if (hasLow)
    {
        os  << *ldum.lowerPtr_;
    }


    forAll(ldum.mesh().interfaces(), i)
    {
        if (ldum.mesh().interfaces().set(i))
        {
            os << ldum.coupleUpper_[i];
            os << ldum.coupleLower_[i];
        }
    }

    os.check(FUNCTION_NAME);

    return os;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<BlockLduMatrix<Type>>& ip
)
{
    const BlockLduMatrix<Type>& ldum = ip.t_;
    if (ldum.diagPtr_)
    {
        os.writeEntry("diagonal", *ldum.diagPtr_);
    }
    else
    {
        // Dummy write for consistency
        os.writeEntry
        (
            "diagonal",
            typename BlockLduMatrix<Type>::TypeCoeffField
            (
                ldum.lduAddr().size()
            )
        );
    }

    if (ldum.upperPtr_)
    {
        os.writeEntry("upper", *ldum.upperPtr_);
    }
    else
    {
        // Dummy write for consistency
        os.writeEntry
        (
            "upper",
            typename BlockLduMatrix<Type>::TypeCoeffField
            (
                ldum.lduAddr().lowerAddr().size()
            )
        );
    }

    if (ldum.lowerPtr_)
    {
        os.writeEntry("lower", *ldum.lowerPtr_);
    }
    else
    {
        // Dummy write for consistency
        os.writeEntry
        (
            "lower",
            typename BlockLduMatrix<Type>::TypeCoeffField
            (
                ldum.lduAddr().lowerAddr().size()
            )
        );
    }

    os.writeEntry("coupleUpper", ldum.coupleUpper_);
    os.writeEntry("coupleLower", ldum.coupleLower_);

    os.check("Ostream& operator<<(Ostream&, const BlockLduMatrix<Type>&");

    return os;
}




// ************************************************************************* //
