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

\*---------------------------------------------------------------------------*/

#include "parcels/Templates/CollidingParcel/CollisionRecordList/PairCollisionRecord/PairCollisionRecord.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord()
:
    origProcOfOther_(0),
    origIdOfOther_(-1),
    data_(Zero)
{}


template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord
(
    bool accessed,
    label origProcOfOther,
    label origIdOfOther,
    const Type& data
)
:
    origProcOfOther_(origProcOfOther + 1),
    origIdOfOther_(origIdOfOther),
    data_(data)
{
    // Default assignment to origProcOfOther_ assumes accessed is true

    if (!accessed)
    {
        setUnaccessed();
    }
}


template<class Type>
Foam::PairCollisionRecord<Type>::PairCollisionRecord
(
    const PairCollisionRecord<Type>& pCR
)
:
    origProcOfOther_(pCR.origProcOfOther_),
    origIdOfOther_(pCR.origIdOfOther_),
    data_(pCR.data_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PairCollisionRecord<Type>::~PairCollisionRecord()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::PairCollisionRecord<Type>::operator=
(
    const PairCollisionRecord<Type>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    origProcOfOther_ = rhs.origProcOfOther_;
    origIdOfOther_ = rhs.origIdOfOther_;
    data_ = rhs.data_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "parcels/Templates/CollidingParcel/CollisionRecordList/PairCollisionRecord/PairCollisionRecordIO.C"


// ************************************************************************* //
