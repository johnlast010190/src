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
    (c) 2011 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::PairCollisionRecord<Type>::match
(
    label queryOrigProcOfOther,
    label queryOrigIdOfOther
) const
{
    return
    (
        queryOrigProcOfOther == origProcOfOther()
     && queryOrigIdOfOther == origIdOfOther()
    );
}


template<class Type>
inline Foam::label Foam::PairCollisionRecord<Type>::origProcOfOther() const
{
    return mag(origProcOfOther_) - 1;
}


template<class Type>
inline Foam::label Foam::PairCollisionRecord<Type>::origIdOfOther() const
{
    return origIdOfOther_;
}


template<class Type>
inline const Type&
Foam::PairCollisionRecord<Type>::collisionData() const
{
    return data_;
}


template<class Type>
inline Type& Foam::PairCollisionRecord<Type>::collisionData()
{
    return data_;
}


template<class Type>
inline bool Foam::PairCollisionRecord<Type>::accessed() const
{
    return pos0(origProcOfOther_);
}


template<class Type>
inline void Foam::PairCollisionRecord<Type>::setAccessed()
{
    origProcOfOther_ = origProcOfOther() + 1;
}


template<class Type>
inline void Foam::PairCollisionRecord<Type>::setUnaccessed()
{
    origProcOfOther_ = -(origProcOfOther() + 1);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::operator==
(
    const PairCollisionRecord<Type>& a,
    const PairCollisionRecord<Type>& b
)
{
    return
    (
        a.origProcOfOther_ == b.origProcOfOther_
     && a.origIdOfOther_ == b.origIdOfOther_
     && a.data_ == b.data_
    );
}


template<class Type>
inline bool Foam::operator!=
(
    const PairCollisionRecord<Type>& a,
    const PairCollisionRecord<Type>& b
)
{
    return !(a == b);
}


// ************************************************************************* //
