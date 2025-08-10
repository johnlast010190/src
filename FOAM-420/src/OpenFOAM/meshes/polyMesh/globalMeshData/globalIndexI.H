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

#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex()
{}


Foam::globalIndex::globalIndex(const Xfer<labelList>& offsets)
:
    offsets_(offsets)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::globalIndex::empty() const
{
    return offsets_.empty() || offsets_.last() == 0;
}

inline Foam::labelList& Foam::globalIndex::offsets()
{
    return offsets_;
}


inline Foam::label Foam::globalIndex::offset(const label proci) const
{
    return offsets_[proci];
}


inline Foam::label Foam::globalIndex::localSize(const label proci) const
{
    return offsets_[proci+1] - offsets_[proci];
}


inline Foam::label Foam::globalIndex::localSize() const
{
    return localSize(Pstream::myProcNo());
}


inline Foam::label Foam::globalIndex::size() const
{
    return offsets_.last();
}


inline Foam::label Foam::globalIndex::toGlobal
(
    const label proci,
    const label i
) const
{
    return i + offsets_[proci];
}


inline Foam::label Foam::globalIndex::toGlobal(const label i) const
{
    return toGlobal(Pstream::myProcNo(), i);
}


//- Is on local processor
inline bool Foam::globalIndex::isLocal(const label proci, const label i) const
{
    return i >= offsets_[proci] && i < offsets_[proci+1];
}


inline bool Foam::globalIndex::isLocal(const label i) const
{
    return isLocal(Pstream::myProcNo(), i);
}


inline Foam::label Foam::globalIndex::toLocal(const label proci, const label i)
const
{
    label localI = i - offsets_[proci];

    if (localI < 0 || i >= offsets_[proci+1])
    {
        FatalErrorInFunction
            << "Global " << i << " does not belong on processor "
            << proci << endl << "Offsets:" << offsets_
            << abort(FatalError);
    }
    return localI;
}


inline Foam::label Foam::globalIndex::toLocal(const label i) const
{
    return toLocal(Pstream::myProcNo(), i);
}


inline void Foam::globalIndex::inplaceToGlobal
(
    const label proci,
    labelList& labels
) const
{
    const label off = offsets_[proci];

    for (label& val : labels)
    {
        val += off;
    }
}


inline void Foam::globalIndex::inplaceToGlobal(labelList& labels) const
{
    inplaceToGlobal(Pstream::myProcNo(), labels);
}


inline Foam::label Foam::globalIndex::whichProcID(const label i) const
{
    if (i < 0 || i >= size())
    {
        FatalErrorInFunction
            << "Global " << i << " does not belong on any processor."
            << " Offsets:" << offsets_
            << abort(FatalError);
    }

    return findLower(offsets_, i+1);
}


// ************************************************************************* //
