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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::fv::cellSetOption::timeStart() const
{
    return timeStart_;
}


inline Foam::scalar Foam::fv::cellSetOption::duration() const
{
    return duration_;
}


inline bool Foam::fv::cellSetOption::inTimeLimits(const scalar time) const
{
    return
    (
        (timeStart_ < 0)
     ||
        (
            (mesh_.time().value() >= timeStart_)
         && (mesh_.time().value() <= (timeStart_ + duration_))
        )
    );
}


inline const Foam::fv::cellSetOption::selectionModeType&
Foam::fv::cellSetOption::selectionMode() const
{
    return selectionMode_;
}


inline const Foam::word& Foam::fv::cellSetOption::cellSetName() const
{
    return cellSetName_;
}


inline Foam::scalar Foam::fv::cellSetOption::V() const
{
    return V_;
}


inline const Foam::labelList& Foam::fv::cellSetOption::cells() const
{
    return cells_;
}


inline Foam::scalar& Foam::fv::cellSetOption::timeStart()
{
    return timeStart_;
}


inline Foam::scalar& Foam::fv::cellSetOption::duration()
{
    return duration_;
}


// ************************************************************************* //
