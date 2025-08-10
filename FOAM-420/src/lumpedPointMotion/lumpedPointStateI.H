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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

inline bool Foam::lumpedPointState::valid() const
{
    return points_.size() && points_.size() == angles_.size();
}


inline bool Foam::lumpedPointState::empty() const
{
    return points_.empty();
}


inline Foam::label Foam::lumpedPointState::size() const
{
    return points_.size();
}


inline const Foam::pointField& Foam::lumpedPointState::points() const
{
    return points_;
}


inline const Foam::vectorField& Foam::lumpedPointState::angles() const
{
    return angles_;
}


inline const Foam::tensorField& Foam::lumpedPointState::rotations() const
{
    if (!rotationPtr_)
    {
        calcRotations();
    }

    return *rotationPtr_;
}


// ************************************************************************* //
