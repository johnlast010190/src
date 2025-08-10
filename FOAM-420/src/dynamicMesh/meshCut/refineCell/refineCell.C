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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshCut/refineCell/refineCell.H"
#include "db/error/error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refineCell::refineCell()
:
    cellNo_(-1),
    direction_(Zero)
{}


Foam::refineCell::refineCell(const label celli, const vector& direction)
:
    cellNo_(celli),
    direction_(direction)
{
    scalar magDir = mag(direction_);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "(almost)zero vector as direction for cell " << cellNo_
            << abort(FatalError);
    }
    else if (mag(magDir - 1) > SMALL)
    {
        // Normalize
        direction_ /= mag(direction_);
    }
}


Foam::refineCell::refineCell(Istream& is)
:
    cellNo_(readLabel(is)),
    direction_(is)
{
    scalar magDir = mag(direction_);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "(almost)zero vector as direction for cell " << cellNo_
            << abort(FatalError);
    }
    else if (mag(magDir - 1) > SMALL)
    {
        // Normalize
        direction_ /= mag(direction_);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const refineCell& r)
{
    if (os.format() == IOstream::ASCII)
    {
        os << r.cellNo() << token::SPACE << r.direction();
    }
    else
    {
        os << r.cellNo() << r.direction();
    }
    return os;
}


// ************************************************************************* //
