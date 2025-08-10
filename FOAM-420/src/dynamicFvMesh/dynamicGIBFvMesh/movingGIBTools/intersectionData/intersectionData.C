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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/movingGIBTools/intersectionData/intersectionData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::intersectionData::intersectionData
(
    const label& pp,
    const label& pI,
    const label& pAdd,
    const label& fromProc,
    const label& toPatch,
    const label& neiProcPatch,
    const point& pS,
    const point& pE,
    const bool& hitWall,
    const bool& hitProc,
    const bool& hitNone,
    const point& hitP
)
:
    pp_(pp),
    pI_(pI),
    pAdd_(pAdd),
    fromProc_(fromProc),
    toProc_(toPatch),
    neiProcPatch_(neiProcPatch),
    pS_(pS),
    pE_(pE),
    hitWall_(hitWall),
    hitProc_(hitProc),
    hitNone_(hitNone),
    hitP_(hitP)
{
}

Foam::intersectionData::intersectionData()
:
    pp_(-1),
    pI_(-1),
    pAdd_(-1),
    fromProc_(-1),
    toProc_(-1),
    neiProcPatch_(-1),
    pS_(vector::zero),
    pE_(vector::zero),
    hitWall_(false),
    hitProc_(false),
    hitNone_(false),
    hitP_()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::intersectionData::operator==
(
    const Foam::intersectionData& rhs
) const
{
    return
    (
        pp() == rhs.pp()
     && pI() == rhs.pI()
     && pAdd() == rhs.pAdd()
     && fromProc() == rhs.fromProc()
     && toProc() == rhs.toProc()
     && neiProcPatch() == rhs.neiProcPatch()
     && pS() == rhs.pS()
     && pE() == rhs.pE()
     && hitWall() == rhs.hitWall()
     && hitProc() == rhs.hitProc()
     && hitNone() == rhs.hitNone()
     && hitP() == rhs.hitP()
    );
}


bool Foam::intersectionData::operator!=
(
    const Foam::intersectionData& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * *


Foam::Istream& Foam::operator>>(Istream& is, intersectionData& iD)
{
    is  >> iD.pp()
        >> iD.pI()
        >> iD.pAdd()
        >> iD.fromProc()
        >> iD.toProc()
        >> iD.neiProcPatch()
        >> iD.pS()
        >> iD.pE()
        >> iD.hitWall()
        >> iD.hitProc()
        >> iD.hitNone()
        >> iD.hitP();

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>(Foam::Istream&, "
        "Foam::intersectionData&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const intersectionData& iD)
{
    os  << iD.pp() << token::SPACE
        << iD.pI() << token::SPACE
        << iD.pAdd() << token::SPACE
        << iD.fromProc() << token::SPACE
        << iD.toProc() << token::SPACE
        << iD.neiProcPatch() << token::SPACE
        << iD.pS() << token::SPACE
        << iD.pE() << token::SPACE
        << iD.hitWall() << token::SPACE
        << iD.hitProc() << token::SPACE
        << iD.hitNone() << token::SPACE
        << iD.hitP() << token::SPACE
        << endl;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::intersectionData&)"
    );

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
