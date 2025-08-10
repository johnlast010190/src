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

#include "db/dictionary/entry/entry.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/StringStreams/StringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::entry::disableFunctionEntries
(
    Foam::debug::infoSwitch("disableFunctionEntries", 0)
);


Foam::entry::inputMode Foam::entry::globalInputMode = inputMode::MERGE;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::entry::resetInputMode()
{
    globalInputMode = inputMode::MERGE;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entry::entry(const keyType& keyword)
:
    IDLList<entry>::link(),
    keyword_(keyword)
{}


Foam::entry::entry(const entry& e)
:
    IDLList<entry>::link(),
    keyword_(e.keyword_)
{}


Foam::autoPtr<Foam::entry> Foam::entry::clone() const
{
    return clone(dictionary::null);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::entry::operator=(const entry& e)
{
    // check for assignment to self
    if (this == &e)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    keyword_ = e.keyword_;
}


bool Foam::entry::operator==(const entry& e) const
{
    if (keyword_ != e.keyword_)
    {
        return false;
    }
    else
    {
        OStringStream oss1;
        oss1 << *this;

        OStringStream oss2;
        oss2 << e;

        return oss1.str() == oss2.str();
    }
}


bool Foam::entry::operator!=(const entry& e) const
{
    return !operator==(e);
}


// ************************************************************************* //
