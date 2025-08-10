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

Description

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline objectMap::objectMap()
:
    index_(-1),
    masterObjects_(0)
{}


inline objectMap::objectMap(const label index, const labelList& master)
:
    index_(index),
    masterObjects_(master)
{}


inline objectMap::objectMap(Istream& is)
{
    // Read beginning of objectMap
    is.readBegin("objectMap");

    is >> index_ >> static_cast<labelList&>(masterObjects_);

    // Read master of objectMap
    is.readEnd("objectMap");

    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label& objectMap::index()
{
    return index_;
}


inline label objectMap::index() const
{
    return index_;
}


inline labelList& objectMap::masterObjects()
{
    return masterObjects_;
}


inline const labelList& objectMap::masterObjects() const
{
    return masterObjects_;
}


// * * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * //

inline bool operator==(const objectMap& a, const objectMap& b)
{
    return
    (
        (a.index_ == b.index_) && (a.masterObjects_ == b.masterObjects_)
    );
}


inline bool operator!=(const objectMap& a, const objectMap& b)
{
    return (!(a == b));
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const objectMap& a)
{
    os  << token::BEGIN_LIST
        << a.index_ << token::SPACE
        << a.masterObjects_
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


inline Istream& operator>>(Istream& is, objectMap& a)
{
    is.readBegin("objectMap");
    is  >> a.index_ >> a.masterObjects_;
    is.readEnd("objectMap");

    is.check(FUNCTION_NAME);
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // Master namespace Foam

// ************************************************************************* //
