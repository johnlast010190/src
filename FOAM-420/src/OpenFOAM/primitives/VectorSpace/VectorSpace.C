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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/VectorSpace/VectorSpace.H"
#include "db/IOstreams/IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::VectorSpace<Form, Cmpt, Ncmpts>::VectorSpace
(
    Istream& is
)
{
    // Read beginning of VectorSpace<Cmpt>
    is.readBegin("VectorSpace<Form, Cmpt, Ncmpts>");

    for (direction i=0; i<Ncmpts; i++)
    {
        is >> v_[i];
    }

    // Read end of VectorSpace<Cmpt>
    is.readEnd("VectorSpace<Form, Cmpt, Ncmpts>");

    // Check state of Istream
    is.check(FUNCTION_NAME);
}


template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::word Foam::name
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    std::ostringstream buf;

    buf << '(' << vs.v_[0];

    for (direction i=1; i<Ncmpts; i++)
    {
        buf << ',' << vs.v_[i];
    }

    buf << ')';

    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    // Read beginning of VectorSpace<Cmpt, Ncmpts>
    is.readBegin("VectorSpace<Form, Cmpt, Ncmpts>");

    for (direction i=0; i<Ncmpts; i++)
    {
        is >> vs.v_[i];
    }

    // Read end of VectorSpace<Cmpt, Ncmpts>
    is.readEnd("VectorSpace<Form, Cmpt, Ncmpts>");

    // Check state of Istream
    is.check(FUNCTION_NAME);

    return is;
}


template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    os << token::BEGIN_LIST << token::SPACE << vs.v_[0];

    for (direction i=1; i<Ncmpts; i++)
    {
        os << token::SPACE << vs.v_[i];
    }

    os << token::SPACE << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
