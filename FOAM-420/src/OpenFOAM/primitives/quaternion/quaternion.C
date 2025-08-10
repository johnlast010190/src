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
    (c) 2019-2020 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/quaternion/quaternion.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOstreams/StringStreams/OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::quaternion::typeName = "quaternion";
const Foam::quaternion Foam::quaternion::zero(0, vector(0, 0, 0));
const Foam::quaternion Foam::quaternion::I(1, vector(0, 0, 0));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quaternion::quaternion(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const quaternion& q)
{
    OStringStream buf;
    buf << '(' << q.w() << ',' << q.v() << ')';
    return buf.str();
}


Foam::quaternion Foam::slerp
(
    const quaternion& qa,
    const quaternion& qb,
    const scalar t
)
{
    label sign = 1;

    if ((qa & qb) < 0)
    {
        sign = -1;
    }

    return qa*pow((inv(qa)*sign*qb), t);
}


Foam::quaternion Foam::average
(
    const UList<quaternion>& qs,
    const UList<scalar> w
)
{
    quaternion qa(w[0]*qs[0]);

    for (label i=1; i<qs.size(); i++)
    {
        // Invert quaternion if it has the opposite sign to the average
        if ((qa & qs[i]) > 0)
        {
            qa += w[i]*qs[i];
        }
        else
        {
            qa -= w[i]*qs[i];
        }
    }

    return qa;
}


Foam::quaternion Foam::exp(const quaternion& q)
{
    const scalar magV = mag(q.v());

    if (magV == 0)
    {
        return quaternion(1, Zero);
    }

    const scalar expW = exp(q.w());

    return quaternion
    (
        expW*cos(magV),
        expW*sin(magV)*q.v()/magV
    );
}


Foam::quaternion Foam::pow(const quaternion& q, const label power)
{
    const scalar magQ = mag(q);
    const scalar magV = mag(q.v());

    quaternion powq(q.v());

    if (magV != 0 && magQ != 0)
    {
        powq /= magV;
        powq *= power*acos(q.w()/magQ);
    }

    return pow(magQ, power)*exp(powq);
}


Foam::quaternion Foam::pow(const quaternion& q, const scalar power)
{
    const scalar magQ = mag(q);
    const scalar magV = mag(q.v());

    quaternion powq(q.v());

    if (magV != 0 && magQ != 0)
    {
        powq /= magV;
        powq *= power*acos(q.w()/magQ);
    }

    return pow(magQ, power)*exp(powq);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, quaternion& q)
{
    // Read beginning of quaternion
    is.readBegin("quaternion");

    is  >> q.w() >> q.v();

    // Read end of quaternion
    is.readEnd("quaternion");

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const quaternion& q)
{
    os  << token::BEGIN_LIST
        << q.w() << token::SPACE << q.v()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
