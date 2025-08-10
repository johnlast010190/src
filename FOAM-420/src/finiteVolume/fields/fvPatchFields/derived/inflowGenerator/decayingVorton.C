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
    (c) 2017-2022 Esi Ltd.
    (c) 2014 Hannes Kroeger

\*---------------------------------------------------------------------------*/


#include "fields/fvPatchFields/derived/inflowGenerator/decayingVorton.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "primitives/random/random11/random11.H"
extern Foam::random11 ranGen;
#else
#include "primitives/random/Random/Random.H"
extern Foam::Random ranGen;
#endif

namespace Foam
{

decayingVorton::decayingVorton()
:
    length_  (0.0),
    location_(pTraits<vector>::zero),
    omega_   (pTraits<vector>::zero),
    velocity_(pTraits<vector>::zero),
    xmax_    (0.0)
{
}

decayingVorton::decayingVorton(Istream& s)
:
    length_  (readScalar(s)),
    location_(s),
    omega_   (s),
    velocity_(s),
    xmax_    (readScalar(s))
{
}

decayingVorton::decayingVorton(scalar length, const vector& location, const vector& velocity, scalar xmax)
:
    length_  (length),
    location_(location),
    omega_   (pTraits<vector>::zero),
    velocity_(velocity),
    xmax_    (xmax)
{
    vector omg = 2*ranGen.vector01() - pTraits<vector>::one;
    omega_     = omg/mag(omg);
}

autoPtr<decayingVorton> decayingVorton::New(Istream& s)
{
    Info<<"reading"<<endl;
    return autoPtr<decayingVorton>(new decayingVorton(s));
}


bool decayingVorton::operator!=(const decayingVorton& vt) const
{
    return (length_   != vt.length_)   ||
           (location_ != vt.location_) ||
           (omega_    != vt.omega_)    ||
           (velocity_ != vt.velocity_) ||
           (xmax_     != vt.xmax_);
}

Ostream& operator<<(Ostream& s, const decayingVorton& vt)
{
    s<<vt.length_<<endl;
    s<<vt.location_<<endl;
    s<<vt.omega_<<endl;
    s<<vt.velocity_<<endl;
    s<<vt.xmax_<<endl;

    return s;
}

Istream& operator>>(Istream& s, decayingVorton& vt)
{
    scalar len(readScalar(s));
    vector loc(s);
    vector omg(s);
    vector vel(s);
    scalar xmax(readScalar(s));

    vt.length_   = len;
    vt.location_ = loc;
    vt.omega_    = omg;
    vt.velocity_ = vel;
    vt.xmax_     = xmax;

    return s;
}

}
