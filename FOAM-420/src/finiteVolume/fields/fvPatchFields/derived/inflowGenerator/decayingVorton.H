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

#pragma once

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

namespace Foam
{

class decayingVorton
{
public:

    decayingVorton();

    decayingVorton(Istream& s);

    decayingVorton(scalar length, const vector& location, const vector& velocity, scalar xmax);

    static autoPtr<decayingVorton> New(Istream& s);

    const vector& location() const { return location_; }

    const vector& velocity() const { return velocity_; }

    bool operator!=(const decayingVorton& vt) const;

    autoPtr<decayingVorton> clone() const;

    scalar length() const { return length_; }

    scalar lspot() const { return 3*length_; }

    scalar xmax() const { return xmax_; }

    vector velocityAt(const vector& v) const
    {
        vector dv = v - location_;
        scalar nrm2 = magSqr(dv);
        vector t = dv ^ omega_;

        return (1/length_)*exp(-(M_PI/2)*nrm2/(length_*length_))*t;
    }

    void move(scalar dt) { location_ += vector(velocity_.x(), 0, 0)*dt; }


    friend Ostream& operator<<(Ostream& s, const decayingVorton& vt);

    friend Istream& operator>>(Istream& s, decayingVorton& vt);

private:

    scalar length_;

    vector location_;

    vector omega_;

    vector velocity_;

    scalar xmax_;

};

}
