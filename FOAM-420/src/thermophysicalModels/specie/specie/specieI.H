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

#include "specie/specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline specie::specie
(
    const fvMesh& mesh,
    const word& name,
    const scalar Y,
    const scalar molWeight
)
:
    mesh_(mesh),
    name_(name),
    Y_(Y),
    molWeight_(molWeight)
{}


inline specie::specie
(
    const fvMesh& mesh,
    const scalar Y,
    const scalar molWeight
)
:
    mesh_(mesh),
    Y_(Y),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline specie::specie(const word& name, const specie& st)
:
    mesh_(st.mesh_),
    name_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& specie::name() const
{
    return name_;
}


inline scalar specie::W() const
{
    return molWeight_;
}


inline scalar specie::Y() const
{
    return Y_;
}


inline scalar specie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void specie::operator=(const specie& st)
{
    //mesh_ = st.mesh_;
    //name_ = st.name_;
    Y_ = st.Y_;
    molWeight_ = st.molWeight_;
}


inline void specie::operator+=(const specie& st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > SMALL)
    {
        molWeight_ = sumY/(Y_/molWeight_ + st.Y_/st.molWeight_);
    }

    Y_ = sumY;
}


inline void specie::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline specie operator+(const specie& st1, const specie& st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, SMALL);

    if (mag(sumY) > SMALL)
    {
        return specie
        (
            st1.mesh_,
            sumY,
            sumY/(st1.Y_/st1.molWeight_ + st2.Y_/st2.molWeight_)
        );
    }
    else
    {
        return st1;
    }
}


inline specie operator*(const scalar s, const specie& st)
{
    return specie
    (
        st.mesh_,
        s*st.Y_,
        st.molWeight_
    );
}


inline specie operator==(const specie& st1, const specie& st2)
{
    scalar diffY = st2.Y_ - st1.Y_;
    if (mag(diffY) < SMALL)
    {
        diffY = SMALL;
    }

    const scalar diffRW = st2.Y_/st2.molWeight_ - st1.Y_/st1.molWeight_;

    // if (mag(diffRW) > SMALL)
    // {
    //     molWeight = diffY/diffRW;
    // }

    // Using intermediate volatile bool to prevent compiler optimising out the
    // if block (above) - CLANG 3.7.1
    volatile const bool valid = (mag(diffRW) > SMALL);
    const scalar molWeight = valid ? diffY/diffRW : GREAT;

    return specie(st1.mesh_, diffY, molWeight);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
