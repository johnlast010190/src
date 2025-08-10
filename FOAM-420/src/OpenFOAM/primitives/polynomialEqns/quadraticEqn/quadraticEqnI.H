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
    (c) 2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::quadraticEqn::quadraticEqn()
{}


inline Foam::quadraticEqn::quadraticEqn(const Foam::zero)
:
    VectorSpace<quadraticEqn, scalar, 3>(Foam::zero())
{}


inline Foam::quadraticEqn::quadraticEqn
(
    const scalar a,
    const scalar b,
    const scalar c
)
{
    this->v_[A] = a;
    this->v_[B] = b;
    this->v_[C] = c;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::quadraticEqn::a() const
{
    return this->v_[A];
}


inline Foam::scalar Foam::quadraticEqn::b() const
{
    return this->v_[B];
}


inline Foam::scalar Foam::quadraticEqn::c() const
{
    return this->v_[C];
}


inline Foam::scalar& Foam::quadraticEqn::a()
{
    return this->v_[A];
}


inline Foam::scalar& Foam::quadraticEqn::b()
{
    return this->v_[B];
}


inline Foam::scalar& Foam::quadraticEqn::c()
{
    return this->v_[C];
}


inline Foam::scalar Foam::quadraticEqn::value(const scalar x) const
{
    return x*(x*a() + b()) + c();
}


inline Foam::scalar Foam::quadraticEqn::error(const scalar x) const
{
    return mag(SMALL*x*(x*2*a() + b()));
}


// ************************************************************************* //
