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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "thermo/hFunctionT/hFunctionTThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState>::hFunctionTThermo
(
    const EquationOfState& pt,
    const scalar Hf,
    const scalar Sf,
    const Function1<scalar>& Cp
)
:
    EquationOfState(pt),
    Hf_(Hf),
    Sf_(Sf),
    Cp_(Cp.clone().ptr())
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState>::hFunctionTThermo
(
    const word& name,
    const hFunctionTThermo& pt
)
:
    EquationOfState(name, pt),
    Hf_(pt.Hf_),
    Sf_(pt.Sf_),
    Cp_(pt.Cp_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::Cp
(
    const scalar p, const scalar T
) const
{
    return Cp_().value(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Cp_().integrate(0, T) + Hf_ + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::Hc()
const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hFunctionTThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return Cp_().integrateYoverX(Tstd, T) + Sf_ + EquationOfState::S(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hFunctionTThermo<EquationOfState>::operator=
(
    const hFunctionTThermo<EquationOfState>& pt
)
{
    EquationOfState::operator=(pt);

    Hf_ = pt.Hf_;
    Sf_ = pt.Sf_;
    Cp_.reset(pt.Cp_.clone());
}


template<class EquationOfState>
inline void Foam::hFunctionTThermo<EquationOfState>::operator+=
(
    const hFunctionTThermo<EquationOfState>& pt
)
{
    EquationOfState::operator+=(pt);

    if (mag(this->Y()) > SMALL)
    {
        NotImplemented;
    }
}


template<class EquationOfState>
inline void Foam::hFunctionTThermo<EquationOfState>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState> Foam::operator+
(
    const hFunctionTThermo<EquationOfState>& pt1,
    const hFunctionTThermo<EquationOfState>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return hFunctionTThermo<EquationOfState>
        (
            eofs,
            pt1.Hf_,
            pt1.Sf_,
            pt1.Cp_(),
            pt1.h_(),
            pt1.s_()
        );
    }
    {
        NotImplemented;
    }
}


template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hFunctionTThermo<EquationOfState>& pt
)
{
    return hFunctionTThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.Hf_,
        pt.Sf_,
        pt.Cp_(),
        pt.h_(),
        pt.s_()
    );
}


template<class EquationOfState>
inline Foam::hFunctionTThermo<EquationOfState> Foam::operator==
(
    const hFunctionTThermo<EquationOfState>& pt1,
    const hFunctionTThermo<EquationOfState>& pt2
)
{
    NotImplemented;
}


// ************************************************************************* //
