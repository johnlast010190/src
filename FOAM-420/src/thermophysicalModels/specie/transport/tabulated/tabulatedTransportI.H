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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "specie/specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const Thermo& t,
    const dictionary&  muDict,
    const dictionary&  kappaDict,
    const autoPtr<interpolation2DTable<scalar>>& muTable,
    const autoPtr<interpolation2DTable<scalar>>& kappaTable
)
:
    Thermo(t),
    muDict_(muDict),
    kappaDict_(kappaDict),
    muTable_(new interpolation2DTable<scalar>(muTable())),
    kappaTable_(new interpolation2DTable<scalar>(kappaTable()))
{}


template<class Thermo>
inline Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const word& name,
    const tabulatedTransport& st
)
:
    Thermo(name, st),
    muDict_(st.muDict_),
    kappaDict_(st.kappaDict_),
    muTable_(new interpolation2DTable<scalar>(st.muTable_())),
    kappaTable_(new interpolation2DTable<scalar>(st.kappaTable_()))
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo>>
Foam::tabulatedTransport<Thermo>::clone() const
{
    return autoPtr<tabulatedTransport<Thermo>>
    (
        new tabulatedTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedTransport<Thermo>>
Foam::tabulatedTransport<Thermo>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<tabulatedTransport<Thermo>>
    (
        new tabulatedTransport<Thermo>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return muTable_()(p,T);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    return kappaTable_()(p, T);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cpv(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::tabulatedTransport<Thermo>::operator=
(
    const tabulatedTransport<Thermo>& st
)
{
    Thermo::operator=(st);
    muDict_            = st.muDict_;
    kappaDict_        = st.kappaDict_;
    muTable_.reset
    (
        new interpolation2DTable<scalar>(st.muTable_())
    );

    kappaTable_.reset
    (
        new interpolation2DTable<scalar>(st.kappaTable_())
    );

}


template<class Thermo>
inline void Foam::tabulatedTransport<Thermo>::operator+=
(
    const tabulatedTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        muTable_().idaxpy(Y1,Y2,st.muTable_());

        kappaTable_().idaxpy(Y1,Y2,st.kappaTable_());
    }
}


template<class Thermo>
inline void Foam::tabulatedTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulatedTransport<Thermo> Foam::operator+
(
    const tabulatedTransport<Thermo>& st1,
    const tabulatedTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return tabulatedTransport<Thermo>
        (
            t,
            st1.muDict_,
            st1.kappaDict_,
            st1.muTable_,
            st1.kappaTable_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        interpolation2DTable<scalar> iMT;
        iMT = st1.muTable_();
        iMT.idaxpy(Y1, Y2, st2.muTable_());

        interpolation2DTable<scalar> iKT;
        iKT = st1.kappaTable_();
        iKT.idaxpy(Y1, Y2, st2.kappaTable_());

        return tabulatedTransport<Thermo>
        (
            t,
            dictionary("emptyMuDictPlus"),
            dictionary("emptyKappaDictPlus"),
            iMT,
            iKT
        );
    }
}


template<class Thermo>
inline Foam::tabulatedTransport<Thermo> Foam::operator*
(
    const scalar s,
    const tabulatedTransport<Thermo>& st
)
{
    return tabulatedTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.muDict_,
        st.kappaDict_,
        st.muTable_,
        st.kappaTable_
    );
}


// ************************************************************************* //
