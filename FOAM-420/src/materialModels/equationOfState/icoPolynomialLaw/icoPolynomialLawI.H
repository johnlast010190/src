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
    (c) 2016-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "icoPolynomialLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::icoPolynomialLaw>
Foam::icoPolynomialLaw::clone() const
{
    return autoPtr<icoPolynomialLaw>(new icoPolynomialLaw(*this));
}


inline Foam::autoPtr<Foam::icoPolynomialLaw> Foam::icoPolynomialLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<icoPolynomialLaw>
    (
        new icoPolynomialLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::icoPolynomialLaw::rho(const scalar T) const
{
    return max(rhoMin_, min(rhoMax_, rhoCoeffs_.value(T)));
}


inline Foam::scalar Foam::icoPolynomialLaw::CpMCv
(
    const scalar p,
    const scalar T
) const
{
    scalar drdt = rhoCoeffs_.derivative(T);
    scalar rw = rhoCoeffs_.value(T);

    if (rw > rhoMax_)
    {
        rw = rhoMax_;
        drdt = 0;
    }
    else if (rw < rhoMin_)
    {
        rw = rhoMin_;
        drdt = 0;
    }

    return -(p/sqr(rw))*drdt;
}


// ************************************************************************* //
