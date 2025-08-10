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
    (c) 2013-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/adiabaticPerfectFluid/adiabaticPerfectFluid.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    Specie(obr, dict),
    p0_(readScalar(dict.subDict("equationOfState").lookup("p0"))),
    rho0_(readScalar(dict.subDict("equationOfState").lookup("rho0"))),
    gamma_(readScalar(dict.subDict("equationOfState").lookup("gamma"))),
    B_(readScalar(dict.subDict("equationOfState").lookup("B")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::adiabaticPerfectFluid<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("p0", p0_);
    dict.add("rho0", rho0_);
    dict.add("gamma", gamma_);
    dict.add("B", B_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const adiabaticPerfectFluid<Specie>& pf
)
{
    pf.write(os);
    return os;
}


// ************************************************************************* //
