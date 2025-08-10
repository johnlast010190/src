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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "AntoineExtended.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(AntoineExtended, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        AntoineExtended,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::AntoineExtended::AntoineExtended
(
    const dictionary& dict
)
:
    Antoine(dict),
    D_("D", dimless, dict),
    F_("F", dimless, dict),
    E_("E", dimless/pow(dimTemperature, F_), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::AntoineExtended::~AntoineExtended()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineExtended::pSat
(
    const volScalarField& T
) const
{
    return
        dimensionedScalar("one", dimPressure/pow(dimTemperature, D_), 1)
       *exp(A_ + B_/(C_ + T) + E_*pow(T, F_))
       *pow(T, D_);
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineExtended::pSatPrime
(
    const volScalarField& T
) const
{
    return pSat(T)*((D_ + E_*F_*pow(T, F_))/T - B_/sqr(C_ + T));
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineExtended::lnPSat
(
    const volScalarField& T
) const
{
    return
        A_
      + B_/(C_ + T)
      + D_*log(T*dimensionedScalar("one", dimless/dimTemperature, 1))
      + E_*pow(T, F_);
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineExtended::Tsat
(
    const volScalarField& p
) const
{
    NotImplemented;
}


// ************************************************************************* //
