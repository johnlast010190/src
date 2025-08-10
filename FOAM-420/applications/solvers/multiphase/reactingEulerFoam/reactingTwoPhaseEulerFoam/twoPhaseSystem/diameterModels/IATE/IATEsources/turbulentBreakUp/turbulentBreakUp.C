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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "turbulentBreakUp.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(turbulentBreakUp, 0);
    addToRunTimeSelectionTable(IATEsource, turbulentBreakUp, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::turbulentBreakUp::
turbulentBreakUp
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    Cti_("Cti", dimless, dict),
    WeCr_("WeCr", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::IATEsources::turbulentBreakUp::R
(
    const volScalarField& alphai,
    volScalarField& kappai
) const
{
    volScalarField::Internal R
    (
        IOobject
        (
            "turbulentBreakUp:R",
            iate_.phase().time().timeName(),
            iate_.phase().mesh()
        ),
        iate_.phase().mesh(),
        dimensionedScalar("R", kappai.dimensions()/dimTime, 0)
    );

    const scalar Cti = Cti_.value();
    const scalar WeCr = WeCr_.value();
    const volScalarField Ut(this->Ut());
    const volScalarField We(this->We());

    forAll(R, celli)
    {
        if (We[celli] > WeCr)
        {
            R[celli] =
               2*Cti*Ut[celli]*sqrt(1 - WeCr/We[celli])*exp(-WeCr/We[celli]);
        }
    }

    return fvm::Su(R, kappai);
}


// ************************************************************************* //
