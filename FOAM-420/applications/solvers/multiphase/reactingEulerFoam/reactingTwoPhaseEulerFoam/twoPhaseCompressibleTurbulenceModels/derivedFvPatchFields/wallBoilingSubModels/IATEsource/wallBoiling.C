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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "wallBoiling.H"
#include "phaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.H"
#include "derivedFvPatchFields/alphatWallBoilingWallFunction/alphatWallBoilingWallFunctionFvPatchScalarField.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(wallBoiling, 0);
    addToRunTimeSelectionTable(IATEsource, wallBoiling, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::wallBoiling::wallBoiling
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::IATEsources::wallBoiling::R
(
    const volScalarField& alphai,
    volScalarField& kappai
) const
{
    volScalarField::Internal R
    (
        IOobject
        (
            "wallBoiling:R",
            phase().time().timeName(),
            phase().mesh()
        ),
        phase().mesh(),
        dimensionedScalar("R", dimless/dimTime, 0)
    );

    volScalarField::Internal Rdk
    (
        IOobject
        (
            "wallBoiling:Rdk",
            phase().time().timeName(),
            phase().mesh()
        ),
        phase().mesh(),
        dimensionedScalar("Rdk", kappai.dimensions()/dimTime, 0)
    );

    const phaseCompressibleTurbulenceModel& turbulence =
        const_cast<phaseCompressibleTurbulenceModel&>
        (
            phase().db().lookupObject<phaseCompressibleTurbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    otherPhase().name()
                )
            )
        );

    const tmp<volScalarField> talphat(turbulence.alphat());
    const volScalarField::Boundary& alphatBf = talphat().boundaryField();

    tmp<scalarField> trho(phase().rho());
    const scalarField& rho = trho();

    typedef compressible::alphatWallBoilingWallFunctionFvPatchScalarField
        alphatWallBoilingWallFunction;

    forAll(alphatBf, patchi)
    {
        if
        (
            isA<alphatWallBoilingWallFunction>(alphatBf[patchi])
        )
        {
            const alphatWallBoilingWallFunction& alphatw =
                refCast<const alphatWallBoilingWallFunction>(alphatBf[patchi]);

            const scalarField& dmdt = alphatw.dmdt();
            const scalarField& dDep = alphatw.dDeparture();

            const labelList& faceCells = alphatw.patch().faceCells();

            forAll(alphatw, facei)
            {
                if (dmdt[facei] > SMALL)
                {
                    const label faceCelli = faceCells[facei];
                    R[faceCelli] =
                        dmdt[facei]/(alphai[faceCelli]*rho[faceCelli]);
                    Rdk[faceCelli] = R[faceCelli]*(6.0/dDep[facei]);
                }
            }
        }
    }

    return Rdk - fvm::Sp(R, kappai);
}


// ************************************************************************* //
