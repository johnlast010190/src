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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Lee.H"
#include "disperseEulerian/phase/phase.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{
namespace massTransferModels
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable(massTransferModel, Lee, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decoupledEulerian::massTransferModels::Lee::Lee
(
    const dictionary& dict,
    const phase& phase
)
:
    massTransferModel(dict, phase)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decoupledEulerian::massTransferModels::Lee::~Lee()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::decoupledEulerian::massTransferModels::Lee::mDot() const
{
    volScalarField lambda
    (
        6/phase_.diam()*phase_.L()
        * (phase_.rhoc()*phase_.rhod()/(phase_.rhod() - phase_.rhoc()))
        * sqrt( phase_.Md()/(2 * constant::mathematical::pi * constant::thermodynamic::RR * phase_.Tsat()) )
    );

    //dimensionedScalar lambdac("lambdac", dimDensity/dimTime, phase_.dict().lookupOrDefault<scalar>("lambdac", lambda.value()));
    volScalarField lambdac
    (
        IOobject
        (
            "lambdac",
            phase_.alphad().mesh().time().timeName(),
            phase_.alphad().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase_.alphad().mesh(),
        dimensionedScalar
        (
            "Dc",
            dimDensity/dimTime,
            phase_.dict().lookupOrDefault<scalar>("lambdac", 0)
        )
    );

    if (!phase_.dict().found("lambdac"))
    {
        lambdac = lambda;
    }

    if (debug && phase_.dict().found("lambdac"))
    {
        Info<< "Lee::mDot() : computing mass transfer due to phase change based on lambdac = "
            << phase_.dict().lookup("lambdac") << endl;
    }

    const volScalarField& Tc = phase_.alphad().mesh().lookupObject<volScalarField>("T");
    const dimensionedScalar& Tsat = phase_.Tsat();

    // compute mDot (condensation minus evaporation; assumes disperse phase cools continuous phase)
    return lambdac * phase_.alphad()*(1.-phase_.alphad()) * ((Tsat - phase_.Td())/Tsat - (Tc - Tsat)/Tsat);
}


Foam::tmp<Foam::volScalarField>
Foam::decoupledEulerian::massTransferModels::Lee::Sh() const
{
    return phase_.Re() * 0;
}

// ************************************************************************* //
