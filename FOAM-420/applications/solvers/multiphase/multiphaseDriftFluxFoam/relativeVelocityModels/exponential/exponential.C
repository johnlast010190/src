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
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "exponential.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, exponential, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::exponential::exponential
(
    const word& phaseName,
    const dictionary& dict,
    const volScalarField& alpha,
    const scalar& alphaMax
)
:
    relativeVelocityModel(phaseName, dict, alpha, alphaMax),
    a_("a", dimless, dict),
    V0_("V0", dimVelocity, dict),
    Vstokes_("Vstokes", dimVelocity, dict),
    residualAlpha_("residualAlpha", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::exponential::~exponential()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::exponential::correct()
{
    scalar magVstokes(mag(Vstokes_.value()));
    scalar magV0(mag(V0_.value()));

    Udm_ = Vstokes_/magVstokes *
    min
    (
        magVstokes,
        magV0*exp(-a_*max(alpha_, scalar(0)))
    );
}


// ************************************************************************* //
