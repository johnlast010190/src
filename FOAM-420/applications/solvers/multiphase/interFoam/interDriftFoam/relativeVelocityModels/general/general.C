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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "general/general.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(general, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, general, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::general::general
(
    const word& phaseName,
    const dictionary& dict,
    const volScalarField& alpha
)
:
    relativeVelocityModel(phaseName, dict, alpha),
    a_("a", dimless, dict),
    a1_("a1", dimless, dict),
    V0_("V0", dimVelocity, dict),
    residualAlpha_("residualAlpha", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::general::~general()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::general::correct()
{
    Udm_ =
       V0_
       *(
            exp(-a_*max(alpha_ - residualAlpha_, scalar(0)))
          - exp(-a1_*max(alpha_ - residualAlpha_, scalar(0)))
        );
}


scalar Foam::relativeVelocityModels::general::alpha23() const
{
    return ((Foam::log(a1_.value()) - Foam::log(a_.value()))
         / (a1_.value() - a_.value()) + residualAlpha_.value());
}


// ************************************************************************* //
