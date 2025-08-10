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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cosine.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{
    defineTypeNameAndDebug(cosine, 0);
    addToRunTimeSelectionTable
    (
        partitioningModel,
        cosine,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
cosine::cosine(const dictionary& dict)
:
    partitioningModel(),
    alphaLiquid1_(readScalar(dict.lookup("alphaLiquid1"))),
    alphaLiquid0_(readScalar(dict.lookup("alphaLiquid0")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
cosine::~cosine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::partitioningModels::
cosine::fLiquid
(
    const scalarField& alphaLiquid
) const
{
    return
        pos0(alphaLiquid1_ - alphaLiquid)
       *(
            neg(alphaLiquid0_ - alphaLiquid)
           *(
                0.5
               *(
                    1 - cos
                    (
                        constant::mathematical::pi
                       *(alphaLiquid - alphaLiquid0_)
                       /(alphaLiquid1_ - alphaLiquid0_)
                    )
                )
            )
        )
      + neg(alphaLiquid1_ - alphaLiquid);
}


void Foam::wallBoilingModels::partitioningModels::
cosine::write(Ostream& os) const
{
    partitioningModel::write(os);
    os.writeEntry("alphaLiquid1", alphaLiquid1_);
    os.writeEntry("alphaLiquid0", alphaLiquid0_);
}


// ************************************************************************* //
