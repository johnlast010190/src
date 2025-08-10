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

\*---------------------------------------------------------------------------*/

#include "submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(phaseChangeModel, 0);
defineRunTimeSelectionTable(phaseChangeModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseChangeModel::phaseChangeModel
(
    surfaceFilmModel& film
)
:
    filmSubModelBase(film),
    latestMassPC_(0.0),
    totalMassPC_(0.0)
{}


phaseChangeModel::phaseChangeModel
(
    const word& modelType,
    surfaceFilmModel& film,
    const dictionary& dict
)
:
    filmSubModelBase(film, dict, typeName, modelType),
    latestMassPC_(0.0),
    totalMassPC_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

phaseChangeModel::~phaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void phaseChangeModel::correct
(
    const scalar dt,
    scalarField& availableMass,
    volScalarField& dMass,
    volScalarField& dEnergy
)
{
    if (!active())
    {
        return;
    }

    correctModel
    (
        dt,
        availableMass,
        dMass,
        dEnergy
    );

    latestMassPC_ = sum(dMass.primitiveField());
    totalMassPC_ += latestMassPC_;

    availableMass -= dMass;
    dMass.correctBoundaryConditions();

    if (writeTime())
    {
        scalar phaseChangeMass = getModelProperty<scalar>("phaseChangeMass");
        phaseChangeMass += returnReduce(totalMassPC_, sumOp<scalar>());
        setModelProperty<scalar>("phaseChangeMass", phaseChangeMass);
        totalMassPC_ = 0.0;
    }
}


void phaseChangeModel::info(Ostream& os) const
{
    const scalar massPCRate =
        returnReduce(latestMassPC_, sumOp<scalar>())
       /filmModel_.time().deltaTValue();

    scalar phaseChangeMass = getModelProperty<scalar>("phaseChangeMass");
    phaseChangeMass += returnReduce(totalMassPC_, sumOp<scalar>());

    os  << indent << "mass phase change  = " << phaseChangeMass << nl
        << indent << "vapourisation rate = " << massPCRate << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam

// ************************************************************************* //
