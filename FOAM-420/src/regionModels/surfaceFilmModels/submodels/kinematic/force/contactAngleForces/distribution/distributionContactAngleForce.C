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
    (c) 2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/kinematic/force/contactAngleForces/distribution/distributionContactAngleForce.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(distributionContactAngleForce, 0);
addToRunTimeSelectionTable(force, distributionContactAngleForce, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distributionContactAngleForce::distributionContactAngleForce
(
    surfaceFilmModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    rndGen_(),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

distributionContactAngleForce::~distributionContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> distributionContactAngleForce::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":theta",
                filmModel_.time().timeName(),
                filmModel_.regionMesh()
            ),
            filmModel_.regionMesh(),
            dimensionedScalar("0", dimless, 0)
        )
    );

    volScalarField& theta = ttheta.ref();
    volScalarField::Internal& thetai = theta.ref();

    forAll(thetai, celli)
    {
        thetai[celli] = distribution_->sample();
    }

    forAll(theta.boundaryField(), patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            fvPatchField<scalar>& thetaf = theta.boundaryFieldRef()[patchi];

            forAll(thetaf, facei)
            {
                thetaf[facei] = distribution_->sample();
            }
        }
    }

    return ttheta;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
