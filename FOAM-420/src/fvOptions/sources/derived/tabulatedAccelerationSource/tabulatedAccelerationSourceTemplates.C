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

#include "tabulatedAccelerationSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::tabulatedAccelerationSource::addSup
(
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    Vector<vector> acceleration(motion_.acceleration());
    point origin(motion_.origin());

    // If gravitational force is present combine with the linear acceleration
    if (obr_.foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField& g =
            const_cast<uniformDimensionedVectorField&>
            (
                obr_.lookupObject<uniformDimensionedVectorField>("g")
            );

        const uniformDimensionedScalarField& hRef =
            obr_.lookupObject<uniformDimensionedScalarField>("hRef");

        g = g0_ - dimensionedVector("a", dimAcceleration, acceleration.x());

        dimensionedScalar ghRef
        (
            mag(g.value()) > SMALL
          ? g & (cmptMag(g.value())/mag(g.value()))*hRef
          : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
        );

        const_cast<volScalarField&>
        (
            obr_.lookupObject<volScalarField>("gh")
        ) = (g & mesh_.C()) - ghRef;

        const_cast<surfaceScalarField&>
        (
            obr_.lookupObject<surfaceScalarField>("ghf")
        ) = (g & mesh_.Cf()) - ghRef;
    }
    // ... otherwise include explicitly in the momentum equation
    else
    {
        eqn -= rho*dimensionedVector("a", dimAcceleration, acceleration.x());
    }

    dimensionedVector Omega
    (
        "Omega",
        dimensionSet(0, 0, -1, 0, 0),
        acceleration.y()
    );

    dimensionedVector dOmegaDT
    (
        "dOmegaDT",
        dimensionSet(0, 0, -2, 0, 0),
        acceleration.z()
    );

    //Coriolis + Centrifugal + Angular tabulatedAcceleration force
    eqn -=
    (
        rho*(2*Omega ^ eqn.psi())
        + rho*(Omega ^ (Omega ^ (mesh_.C()-origin)))
        + rho*(dOmegaDT ^ (mesh_.C()-origin))
    );
}


// ************************************************************************* //
