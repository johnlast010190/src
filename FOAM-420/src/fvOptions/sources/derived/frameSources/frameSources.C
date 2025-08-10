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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "frameSources.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::frameSources::dimensionCheck
(
    const dimensionSet& dim1,
    const dimensionSet& dim2
) const
{
    if (dim1 != dim2)
    {
        FatalErrorInFunction
            << "UEqn dimensions do not match " << dim2
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::frameSources::frameSources
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    coorFramePtr_(nullptr),
    frameSourceFaces_
    (
        this->name(),
        this->mesh(),
        cells_,
        active_,
        dict
    ),
    frameSwitches_(coeffs_)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * //

Foam::fv::frameSources::~frameSources()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::frameSourceFaces&
Foam::fv::frameSources::getFrameSourceFaces() const
{
    return frameSourceFaces_;
}


const Foam::frameSwitches& Foam::fv::frameSources::getFrameSwitches() const
{
    return frameSwitches_;
}


Foam::tmp<Foam::vectorField> Foam::fv::frameSources::frameAcceleration
(
    const vectorField& U,
    bool addCoriolis
) const
{
    const vectorField cellCentres(UIndirectList<vector>(mesh_.C(), cells_));
    const vectorField Ucells(UIndirectList<vector>(U, cells_));

    // Note this is to provide an acceleration including non-inertial rotating
    // frame for absolute velocity U. The acceleration is not the same as
    // the acceleration of the frame itself.
    tmp<vectorField> tAccel(new vectorField(cellCentres.size(), Zero));
    vectorField& accel = tAccel.ref();
    const vector Omega = coorFrame().Omega();

    // Coriolis acceleration
    if (addCoriolis)
    {
        if (coorFrame().validParentFrame())
        {
            accel +=
                (
                    Omega
                  ^ (
                        Ucells
                      - coorFrame().parentFrame().frameVelocity(cellCentres)
                    )
                );
        }
        else
        {
            accel += (Omega ^ Ucells);
        }
    }

    // Linear acceleration is already included in dU/dt
    // The term is bellow is rotational acceleration but still isn't clear
    // is this should be included in MRF/GRF or not.
    // accel +=
    //     (
    //         coorFrame().acceleration().second()
    //      ^ (cellCentres - coorFrame().CofR())
    //     );

    return tAccel;
}


void Foam::fv::frameSources::addAcceleration
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    if (!accelerationActive())
    {
        return;
    }

    vectorField& ddtUc = ddtU.primitiveFieldRef();
    tmp<vectorField> accel(frameAcceleration(U));

    forAll(cells_, i)
    {
        ddtUc[cells_[i]] += accel()[i];
    }
}


void Foam::fv::frameSources::addAcceleration
(
    fvVectorMatrix& UEqn,
    const bool rhs
) const
{
    if (!accelerationActive())
    {
        return;
    }

    dimensionCheck(UEqn.dimensions(), dimensionSet(0, 4, -2, 0, 0));

    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    tmp<vectorField> accel(frameAcceleration(UEqn.psi()));

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Usource[celli] += V[celli]*(rhs ? accel()[i] : -accel()[i]);
    }
}


void Foam::fv::frameSources::addAcceleration
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    const bool rhs
) const
{
    if (!accelerationActive())
    {
        return;
    }

    dimensionCheck(UEqn.dimensions(), dimensionSet(1, 1, -2, 0, 0));

    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    tmp<vectorField> accel(frameAcceleration(UEqn.psi()));
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Usource[celli] +=
                V[celli]*rho[celli]*(rhs ? accel()[i] : -accel()[i]);
    }
}


void Foam::fv::frameSources::addAcceleration
(
    fvBlockMatrix<vector>& UEqn
) const
{
    if (!accelerationActive())
    {
        return;
    }

    const vector Omega = coorFrame().Omega();
    const scalarField& V = mesh_.V();
    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    forAll(cells_, i)
    {
        const label celli = cells_[i];

        blockDiag[celli](0, 1) -= Omega.component(2)*V[celli];
        blockDiag[celli](0, 2) += Omega.component(1)*V[celli];
        blockDiag[celli](1, 0) += Omega.component(2)*V[celli];
        blockDiag[celli](1, 2) -= Omega.component(0)*V[celli];
        blockDiag[celli](2, 0) -= Omega.component(1)*V[celli];
        blockDiag[celli](2, 1) += Omega.component(0)*V[celli];
    }
    vectorField& source = UEqn.source();

    if (coorFrame().validParentFrame())
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(coorFrame().parentFrameName());

        forAll(cells_, i)
        {
            const label celli = cells_[i];

            source[celli] +=
                V[celli]*(Omega ^ parentFrame.frameVelocity(mesh_.C()[celli]));
        }
    }

    // Add all accleleration except for the Coriolis term
    tmp<vectorField> accel(frameAcceleration(UEqn.psi(), false));

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        source[celli] += V[celli]*accel()[i];
    }
}


void Foam::fv::frameSources::addAcceleration
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& UEqn
) const
{
    if (!accelerationActive())
    {
        return;
    }

    const vector Omega = coorFrame().Omega();
    const scalarField& V = mesh_.V();
    const scalarField& rhoInt = rho.primitiveField();
    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        scalar rhoV = rhoInt[celli]*V[celli];

        blockDiag[celli](0, 1) -= Omega.component(2)*rhoV;
        blockDiag[celli](0, 2) += Omega.component(1)*rhoV;
        blockDiag[celli](1, 0) += Omega.component(2)*rhoV;
        blockDiag[celli](1, 2) -= Omega.component(0)*rhoV;
        blockDiag[celli](2, 0) -= Omega.component(1)*rhoV;
        blockDiag[celli](2, 1) += Omega.component(0)*rhoV;
    }

    if (coorFrame().validParentFrame())
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(coorFrame().parentFrameName());

        vectorField& source = UEqn.source();

        forAll(cells_, i)
        {
            const label celli = cells_[i];

            source[celli] +=
                V[celli]
               *rho[celli]
               *(Omega ^ parentFrame.frameVelocity(mesh_.C()[celli]));
        }
    }
    vectorField& source = UEqn.source();
    tmp<vectorField> accel(frameAcceleration(UEqn.psi(), false));
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        source[celli] += rho[celli]*V[celli]*accel()[i];
    }
}


// ************************************************************************* //
