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
    (c) 2013-2014 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "incompressibleTurbulenceModel.H"
#include "turbulentTransportModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatKaderWallFunction/alphatKaderWallFunctionFvPatchScalarField.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTurbulenceModel, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::alphatReadIfPresent
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    IOobject alphatHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (alphatHeader.typeHeaderOk<volScalarField>(true))
    {
        tmp<volScalarField> talphat(new volScalarField(alphatHeader, mesh));

        return talphat;
    }
    else
    {
        return tmp<volScalarField>(nullptr);
    }
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::alphatAutoCreate
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    Info<< "--> Creating " << fieldName
        << " to employ run-time selectable wall functions"
        << endl;

    const fvBoundaryMesh& bm = mesh.boundary();

    wordList alphatBoundaryTypes(bm.size());

    forAll(bm, patchI)
    {
        if (isA<wallFvPatch>(bm[patchI]))
        {
            alphatBoundaryTypes[patchI] =
                incompressible::
                    alphatKaderWallFunctionFvPatchScalarField::typeName;
        }
        else
        {
            alphatBoundaryTypes[patchI] =
                calculatedFvPatchField<scalar>::typeName;
        }
    }

    tmp<volScalarField> alphat
    (
        new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimArea/dimTime, 0.0),
            alphatBoundaryTypes
        )
    );

    //Only write for restarts
    //Info<< "    Writing new " << fieldName << endl;
    //alphat().write();

    return alphat;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTurbulenceModel::incompressibleTurbulenceModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    turbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    ),
    alphatPtr_(nullptr)
{
    tmp<volScalarField> alphaTemp
    (
        alphatReadIfPresent("alphat", mesh_)
    );

    if (alphaTemp.valid())
    {
        alphatPtr_ = alphaTemp.ptr();
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::incompressibleTurbulenceModel::~incompressibleTurbulenceModel()
{
    deleteDemandDrivenData(alphatPtr_);
}


const Foam::volScalarField& Foam::incompressibleTurbulenceModel::alphat() const
{
    if (!alphatPtr_)
    {
        //create alphat when needed
        alphatPtr_ = alphatAutoCreate("alphat", mesh_).ptr();
        *alphatPtr_ = nut()/Prt();
        alphatPtr_->correctBoundaryConditions();
    }

    //backward compatibility
    if (alphatPtr_->dimensions() == dimMass / dimLength / dimTime)
    {
        *alphatPtr_ /= rho();
    }


    return *alphatPtr_;
}

const Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::alphaEff() const
{
    return (alphat() + alphaLam());
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::mu() const
{
    return nu();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::mu(const label patchi) const
{
    return nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::mut() const
{
    return nut();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::mut(const label patchi) const
{
    return nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::muEff() const
{
    return nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::muEff(const label patchi) const
{
    return nuEff(patchi);
}


void Foam::incompressibleTurbulenceModel::correct()
{
    turbulenceModel::correct();

    // this will be superfluous for laminar model, alphat=0
    if (alphatPtr_)
    {
        *alphatPtr_ = nut()/Prt();
        alphatPtr_->correctBoundaryConditions();
    }
}

Foam::tmp<Foam::scalarField> Foam::incompressibleTurbulenceModel::uWallCoeffs
(
    const label& patchI
)
{
    const label size = mesh_.boundaryMesh()[patchI].size();

    tmp<scalarField> tuWallCoeffs(new scalarField(size));
    scalarField& uWallCoeffs = tuWallCoeffs.ref();

    const scalarField& ypatch = this->y()[patchI];

    uWallCoeffs = muEff(patchI)/ypatch;

    return tuWallCoeffs;
}

// ************************************************************************* //
