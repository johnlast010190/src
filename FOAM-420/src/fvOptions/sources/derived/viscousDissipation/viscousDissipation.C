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
    (c) 2017 OpenCFD Ltd.
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "viscousDissipation.H"
#include "fvMatrices/fvMatrices.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(viscousDissipation, 0);

    addToRunTimeSelectionTable
    (
        option,
        viscousDissipation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::viscousDissipation::rho() const
{
    auto trho = tmp<volScalarField>::New
    (
        IOobject
        (
            "trho",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        rho_
    );

    const basicThermo* thermoPtr =
        obr_.lookupObjectPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        if (rho_.value() > 0)
        {
            return trho;
        }
        else if (rhoName_ != "none")
        {
            trho.ref() = obr_.lookupObject<volScalarField>(rhoName_);
            return trho;
        }

        FatalErrorInFunction
            << "Neither rhoName nor rho are specified."
            << exit(FatalError);
        ::abort();
    }
    else
    {
        if (!obr_.foundObject<transportModel>("transportProperties"))
        {
            FatalErrorInFunction
                << "Could not find transportProperties object to access rho "
                << exit(FatalError);
        }

        const transportModel &transport =
            obr_.lookupObject<transportModel>("transportProperties");

        tmp<volScalarField> rho(transport.rho());

        return rho;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::viscousDissipation::viscousDissipation
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "none")),
    rho_
    (
        coeffs_.lookupOrDefault
        (
            "rhoInf",
            dimensionedScalar("rho", dimDensity, 0)
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::viscousDissipation::initialise()
{
    const basicThermo* thermoPtr =
        obr_.lookupObjectPtr<basicThermo>(basicThermo::dictName);
        //mesh_.foundObject<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        fieldNames_.setSize(1, thermoPtr->he().name());
    }
    else
    {
        fieldNames_.setSize(1, "T");
    }

    if (fieldNames_.empty())
    {
        coeffs_.lookup("fields") >> fieldNames_;
    }

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    return true;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::fv::viscousDissipation::devRhoReff() const
{
    // Incompressible
    {
        const auto* turbPtr =
            obr_.lookupObjectPtr<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        if (turbPtr)
        {
            return tmp<volSymmTensorField>(rho()*turbPtr->devRhoReff());
        }
    }

    // Compressible
    {
        const auto* turbPtr =
            obr_.lookupObjectPtr<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        if (turbPtr)
        {
            return tmp<volSymmTensorField>(turbPtr->devRhoReff());
        }
    }

    FatalErrorInFunction
        << " The turbulence model is not found in the database."
        << exit(FatalError);
    ::abort();
}


void Foam::fv::viscousDissipation::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    typedef typename outerProduct<vector, vector>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const word gradUName("grad(" + UName_ + ')');

    auto tgradU = tmp<GradFieldType>::New
    (
        IOobject
        (
            "gradU",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(inv(dimTime))
    );

    // Cached?
    const GradFieldType* gradUPtr = mesh_.lookupObjectPtr<GradFieldType>(gradUName);

    if (gradUPtr)
    {
        tgradU.ref() = *gradUPtr;
    }
    else
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        tgradU.ref() = fvc::grad(U);
    }

    const auto *compTurbPtr =
        obr_.lookupObjectPtr<compressible::turbulenceModel>(
            turbulenceModel::propertiesName);

    const auto *incompTurbPtr =
        obr_.lookupObjectPtr<incompressible::turbulenceModel>(
            turbulenceModel::propertiesName);

    if (compTurbPtr)
    {
        const volScalarField rhoEpsilon("rhoEpsilon", rho * compTurbPtr->epsilon());
        const volSymmTensorField turbProdTermG(
            "turbProdTermG",
            (rho * compTurbPtr->nut()) * (dev(twoSymm(tgradU.ref()))));

        const volScalarField D(
            "D",
            ((devRhoReff() + turbProdTermG) && tgradU.ref()) - rhoEpsilon);

        eqn -= D;
    }
    else if (incompTurbPtr)
    {
        const volScalarField rhoEpsilon("rhoEpsilon", rho * incompTurbPtr->epsilon());
        const volSymmTensorField devPart("devPart", dev(twoSymm(tgradU.ref())));
        const volScalarField devLamPart(
            "devLamPart",
            -rho*incompTurbPtr->mu() * (devPart && tgradU.ref()));

        const volScalarField D(
            "D",
            devLamPart - rhoEpsilon);

        if (fieldNames_.found("T"))
        {
            if (!obr_.foundObject<transportModel>("transportProperties"))
            {
                FatalErrorInFunction
                    << "Could not find transportProperties object to access Cp "
                    << exit(FatalError);
            }

            const transportModel &transport =
                obr_.lookupObject<transportModel>("transportProperties");
            tmp<volScalarField> Cp(transport.Cp());

            eqn -= D/(rho*Cp);
        }
        else
        {
            eqn -= D;
        }
    }
}

void Foam::fv::viscousDissipation::addSup
(
    fvMatrix<scalar> &eqn,
const label fieldi
)
{
    addSup(rho(), eqn, fieldi);
}


// ************************************************************************* //
