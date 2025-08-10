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
    (c) 2016 OpenCFD Ltd.
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2017 CSIR
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "wallHeatFlux/wallHeatFlux.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcSnGrad.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"
#include "solidThermo/solidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFlux::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeDelimited(os, "patch");
    writeDelimited(os, "min");
    writeDelimited(os, "max");
    writeDelimited(os, "integral");
    os  << endl;
}


void Foam::functionObjects::wallHeatFlux::calcHeatFlux
(
    const volScalarField& T,
    volScalarField& wallHeatFlux
)
{
    tmp<volVectorField> tgradT = fvc::grad(T);

    volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

    volVectorField::Boundary& gradTBf = tgradT.ref().boundaryFieldRef();

    forAll(wallHeatFluxBf, patchi)
    {
        if (!wallHeatFluxBf[patchi].patch().coupled())
        {
            // Replace wall-normal component with snGrad
            gradTBf[patchi] +=
                (
                    T.boundaryField()[patchi].snGrad()
                  - (gradTBf[patchi] & mesh_.boundary()[patchi].nf())
                )*mesh_.boundary()[patchi].nf();
            wallHeatFluxBf[patchi] =
                -(boundaryKappa_[patchi].nfKappa() & gradTBf[patchi]);
        }
    }

    if (foundObject<volScalarField>("Qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("Qr");

        const volScalarField::Boundary& radHeatFluxBf = qr.boundaryField();

        forAll(wallHeatFluxBf, patchi)
        {
            if (!wallHeatFluxBf[patchi].patch().coupled())
            {
                wallHeatFluxBf[patchi] += radHeatFluxBf[patchi];
            }
        }
    }
}


void Foam::functionObjects::wallHeatFlux::calcHeatFluxIco
(
    const incompressible::turbulenceModel& turb,
    const volScalarField& T,
    volScalarField& wallHeatFlux
)
{
    tmp<volVectorField> tgradT = fvc::grad(T);

    volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

    const volVectorField::Boundary& gradTBf = tgradT().boundaryField();
    forAll(wallHeatFluxBf, patchi)
    {
        const fvPatch& patch = wallHeatFluxBf[patchi].patch();
        if (!patch.coupled())
        {
            vectorField kappaEff
            (
                patch.nf()
              * turb.alphaEff()->boundaryField()[patchi]
              * turb.rho()->boundaryField()[patchi]
              * turb.Cp()->boundaryField()[patchi]
            );
            wallHeatFluxBf[patchi] =
                (kappaEff & gradTBf[patchi]);
        }
    }

    if (foundObject<volScalarField>("Qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("Qr");

        const volScalarField::Boundary& radHeatFluxBf = qr.boundaryField();

        forAll(wallHeatFluxBf, patchi)
        {
            if (!wallHeatFluxBf[patchi].patch().coupled())
            {
                wallHeatFluxBf[patchi] += radHeatFluxBf[patchi];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFlux::wallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr(), name, typeName, dict),
//    mesh_(refCast<const fvMesh>(obr_)),
    patchSet_()
{
    volScalarField* wallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                obr(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimMass/pow3(dimTime), 0)
        )
    );

    obr().objectRegistry::store(wallHeatFluxPtr);

    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFlux::~wallHeatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }
    phaseName_ = dict.lookupOrDefault("phaseName", word::null);

    boundaryKappa_.clear();
    boundaryKappa_.setSize(mesh_.boundary().size());
    dictionary dict2(dict);
    dict2.add("kappaMethod", "fluidThermo");
    forAll(mesh_.boundary(), patchi)
    {
        boundaryKappa_.set
        (
            patchi,
            new boundaryKappa
            (
                obr(),
                mesh_.boundary()[patchi],
                dict2,
                phaseName_
            )
        );
    }

    return true;
}


bool Foam::functionObjects::wallHeatFlux::execute()
{
    volScalarField& wallHeatFlux = const_cast<volScalarField&>
    (
        obr().lookupObject<volScalarField>(type())
    );
    const word turbName =
        IOobject::groupName(turbulenceModel::propertiesName, phaseName_);
    const word thermoName =
        IOobject::groupName(basicThermo::dictName, phaseName_);


    if
    (
        boundaryKappa_[0].KMethod() == boundaryKappa::mtFluidThermo
     || boundaryKappa_[0].KMethod() == boundaryKappa::mtLookup
    )
    {
        if
        (
            obr().foundObject<compressible::turbulenceModel>(turbName)
        )
        {
            const compressible::turbulenceModel& turbModel =
                obr().lookupObject<compressible::turbulenceModel>(turbName);

            calcHeatFlux
            (
                turbModel.transport().T(),
                wallHeatFlux
            );
        }
        else if (obr().foundObject<fluidThermo>(thermoName))
        {
            const fluidThermo& thermo =
                obr().lookupObject<fluidThermo>(thermoName);

            calcHeatFlux
            (
                thermo.T(),
                wallHeatFlux
            );
        }
        else if (obr().foundObject<incompressible::turbulenceModel>(turbName))
        {
            const auto& turb =
                obr().lookupObject<incompressible::turbulenceModel>(turbName);

            calcHeatFluxIco
            (
                turb,
                lookupObject<volScalarField>("T"),
                wallHeatFlux
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unable to find turbulence or thermophysical model in the "
                << "database" << exit(FatalError);
        }
    }
    else if (boundaryKappa_[0].KMethod() == boundaryKappa::mtSolidThermo)
    {
        const solidThermo& thermo =
            obr().lookupObject<solidThermo>(thermoName);

        calcHeatFlux
        (
            thermo.T(),
            wallHeatFlux
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unimplemented method " << boundaryKappa_[0].KMethod() << nl
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallHeatFlux::write()
{
    const volScalarField& wallHeatFlux =
        obr().lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << wallHeatFlux.name() << endl;

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            file() << mesh_.time().value()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }

    return true;
}


// ************************************************************************* //
