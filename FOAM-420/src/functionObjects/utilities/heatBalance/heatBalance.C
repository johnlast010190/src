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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "heatBalance/heatBalance.H"
#include "cfdTools/general/include/fvCFD.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(heatBalance, 0);
    addToRunTimeSelectionTable(functionObject, heatBalance, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::functionObjects::heatBalance::writeFileHeader(Ostream& os) const
{
    writeCommented(os, "Time");

    // Get list of patches
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchI)
    {
        if (!isA<processorFvPatch>(patches[patchI]))
        {
            writeDelimited(os,patches[patchI].name()+" [W]");
            writeDelimited(os,"Net flux [W]");

            if (compressible_)
            {
                writeDelimited(os,"DpDt [W]");
            }
            writeDelimited(os,"Volumetric rate of change [W]");
            writeDelimited(os,"Volumetric heat [J]");
        }
    }
    os << endl;
}


tmp<volScalarField> Foam::functionObjects::heatBalance::diffusivity()
{

    tmp<volScalarField> alphaEff;

    if
    (
        obr_.foundObject<incompressible::turbulenceModel>
        ("turbulenceModel")
    )
    {
        alphaEff =
        (
            obr_.lookupObject<incompressible::turbulenceModel>
            ("turbulenceModel").alphaEff()
        );
    }
    else if
    (
        obr_.foundObject<compressible::turbulenceModel>("turbulenceModel")
    )
    {
        alphaEff =
        (
            obr_.lookupObject<compressible::turbulenceModel>
            ("turbulenceModel").alphaEff()
        );
    }
    else
    {
        //no turb model available, exit with FatalError
        FatalErrorInFunction
            << "A valid turbulence model could not be found in the database."
            << Foam::abort(FatalError);
    }

    return alphaEff;

}

tmp<volScalarField> Foam::functionObjects::heatBalance::density()
{
    tmp<volScalarField> density;

    if (compressible_)
    {
        density =
        (
            obr_.lookupObject<basicThermo>(basicThermo::dictName).rho()
        );
    }
    else if
    (
        obr_.foundObject<incompressible::turbulenceModel>
        ("turbulenceModel")
    )
    {
        density =
        (
            obr_.lookupObject<incompressible::turbulenceModel>
            ("turbulenceModel").rho()
        );
    }
    else
    {
        FatalError << "FunctionObject:" << type()
                   << ": requires basicThermo or incompressible LES/RAS model."
                   << exit(FatalError);
    }

    return density;
}

tmp<volScalarField> Foam::functionObjects::heatBalance::enthalpy()
{
    const volScalarField& T = obr_.lookupObject<volScalarField>("T");
    tmp<volScalarField> enthalpy;

    if (compressible_)
    {
        enthalpy =
        (
            obr_.lookupObject<basicThermo>(basicThermo::dictName).he()
          - obr_.lookupObject<basicThermo>(basicThermo::dictName).Cp()
           *refTemp_
        );
    }
    else if
    (
        obr_.foundObject<incompressible::turbulenceModel>
        ("turbulenceModel")
    )
    {
        enthalpy =
        (
            obr_.lookupObject<incompressible::turbulenceModel>
            ("turbulenceModel").Cp()* (T-refTemp_)
        );
    }
    else
    {
        FatalError << "FunctionObject:" << type()
                   << ": requires basicThermo or incompressible LES/RAS model."
                   << exit(FatalError);
    }

    return enthalpy;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::heatBalance::heatBalance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    compressible_(false),
    phiName_("phi"),
    refTemp_("referenceTemperature", dimensionSet(0,0,0,1,0), 273.15),
    domainHeat_("domainHeat", dimensionSet(1,2,-2,0,0), 1.0),
    domainHeatOld_("domainHeatOld", dimensionSet(1,2,-2,0,0), 1.0),
    DpDt_("DpDt", dimensionSet(1,2,-3,0,0), 0.0)
{
    read(dict);
    writeFileHeader(file());

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::heatBalance::~heatBalance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::heatBalance::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    calculate();

    // Get list of patches
    const fvPatchList& patches = mesh_.boundary();

    // Write on screen and file
    Info<< "\nHeat budget: ";
    scalar totalFlux = 0;

    writeTime(file());

    forAll(heatFlux_->boundaryField(),patchI)
    {
        if (!isA<processorFvPatch>(patches[patchI]))
        {
            scalar fluxI = gSum
            (
                mesh_.magSf().boundaryField()[patchI]
                *heatFlux_->boundaryField()[patchI]
            );
            totalFlux += fluxI;
            file() << token::TAB << fluxI;
        }
    }

    Info<< "net boundary flux = " << totalFlux << " W,";

    file() << token::TAB << totalFlux;

    if (compressible_)
    {
        Info<< " contribution of DpDt = " << DpDt_.value()<< " W,";
        file() << token::TAB << DpDt_.value();
    }

    scalar ddtHeat =(domainHeat_.value() - domainHeatOld_.value())
        /obr_.time().deltaT().value();
    Info<< " volumetric rate of change = " << ddtHeat  << " W,";
    Info<< " volumetric heat = " << domainHeat_.value() << " J" << endl;

    file() << token::TAB << ddtHeat << token::TAB << domainHeat_.value();
    file() << endl;

    // Copy the present domain heat to domainHeatOld to be used
    //on next timestep
    domainHeatOld_ = domainHeat_;

    Info<< endl;

    return true;
}

void Foam::functionObjects::heatBalance::calculate()
{
    const surfaceScalarField& phi =
        obr_.lookupObject<surfaceScalarField>(phiName_);

    tmp<volScalarField> alphaEffPtr(diffusivity());
    tmp<volScalarField> hPtr(enthalpy());
    tmp<volScalarField> rhoPtr(density());

    heatFlux_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "heatFlux",
                obr_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(alphaEffPtr())*fvc::snGrad(hPtr())
        )
    );


    // Add flow flux induced part of heat flux (for inlets and outlets)
    const fvPatchList& patches = mesh_.boundary();
    surfaceScalarField::Boundary& patchHeatFlux =
        heatFlux_->boundaryFieldRef();

    forAll(patchHeatFlux, patchI)
    {
        if (!isA<processorFvPatch>(patches[patchI]))
        {
            if (compressible_)
            {
                patchHeatFlux[patchI] += -phi.boundaryField()[patchI]
                    *hPtr->boundaryField()[patchI]
                    /mesh_.magSf().boundaryField()[patchI];
            }
            else
            {
                patchHeatFlux[patchI]
                    += -rhoPtr->boundaryField()[patchI]
                    *phi.boundaryField()[patchI]
                    *hPtr->boundaryField()[patchI]
                    /mesh_.magSf().boundaryField()[patchI];
            }
        }
    }

    // Overall heat in domain
    domainHeat_ = fvc::domainIntegrate(rhoPtr()*hPtr());

    if (compressible_)
    {
        const volScalarField& p = obr_.lookupObject<volScalarField>("p");
        DpDt_ = fvc::domainIntegrate
        (
            fvc::DDt
            (
                surfaceScalarField("phiU", phi/fvc::interpolate(rhoPtr())),
                p
            )
        );
    }
}

dimensionedScalar Foam::functionObjects::heatBalance::calcInitialDomainHeat()
{
    tmp<volScalarField> hPtr(enthalpy());
    tmp<volScalarField> rhoPtr(density());

    // Overall heat in domain
    return(fvc::domainIntegrate(rhoPtr()*hPtr()));
}

bool Foam::functionObjects::heatBalance::write()
{
    return true;
}


bool Foam::functionObjects::heatBalance::read(const dictionary& dict)
{
    Log << type() << " " << name() <<  " read:" << nl;

    if (obr_.foundObject<basicThermo>(basicThermo::dictName))
    {
        compressible_ = true;
    }
    else if
    (
        obr_.foundObject<incompressible::turbulenceModel>
        ("turbulenceModel")
    )
    {
        compressible_= false;
    }

    // Check reference temp level
    if (dict.found("referenceTemperature"))
    {
        refTemp_ = dict.lookup("referenceTemperature");
        Info<< "    Using " << refTemp_
             <<" K as zero level for heat at inlets and outlets." << endl;
    }
    else
    {
        Info<< "    Using 0 deg Celsius as zero level for heat at inlets"
             << " and outlets." << nl
             << "    Set user defined reference temperature with keyword"
             << " referenceTemperature."
             << endl;
    }

    // Setup phi
    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }
    else
    {
        phiName_="phi";
    }

    // Calculate initial value for domain heat
    domainHeatOld_ = calcInitialDomainHeat();

    return true;
}


// ************************************************************************* //
