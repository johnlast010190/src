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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "tabulatedEffectivenessHeatTransfer.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(tabulatedEffectivenessHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            option,
            tabulatedEffectivenessHeatTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::interpolation2DTable<Foam::scalar>&
Foam::fv::tabulatedEffectivenessHeatTransfer::effectivenessTable()
{
    if (!effectivenessTable_.valid())
    {
        effectivenessTable_.reset(new interpolation2DTable<scalar>(coeffs_));
    }

    return effectivenessTable_();
}


const Foam::basicThermo& Foam::fv::tabulatedEffectivenessHeatTransfer::thermo
(
    const fvMesh& mesh
) const
{
    if (!mesh.foundObject<basicThermo>(basicThermo::dictName))
    {
        FatalErrorInFunction
            << " on mesh " << mesh.name()
            << " could not find " << basicThermo::dictName
            << exit(FatalError);
    }

    return mesh.lookupObject<basicThermo>(basicThermo::dictName);
}


void Foam::fv::tabulatedEffectivenessHeatTransfer::initialiseGeometry()
{
    if (Vcore_ < 0)
    {
        if (
            coeffs_.found("mDot")
         && coeffs_.found("mDotNbr")
        )
        {
            coeffs_.lookup("mDot") >> mDot_;
            coeffs_.lookup("mDotNbr") >> mDotNbr_;
        }
        else if (
            coeffs_.found("inletPatch")
         && coeffs_.found("inletPatchNbr")
        )
        {
            inletPatch_ = word(coeffs_.lookup("inletPatch"));
            inletPatchNbr_ = word(coeffs_.lookup("inletPatchNbr"));
        }
        else
        {
            FatalErrorInFunction
                << "Please specify one of the below entry pairs " << nl
                << " mDot and mDotNbr" << nl
                << " inletPatch and inletPatchNbr" << nl
                << exit(FatalError);
        }

        if (!coeffs_.readIfPresent("Vcore", Vcore_))
        {
            Vcore_ = meshInterp().V();
        }

        if (
            coeffs_.found("Cp")
         && coeffs_.found("CpNbr")
        )
        {
            coeffs_.lookup("Cp") >> Cp_;
            coeffs_.lookup("CpNbr") >> CpNbr_;
        }

        if (
            coeffs_.found("Tin")
         && coeffs_.found("TinNbr")
        )
        {
            coeffs_.lookup("Tin") >> Tin_;
            coeffs_.lookup("TinNbr") >> TinNbr_;
        }
        else if (
            coeffs_.found("inletPatch")
         && coeffs_.found("inletPatchNbr")
        )
        {
            inletPatch_ = word(coeffs_.lookup("inletPatch"));
            inletPatchNbr_ = word(coeffs_.lookup("inletPatchNbr"));
        }
        else
        {
            FatalErrorInFunction
                << "Please specify one of the below entry pairs " << nl
                << " Tin and TinNbr" << nl
                << " inletPatch and inletPatchNbr" << nl
                << exit(FatalError);
        }
    }
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::massFlowRate() const
{
    if (mDot_ < 0)
    {
        label patchI = mesh_.boundary().findPatchID(inletPatch_);

        // Calculate scaled mass flow for primary region
        const fvsPatchVectorField& Sf =
            mesh_.Sf().boundaryField()[patchI];
        const fvPatchVectorField& Uf =
            obr_.lookupObject<volVectorField>(UName_).boundaryField()[patchI];
        const fvPatchScalarField& rhof =
            obr_.lookupObject<volScalarField>(rhoName_).boundaryField()[patchI];

        return gSum(rhof*mag(Uf & Sf));
    }

    return mDot_;
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::massFlowRateNbr() const
{
    if (mDotNbr_ < 0)
    {
        const fvMesh& nbrMesh =
            mesh_.time().lookupObject<fvMesh>(nbrRegionName());

        label patchINbr =
            nbrMesh.boundary().findPatchID(inletPatchNbr_);

        // Calculate scaled mass flow for neighbour region
        const fvsPatchVectorField& SfNbr =
            nbrMesh.Sf().boundaryField()[patchINbr];
        const fvPatchVectorField& UfNbr =
            nbrMesh.lookupObject<volVectorField>(UNbrName_).boundaryField()[patchINbr];
        const fvPatchScalarField& rhofNbr =
            nbrMesh.lookupObject<volScalarField>(rhoNbrName_).boundaryField()[patchINbr];

        return gSum(rhofNbr*mag(UfNbr & SfNbr));
    }

    return mDotNbr_;
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::inletTemperature() const
{
    if (Tin_ < 0)
    {
        label patchI = mesh_.boundary().findPatchID(inletPatch_);

        // Calculate mass flow averaged inlet temperature for primary region
        const fvsPatchVectorField& Sf =
            mesh_.Sf().boundaryField()[patchI];
        const fvPatchVectorField& Uf =
            obr_.lookupObject<volVectorField>(UName_).boundaryField()[patchI];
        const fvPatchScalarField& rhof =
            obr_.lookupObject<volScalarField>(rhoName_).boundaryField()[patchI];

        const fvPatchScalarField& Tinf =
            obr_.lookupObject<volScalarField>(TName_).boundaryField()[patchI];

        return gSum(rhof*mag(Uf & Sf)*Tinf)/gSum(rhof*mag(Uf & Sf));
    }

    return Tin_;
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::inletTemperatureNbr() const
{
    if (TinNbr_ < 0)
    {
        const fvMesh& nbrMesh =
            mesh_.time().lookupObject<fvMesh>(nbrRegionName());

        label patchINbr =
            nbrMesh.boundary().findPatchID(inletPatchNbr_);

        // Calculate mass flow averaged inlet temperature for neighbour region
        const fvsPatchVectorField& SfNbr =
            nbrMesh.Sf().boundaryField()[patchINbr];
        const fvPatchVectorField& UfNbr =
            nbrMesh.lookupObject<volVectorField>(UNbrName_).boundaryField()[patchINbr];
        const fvPatchScalarField& rhofNbr =
            nbrMesh.lookupObject<volScalarField>(rhoNbrName_).boundaryField()[patchINbr];

        const fvPatchScalarField& TinfNbr =
            nbrMesh.lookupObject<volScalarField>(TNbrName_).boundaryField()[patchINbr];

        return gSum(rhofNbr*mag(UfNbr & SfNbr)*TinfNbr)/gSum(rhofNbr*mag(UfNbr & SfNbr));
    }

    return TinNbr_;
}

Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::inletCp() const
{
    if (Cp_ < 0)
    {
        const basicThermo& thermo = this->thermo(mesh_);

        label patchI = mesh_.boundary().findPatchID(inletPatch_);

        // Calculate mass flow averaged inlet temperature for primary region
        const fvsPatchVectorField& Sf =
            mesh_.Sf().boundaryField()[patchI];
        const fvPatchVectorField& Uf =
            obr_.lookupObject<volVectorField>(UName_).boundaryField()[patchI];
        const fvPatchScalarField& rhof =
            obr_.lookupObject<volScalarField>(rhoName_).boundaryField()[patchI];

        const fvPatchScalarField Cpf =
            thermo.Cp()().boundaryField()[patchI];

        return gSum(rhof*mag(Uf & Sf)*Cpf)/gSum(rhof*mag(Uf & Sf));
    }

    return Cp_;
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::inletCpNbr() const
{
    if (CpNbr_ < 0)
    {
        const fvMesh& nbrMesh =
            mesh_.time().lookupObject<fvMesh>(nbrRegionName());

        const basicThermo& thermoNbr = this->thermo(nbrMesh);

        label patchINbr =
            nbrMesh.boundary().findPatchID(inletPatchNbr_);

        // Calculate mass flow averaged inlet temperature for neighbour region
        const fvsPatchVectorField& SfNbr =
            nbrMesh.Sf().boundaryField()[patchINbr];
        const fvPatchVectorField& UfNbr =
            nbrMesh.lookupObject<volVectorField>(UNbrName_).boundaryField()[patchINbr];
        const fvPatchScalarField& rhofNbr =
            nbrMesh.lookupObject<volScalarField>(rhoNbrName_).boundaryField()[patchINbr];

        const fvPatchScalarField CpfNbr =
            thermoNbr.Cp()().boundaryField()[patchINbr];

        return gSum(rhofNbr*mag(UfNbr & SfNbr)*CpfNbr)/gSum(rhofNbr*mag(UfNbr & SfNbr));
    }

    return CpNbr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::tabulatedEffectivenessHeatTransfer::tabulatedEffectivenessHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    interRegionHeatTransferModel(name, modelType, dict, obr),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    UNbrName_(coeffs_.lookupOrDefault<word>("UNbr", "U")),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    rhoNbrName_(coeffs_.lookupOrDefault<word>("rhoNbr", "rho")),
    TName_(coeffs_.lookupOrDefault<word>("T", "T")),
    TNbrName_(coeffs_.lookupOrDefault<word>("TNbr", "T")),
    effectivenessTable_(),
    Vcore_(-1),
    mDot_(-1),
    mDotNbr_(-1),
    Cp_(-1),
    CpNbr_(-1),
    Tin_(-1),
    TinNbr_(-1),
    inletPatch_("none"),
    inletPatchNbr_("none"),
    Qtarget_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::tabulatedEffectivenessHeatTransfer::~tabulatedEffectivenessHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::tabulatedEffectivenessHeatTransfer::calculateHtc(scalar deltaT)
{
    initialiseGeometry();
    Qtarget_ = calcTargetEnergy();

    scalarField& htcc = htc_.primitiveFieldRef();

    // rough estimate on q until deltaT is given as argument
    htcc = mag(Qtarget_/deltaT/Vcore_);
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::calcTargetEnergy()
{
    // calculate mass flow rates
    const scalar mDot(massFlowRate());
    const scalar mDotNbr(massFlowRateNbr());

    // caluclate inlet temperatures
    const scalar Tin(inletTemperature());
    const scalar TinNbr(inletTemperatureNbr());

    // caluclate inlet Cp's
    const scalar Cp(inletCp());
    const scalar CpNbr(inletCpNbr());

    // calculate Cmin
    const scalar Cmin = min(Cp*mDot, CpNbr*mDotNbr);

    // interpolate effectiveness
    const interpolation2DTable<Foam::scalar>& effectivenessTable
        = this->effectivenessTable();
    const scalar effectiveness = effectivenessTable(mDot, mDotNbr);

    Info<< "    Mean mass flow primary region   : " << mDot << nl
        << "    Mean mass flow neighbour region : " << mDotNbr << nl
        << "    Heat exchanger overlap volume   : " << meshInterp().V() << nl
        << "    Neighbour region name           : " << nbrRegionName() << nl
        << "    Effectiveness mean value : " << effectiveness << nl
        << "    cMin mean value          : " << Cmin << nl
        << "    Tin primary region       : " << Tin << nl
        << "    Tin neighbour region     : " << TinNbr << nl
        << "    Cp primary region        : " << Cp << nl
        << "    Cp neighbour region      : " << CpNbr << nl
        << endl;

    // compute q = eff * qmax
    return effectiveness * Cmin * (TinNbr - Tin);
}


Foam::scalar Foam::fv::tabulatedEffectivenessHeatTransfer::getTargetEnergy()
{
    // htc and calcTargetEnergy only ever called on master!
    if (master_)
    {
        return Qtarget_;
    }
    else
    {
        return nbrModel().getTargetEnergy();
    }
}


bool Foam::fv::tabulatedEffectivenessHeatTransfer::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("U", UName_);
        coeffs_.readIfPresent("UNbr", UNbrName_);
        coeffs_.readIfPresent("rho", rhoName_);
        coeffs_.readIfPresent("rhoNbr", rhoNbrName_);
        coeffs_.readIfPresent("T", TName_);
        coeffs_.readIfPresent("TNbr", TNbrName_);

        // Force geometry re-initialisation
        Vcore_ = -1;
        initialiseGeometry();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
