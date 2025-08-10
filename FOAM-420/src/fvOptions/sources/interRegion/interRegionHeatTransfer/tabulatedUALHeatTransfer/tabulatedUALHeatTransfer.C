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

#include "tabulatedUALHeatTransfer.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(tabulatedUALHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            option,
            tabulatedUALHeatTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::interpolation2DTable<Foam::scalar>&
Foam::fv::tabulatedUALHeatTransfer::ualTable()
{
    if (!ualTable_.valid())
    {
        ualTable_.reset(new interpolation2DTable<scalar>(coeffs_));
    }

    return ualTable_();
}


void Foam::fv::tabulatedUALHeatTransfer::initialiseGeometry()
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
            coeffs_.found("Ain")
         && coeffs_.found("AinNbr")
        )
        {
            coeffs_.lookup("Ain") >> Ain_;
            coeffs_.lookup("AinNbr") >> AinNbr_;
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
                << " Ain and AinNbr" << nl
                << " inletPatch and inletPatchNbr" << nl
                << exit(FatalError);
        }

        if (!coeffs_.readIfPresent("Vcore", Vcore_))
        {
            //Vcore_ = meshInterp().V();

            // compute Vcore from actual overlap of interpolation!
            volScalarField T
            (
                IOobject
                (
                    "Ti",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("Ti",dimTemperature,scalar(2))
            );

            tmp<volScalarField> tTmapped
            (
                new volScalarField
                (
                    IOobject
                    (
                        type() + ":Tmapped",
                        mesh_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    T
                )
            );

            volScalarField& Tmapped = tTmapped.ref();

            const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName());
            volScalarField Tnbr
            (
                IOobject
                (
                    "TiNbr",
                    nbrMesh.time().timeName(),
                    nbrMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                nbrMesh,
                dimensionedScalar("TiNbr",dimTemperature,scalar(1))
            );
            interpolate(Tnbr, Tmapped.primitiveFieldRef());

            Vcore_ = fvc::domainIntegrate(T - Tmapped).value();
            Info<< "meshInterp().V() " << meshInterp().V()
                 << " , Vcore computed " << Vcore_ << endl;
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::fv::tabulatedUALHeatTransfer::massFlowRate() const
{
    tmp<scalarField> tmDot
    (
        new scalarField(mesh_.V().size(), mDot_)
    );
    scalarField& mDot = tmDot.ref();

    if (mDot_ < 0)
    {
        if (Ain_ < 0)
        {
            label patchI = mesh_.boundary().findPatchID(inletPatch_);

            // Calculate scaled mass flow for primary region
            const fvsPatchVectorField& Sf =
                mesh_.Sf().boundaryField()[patchI];
            const fvPatchVectorField& Uf =
                obr_.lookupObject<volVectorField>(UName_).boundaryField()[patchI];
            const fvPatchScalarField& rhof =
                obr_.lookupObject<volScalarField>(rhoName_).boundaryField()[patchI];

            mDot = gSum(rhof*mag(Uf & Sf));
        }
        else
        {
            // Calculate scaled mass flow for primary region
            const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
            const volScalarField& rho = obr_.lookupObject<volScalarField>(rhoName_);

            mDot = mag(U)*rho*Ain_;
        }
    }

    return tmDot;
}


Foam::tmp<Foam::scalarField> Foam::fv::tabulatedUALHeatTransfer::massFlowRateNbr() const
{
    const fvMesh& nbrMesh =
        mesh_.time().lookupObject<fvMesh>(nbrRegionName());

    tmp<scalarField> tmDot
    (
        new scalarField(mesh_.V().size(), mDotNbr_)
    );
    scalarField& mDot = tmDot.ref();

    if (mDot_ < 0)
    {
        if (Ain_ < 0)
        {
            label patchINbr =
                nbrMesh.boundary().findPatchID(inletPatchNbr_);

            // Calculate scaled mass flow for neighbour region
            const fvsPatchVectorField& SfNbr =
                nbrMesh.Sf().boundaryField()[patchINbr];
            const fvPatchVectorField& UfNbr =
                nbrMesh.lookupObject<volVectorField>(UNbrName_).boundaryField()[patchINbr];
            const fvPatchScalarField& rhofNbr =
                nbrMesh.lookupObject<volScalarField>(rhoNbrName_).boundaryField()[patchINbr];

            mDot = gSum(rhofNbr*mag(UfNbr & SfNbr));
        }
        else
        {
            // Calculate scaled mass flow for neighbour region
            const volVectorField& UNbr =
                nbrMesh.lookupObject<volVectorField>(UNbrName_);
            const scalarField UMagNbr(mag(UNbr));
            const scalarField UMagNbrMapped(interpolate(UMagNbr));
            const scalarField& rhoNbr =
                nbrMesh.lookupObject<volScalarField>(rhoNbrName_).internalField();
            const scalarField rhoNbrMapped(interpolate(rhoNbr));

            mDot = UMagNbrMapped*rhoNbrMapped*AinNbr_;
        }
    }

    return tmDot;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::tabulatedUALHeatTransfer::tabulatedUALHeatTransfer
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
    ualTable_(),
    Ain_(-1),
    AinNbr_(-1),
    Vcore_(-1),
    mDot_(-1),
    mDotNbr_(-1),
    inletPatch_("none"),
    inletPatchNbr_("none")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::tabulatedUALHeatTransfer::~tabulatedUALHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::tabulatedUALHeatTransfer::calculateHtc(scalar deltaT)
{
    initialiseGeometry();

    // calculate mass flow rates
    const scalarField mDot(massFlowRate());
    const scalarField mDotNbr(massFlowRateNbr());

    scalarField& htcc = htc_.primitiveFieldRef();
    const interpolation2DTable<Foam::scalar>& ualTable = this->ualTable();

    forAll(htcc, cellI)
    {
        scalar mDotc = mDot[cellI];
        scalar mDotcNbr = mDotNbr[cellI];
        scalar ual = ualTable(mDotc, mDotcNbr);

        htcc[cellI] = ual/Vcore_;
    }
}


bool Foam::fv::tabulatedUALHeatTransfer::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("U", UName_);
        coeffs_.readIfPresent("UNbr", UNbrName_);
        coeffs_.readIfPresent("rho", rhoName_);
        coeffs_.readIfPresent("rhoNbr", rhoNbrName_);

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
