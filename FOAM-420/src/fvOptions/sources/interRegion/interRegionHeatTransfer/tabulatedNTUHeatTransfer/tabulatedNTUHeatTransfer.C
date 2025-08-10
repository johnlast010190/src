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

#include "tabulatedNTUHeatTransfer.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(tabulatedNTUHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            option,
            tabulatedNTUHeatTransfer,
            dictionary
        );
    }

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::tabulatedNTUHeatTransfer::geometryModeType,
        2
    >::names[] =
    {
        "calculated",
        "user"
    };
}

const Foam::NamedEnum<Foam::fv::tabulatedNTUHeatTransfer::geometryModeType, 2>
Foam::fv::tabulatedNTUHeatTransfer::geometryModelNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::interpolation2DTable<Foam::scalar>&
Foam::fv::tabulatedNTUHeatTransfer::ntuTable()
{
    if (!ntuTable_.valid())
    {
        ntuTable_.reset(new interpolation2DTable<scalar>(coeffs_));
    }

    return ntuTable_();
}


const Foam::basicThermo& Foam::fv::tabulatedNTUHeatTransfer::thermo
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


void Foam::fv::tabulatedNTUHeatTransfer::initialiseGeometry()
{
    if (Ain_ < 0 && !userGeomIsInit_)
    {
        geometryMode_ =
            geometryModelNames_.read(coeffs_.lookup("geometryMode"));

        Info<< "Region " << mesh_.name() << " " << type() << " " << name_ << " "
            << geometryModelNames_[geometryMode_] << " geometry:" << nl;

        switch (geometryMode_)
        {
            case gmCalculated:
            {
                const fvMesh& nbrMesh =
                    mesh_.time().lookupObject<fvMesh>(nbrRegionName());

                inletPatch_ = word(coeffs_.lookup("inletPatch"));
                inletPatchNbr_ = word(coeffs_.lookup("inletPatchNbr"));

                Info<< "    Inlet patch           : " << inletPatch_ << nl
                    << "    Inlet patch neighbour : " << inletPatchNbr_
                    << nl;

                label patchI = mesh_.boundary().findPatchID(inletPatch_);
                label patchINbr =
                    nbrMesh.boundary().findPatchID(inletPatchNbr_);

                scalar alpha(readScalar(coeffs_.lookup("inletBlockageRatio")));

                if ((alpha < 0) || (alpha > 1))
                {
                    FatalErrorInFunction
                        << "Inlet patch blockage ratio must be between 0 and 1"
                        << ".  Current value: " << alpha
                        << abort(FatalError);
                }

                scalar alphaNbr
                (
                    readScalar(coeffs_.lookup("inletBlockageRatioNbr"))
                );

                if ((alphaNbr < 0) || (alphaNbr > 1))
                {
                    FatalErrorInFunction
                        << "Inlet patch neighbour blockage ratio must be "
                        << "between 0 and 1.  Current value: " << alphaNbr
                        << abort(FatalError);
                }

                Info<< "    Inlet blockage ratio  : " << alpha << nl
                    << "    Inlet blockage ratio neighbour : " << alphaNbr
                    << nl;

                Ain_ =
                    (scalar(1) - alpha)
                   *gSum(mesh_.magSf().boundaryField()[patchI]);

                AinNbr_ =
                    (scalar(1) - alphaNbr)
                   *gSum(nbrMesh.magSf().boundaryField()[patchINbr]);

                scalar beta(readScalar(coeffs_.lookup("coreBlockageRatio")));

                if ((beta < 0) || (beta > 1))
                {
                    FatalErrorInFunction
                        << "Core volume blockage ratio must be between 0 and 1"
                        << ".  Current value: " << beta
                        << abort(FatalError);
                }

                Info<< "    Core volume blockage ratio : " << beta << nl;

                Vcore_ = (scalar(1) - beta)*meshInterp().V();

                break;
            }
            case gmUser:
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
                userGeomIsInit_ = true;

                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unhandled enumeration " << geometryMode_
                    << abort(FatalError);
            }
        }

        Info<< "    Inlet area local      : " << Ain_ << nl
            << "    Inlet area neighbour  : " << AinNbr_ << nl
            << "    Core volume           : " << Vcore_ << nl
            << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::fv::tabulatedNTUHeatTransfer::massFlowRate() const
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


Foam::tmp<Foam::scalarField> Foam::fv::tabulatedNTUHeatTransfer::massFlowRateNbr() const
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

Foam::fv::tabulatedNTUHeatTransfer::tabulatedNTUHeatTransfer
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
    ntuTable_(),
    geometryMode_(gmCalculated),
    Ain_(-1),
    AinNbr_(-1),
    Vcore_(-1),
    mDot_(-1),
    mDotNbr_(-1),
    Cp_(-1),
    CpNbr_(-1),
    inletPatch_("none"),
    inletPatchNbr_("none"),
    userGeomIsInit_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::tabulatedNTUHeatTransfer::~tabulatedNTUHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::tabulatedNTUHeatTransfer::calculateHtc(scalar deltaT)
{
    initialiseGeometry();

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName());

    const basicThermo& thermo = this->thermo(mesh_);
    const basicThermo& thermoNbr = this->thermo(nbrMesh);
    scalarField Cp(thermo.Cp()->internalField());
    scalarField CpNbr(interpolate(thermoNbr.Cp()->internalField()));

    // use dictionary Cp's if specified
    if (Cp_ > 0)
    {
        Cp = neg(-Cp)*Cp_;
        CpNbr = neg(-CpNbr)*CpNbr_;
    }

    // calculate mass flow rates
    const scalarField mDot(massFlowRate());
    const scalarField mDotNbr(massFlowRateNbr());

    scalarField& htcc = htc_.primitiveFieldRef();
    const interpolation2DTable<Foam::scalar>& ntuTable = this->ntuTable();

    forAll(htcc, cellI)
    {
        scalar Cpc = Cp[cellI];
        scalar CpcNbr = CpNbr[cellI];
        scalar mDotc = mDot[cellI];
        scalar mDotcNbr = mDotNbr[cellI];
        scalar Cmin = min(Cpc*mDotc, CpcNbr*mDotcNbr);
        scalar ntu = ntuTable(mDotc, mDotcNbr);

        htcc[cellI] = Cmin*ntu/Vcore_;
    }

    if (reportOn())
    {
        scalar meanMassFlow = 0.0;
        scalar meanMassFlowNbr = 0.0;
        scalar vtot = 0.0;
        scalar ntuMin = GREAT;
        scalar ntuMax = -GREAT;
        scalar ntuMean = 0.0;
        scalar cMinMin = GREAT;
        scalar cMinMax = -GREAT;
        scalar cMinMean = 0.0;

        forAll(htcc, cellI)
        {
            if (htcc[cellI] > 0)
            {
                scalar Cpc = Cp[cellI];
                scalar CpcNbr = CpNbr[cellI];
                scalar mDotc = mDot[cellI];
                scalar mDotcNbr = mDotNbr[cellI];
                scalar Cmin = min(Cpc*mDotc, CpcNbr*mDotcNbr);
                scalar ntu = ntuTable(mDotc, mDotcNbr);

                vtot += mesh_.V()[cellI];
                meanMassFlow += mDotc*mesh_.V()[cellI];
                meanMassFlowNbr += mDotcNbr*mesh_.V()[cellI];
                ntuMin = min(ntuMin, ntu);
                ntuMax = max(ntuMax, ntu);
                ntuMean += ntu*mesh_.V()[cellI];
                cMinMin = min(cMinMin, Cmin);
                cMinMax = max(cMinMax, Cmin);
                cMinMean += Cmin*mesh_.V()[cellI];
            }
        }

        reduce(meanMassFlow, sumOp<scalar>());
        reduce(meanMassFlowNbr, sumOp<scalar>());
        reduce(vtot, sumOp<scalar>());
        reduce(ntuMin, minOp<scalar>());
        reduce(ntuMax, maxOp<scalar>());
        reduce(ntuMean, sumOp<scalar>());
        reduce(cMinMin, minOp<scalar>());
        reduce(cMinMax, maxOp<scalar>());
        reduce(cMinMean, sumOp<scalar>());

        meanMassFlow /= vtot;
        meanMassFlowNbr /= vtot;
        ntuMean /= vtot;
        cMinMean /= vtot;

        Info<< "    Mean mass flow primary region   : " << meanMassFlow << nl
            << "    Mean mass flow neighbour region : " << meanMassFlowNbr << nl
            << "    Heat exchanger region volume    : " << vtot << nl
            << "    Heat exchanger overlap volume   : " << meshInterp().V() << nl
            << "    NTU based on mean flow rate     : " << ntuTable(meanMassFlow, meanMassFlowNbr) << nl
            << "    Neighbour region name           : " << nbrRegionName() << nl
            << "    NTU min value      : " << ntuMin << nl
            << "    NTU max value      : " << ntuMax << nl
            << "    NTU mean value     : " << ntuMean << nl
            << "    cMin min value     : " << cMinMin << nl
            << "    cMin max value     : " << cMinMax << nl
            << "    cMin mean value    : " << cMinMean << nl
            << endl;
    }
}


bool Foam::fv::tabulatedNTUHeatTransfer::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("U", UName_);
        coeffs_.readIfPresent("UNbr", UNbrName_);
        coeffs_.readIfPresent("rho", rhoName_);
        coeffs_.readIfPresent("rhoNbr", rhoNbrName_);

        // Force geometry re-initialisation
        Ain_ = -1;
        userGeomIsInit_ = false;
        initialiseGeometry();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
