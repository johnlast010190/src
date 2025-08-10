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
    (c) 2010-2022 Esi Ltd.
    (c) 2013-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermalSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "transportModel/transportModel.H"
#include "basicThermo/basicThermo.H"
#include "solidThermo/solidThermo.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(thermalSource, 0);
        addToRunTimeSelectionTable
        (
            option,
            thermalSource,
            dictionary
        );
    }

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::thermalSource::volumeModeType,
        2
    >::names[] =
    {
        "absolute",
        "specific"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::thermalSource::heatTransferModeType,
        4
    >::names[] =
    {
        "uniformFlux",
        "fixedHTC",
        "scaledHTC",
        "scaledTinf"
    };
}

const Foam::NamedEnum<Foam::fv::thermalSource::volumeModeType, 2>
Foam::fv::thermalSource::volumeModeTypeNames_;

const Foam::NamedEnum<Foam::fv::thermalSource::heatTransferModeType, 4>
Foam::fv::thermalSource::heatTransferModeTypeNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::fv::thermalSource::htcScaled(const volScalarField& he) const
{
    if (he.dimensions() == dimEnergy/dimMass)
    {
        if (obr_.foundObject<basicThermo>(basicThermo::dictName))
        {
            const basicThermo& thermo =
               obr_.lookupObject<basicThermo>(basicThermo::dictName);

            return tmp<scalarField>
            (
                new scalarField(htc_/thermo.Cpv(), cells_)
            );
        }
        else if (obr_.foundObject<basicThermo>(solidThermo::dictName))
        {
            const solidThermo& thermo =
               obr_.lookupObject<solidThermo>(solidThermo::dictName);

            return tmp<scalarField>
            (
                new scalarField(htc_/thermo.Cpv(), cells_)
            );
        }
        else
        {
            FatalErrorInFunction
                << " On mesh " << mesh_.name()
                << " could not find object basicThermo or solidThermo."
                << " The available objects are: "
                << mesh_.names()
                << exit(FatalError);
        }
    }
    else if (he.dimensions() == dimTemperature)
    {
        return tmp<scalarField>
        (
            new scalarField(cells_.size(), htc_)
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported dimensions for field " << he.name() << "." << nl
            << "Must be either dimTemperature or dimEnergy/dimMass."
            << exit(FatalError);
    }

    return tmp<scalarField>();
}


void Foam::fv::thermalSource::updateRate(const fvMatrix<scalar>& eqn)
{
    //update variable htc or Tinf for fixed Q
    const volScalarField& T(obr_.lookupObject<volScalarField>(TName_));
    scalar sumDT = 0; // used for reporting

    //get new Function1 values
    const scalar t = obr_.time().value();
    if (heatTransferMode_ == hmFixedHTC)
    {
        htc_ = htcF1_->value(t);
        Tinf_ = TinfF1_->value(t);

        // input sanity check
        if (htc_ < 0)
        {
            FatalErrorInFunction
                << "HTC must be >= 0"
                << ". Current value: " << htc_ << endl;
        }
    }
    else if (heatTransferMode_ == hmScaledHTC)
    {
        Q_ = QF1_->value(t);
        Tinf_ = TinfF1_->value(t);
    }
    else if (heatTransferMode_ == hmScaledTinf)
    {
        Q_ = QF1_->value(t);
        htc_ = htcF1_->value(t);

        // input sanity check
        if (htc_ < SMALL)
        {
            FatalErrorInFunction
                << "HTC must be positive and larger than " << SMALL
                << ". Current value: " << htc_ << endl;
        }
    }
    else if (heatTransferMode_ == hmUniformFlux)
    {
        Q_ = QF1_->value(t);
    }


    if (heatTransferMode_ == hmScaledHTC)
    {
        //scale htc such that Q is satisfied given Tinf
        //Note: htc will be in the same volumeMode as the given Q
        forAll(cells_, celli)
        {
            label i = cells_[celli];
            sumDT += (Tinf_ - T[i])*mesh_.V()[i];
        }

        reduce(sumDT, sumOp<scalar>());
        sumDT /= V_;

        if ((mag(sumDT) > DTmin_) && (sign(sumDT) == sign(Q_)))
        {
            htc_ = Q_ / sumDT;
        }
        else
        {
            htc_ = mag(Q_) / DTmin_;

            WarningInFunction << "Aggregate temperature differential "
                << "in zone " << cellSetName() << " is insufficient to provide"
                << " specified heat transfer." << nl
                << "    " << (Q_ >= 0 ? "Increase" : "Decrease")
                << " mean bulk temperature Tinf or switch to alternative"
                << " heat transfer mode."
                << endl;
        }
    }

    if (heatTransferMode_ == hmScaledTinf)
    {
        //scale Tinf such that Q is satisfied given htc
        //Note: htc must be specified in the same volumeMode as Q
        scalar Tmean = 0;

        forAll(cells_, celli)
        {
            label i = cells_[celli];
            Tmean += T[i]*mesh_.V()[i];
        }

        reduce(Tmean, sumOp<scalar>());
        Tmean /= V_;

        Tinf_ = Q_/htc_ + Tmean;

        // Tinf is allowed to be Tmean but only if Q is (close to) zero
        // to support turning off the heat source.
        if (mag(Q_/htc_) < DTmin_ && mag(Q_) < SMALL)
        {
            Tinf_ = Tmean + sign(Q_)*DTmin_;
            WarningInFunction << "Aggregate temperature differential "
                << "in zone " << cellSetName() << " is insufficient to provide"
                << " specified heat transfer." << nl
                << "     Decrease heat transfer coefficient htc or switch to"
                << " alternative heat transfer mode."
                << endl;
        }

        if (Tinf_ < 0)
        {
            Tinf_ = 0;
            WarningInFunction << "To produce the specified heat transfer "
                << "in zone " << cellSetName() << " requires negative Tinf." << nl
                << "     Increase heat transfer coefficient htc or switch to"
                << " alternative heat transfer mode."
                << endl;
        }
        sumDT = Q_/htc_;
    }

    if (reportOn_ && heatTransferMode_ != hmUniformFlux)
    {
        Info<< "    " << "Heat transfer coefficient: " << htc_
             << (volumeMode_ == vmAbsolute ? " [W/K]" : " [W/K/m3]") << nl
             << "    " << "Mean temperature differential: "
             << sumDT << " [K]" << nl
             << "    " << "Current heat transfer rate: " << htc_ * sumDT
             << (volumeMode_ == vmAbsolute ? " [W]" : " [W/m3]") << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::thermalSource::thermalSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
    heatTransferMode_(hmUniformFlux),
    QF1_(nullptr),
    Q_(0),
    htcF1_(nullptr),
    htc_(0),
    TinfF1_(nullptr),
    Tinf_(0),
    TName_("T"),
    reportOn_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::thermalSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() != dimEnergy/dimTime)
    {
        FatalErrorInFunction
            << "Unsupported dimensions for " << eqn.psi().name() << " equation." << nl
            << "The dimensions of the energy equation must be dimEnergy/dimTime."
            << exit(FatalError);
    }
    updateRate(eqn);

    const volScalarField& he = eqn.psi();

    volScalarField::Internal Su
    (
        IOobject
        (
            "Su",
            mesh_.time().timeName(),
            obr_
        ),
        mesh_,
        dimensionedScalar("Su", eqn.dimensions()/dimVolume, 0)
    );
    volScalarField::Internal Sp
    (
        IOobject
        (
            "Sp",
            mesh_.time().timeName(),
            obr_
        ),
        mesh_,
        dimensionedScalar("Sp", Su.dimensions()/he.dimensions(), 0)
    );

    //take into account volumeMode_ & heatTransferMode_
    if (heatTransferMode_ == hmUniformFlux)
    {
        UIndirectList<scalar>(Su, cells_) = Q_/VDash_;
    }
    else
    {
        tmp<scalarField> htc(htcScaled(he));
        const scalarField hei(he.primitiveField(), cells_);
        const scalarField& Ti =
            obr_.lookupObject<volScalarField>(TName_).primitiveField();
        const scalarField TiCells(Ti, cells_);

        UIndirectList<scalar>(Su, cells_) = (htc_*(Tinf_-TiCells) + htc()*hei)/VDash_;
        UIndirectList<scalar>(Sp, cells_) = -htc/VDash_;
    }

    if (reportOn_)
    {
        const volScalarField& T =
            obr_.lookupObject<volScalarField>(TName_);

        if (heatTransferMode_ != hmUniformFlux)
        {
            volScalarField::Internal htc
            (
                IOobject
                (
                    "htc",
                    mesh_.time().timeName(),
                    obr_
                ),
                mesh_,
                dimensionedScalar("htc", dimless, 0)
            );
            UIndirectList<scalar>(htc, cells_) = htc_;

            const dimensionedScalar Tinf("Tinf", T.dimensions(), Tinf_);
            const dimensionedScalar energy =
                fvc::domainIntegrate(htc_*(Tinf - T)/VDash_);

            Info<< "Total energy source [J/s] : " <<  energy.value() << endl;
        }
        else
        {
            Info<< "Total energy source [J/s] : " <<  Q_/VDash_ << endl;
        }
    }

    eqn += Su + fvm::SuSp(Sp, he);
}


void Foam::fv::thermalSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn, fieldI);
}


bool Foam::fv::thermalSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        volumeMode_ = volumeModeTypeNames_.read(coeffs_.lookup("volumeMode"));
        heatTransferMode_
            = heatTransferModeTypeNames_.read(coeffs_.lookup("heatTransferMode"));

        fieldNames_.setSize(1);
        fieldNames_[0] = word(coeffs_.lookup("fieldName"));

        TName_ = coeffs_.lookupOrDefault<word>("TName", "T");
        reportOn_ = coeffs_.lookupOrDefault<Switch>("reportOn", false);

        applied_.setSize(1, false);

        // Set volume normalisation
        if (volumeMode_ == vmAbsolute)
        {
            VDash_ = V_;
        }

        if (heatTransferMode_ == hmFixedHTC)
        {
            htcF1_ = Function1<scalar>::New("htc", coeffs_);
            TinfF1_ = Function1<scalar>::New("Tinf", coeffs_);
        }
        else if (heatTransferMode_ == hmScaledHTC)
        {
            QF1_ = Function1<scalar>::New("Q", coeffs_);
            TinfF1_ = Function1<scalar>::New("Tinf", coeffs_);
            DTmin_ = coeffs_.lookupOrDefault<scalar>("dTmin", 1);
        }
        else if (heatTransferMode_ == hmScaledTinf)
        {
            QF1_ = Function1<scalar>::New("Q", coeffs_);
            htcF1_ = Function1<scalar>::New("htc", coeffs_);
            DTmin_ = coeffs_.lookupOrDefault<scalar>("dTmin", 1);
        }
        else if (heatTransferMode_ == hmUniformFlux)
        {
            QF1_ = Function1<scalar>::New("Q", coeffs_);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown heatTransferMode type "
                << heatTransferModeTypeNames_[heatTransferMode_]
                << ". Valid heatTransferMode types are:" << nl
                << heatTransferModeTypeNames_
                << exit(FatalError);
        }

        if
        (
            heatTransferMode_ == hmScaledHTC
         || heatTransferMode_ == hmScaledTinf
        )
        {
            if (DTmin_ < SMALL)
            {
                FatalErrorInFunction
                    << "DTmin must be positive and larger than " << SMALL
                    << ". Current value: " << DTmin_ << endl;
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
