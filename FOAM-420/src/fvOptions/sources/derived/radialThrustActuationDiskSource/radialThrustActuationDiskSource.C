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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2013 Kevin Maki
    (c) 2023 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "radialThrustActuationDiskSource.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
template<>
const char* NamedEnum
<
    fv::radialThrustActuationDiskSource::diskModeType,
    2
>::names[] =
{
    "fixedThrust",
    "owCurves"
};

template<>
const char* NamedEnum
<
    fv::radialThrustActuationDiskSource::propulsionModeType,
    2
>::names[] =
{
    "speed",
    "force"
};

namespace fv
{
    defineTypeNameAndDebug(radialThrustActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        radialThrustActuationDiskSource,
        dictionary
    );
}
}

const Foam::NamedEnum<Foam::fv::radialThrustActuationDiskSource::diskModeType, 2>
    Foam::fv::radialThrustActuationDiskSource::diskModeTypeNames_;

const Foam::NamedEnum<Foam::fv::radialThrustActuationDiskSource::propulsionModeType, 2>
    Foam::fv::radialThrustActuationDiskSource::propulsionModeTypeNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::radialThrustActuationDiskSource::radialThrustActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    constantThrustActuationDiskSource(name, modelType, dict, obr),
    sampledSurfaces_
    (
        (!dict.found("sampleDisk"))
        ?
        sampledSurfaces(name, mesh_, dict)
        :
        sampledSurfaces(name, mesh_, dict.subDict("sampleDisk"))
    ),
    controller_(mesh_, dict),
    PoverD_(readScalar(coeffs_.lookup("PoverD"))),
    R_(readScalar(coeffs_.lookup("R"))),
    rh_(readScalar(coeffs_.lookup("rh"))),
    propPosition_(coeffs_.lookup("propPosition")),
    propOrientation_(coeffs_.lookup("propOrientation")),
    radialCoeffs_(coeffs_.lookup("coeffs")),
    shipResistanceThrust_(coeffs_.lookupOrDefault("shipResistanceThrust", false)),
    tr_(0.0),
    Kts_(0.0),
    n_(0.0),
    updateWake_(false),
    sampleWake_(coeffs_.lookupOrDefault<bool>("sampleWake", false)),
    uDisk_(Zero)
{
    Info<< "Creating radial actuation disk zone: " << name_ << endl;

    IFstream propertiesFile
    (
        mesh_.time().timePath()/"uniform"/(this->name() + "Properties")
    );

    if (propertiesFile.good())
    {
        Info<< "    Reading propeller position from file" << endl;
        dictionary propertiesDict(dictionary::null, propertiesFile);
        propertiesDict.lookup("propPosition") >> propPosition_;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::radialThrustActuationDiskSource::checkData() const
{
    if (R_ <= VSMALL)
    {
        FatalErrorInFunction
           << "disk radius is approximately zero"
           << exit(FatalIOError);
    }
}

bool Foam::fv::radialThrustActuationDiskSource::updatePropellerPosition(const vector& propPosition)
{
    propPosition_ = propPosition;
    return true;
}

bool Foam::fv::radialThrustActuationDiskSource::updatePropellerOrientation(const quaternion& propOrientation)
{
    propOrientation_ = propOrientation;
    return true;
}


bool Foam::fv::radialThrustActuationDiskSource::calcSelfPropulsionPoint() const
{
    if (updateWake_)
    {
        return true;
    }

    vector dir = coorFramePtr_->coorSys().e1();
    scalar va = (1-w_)*(uShip_ & dir);
    scalar c = T_/(rho_*pow(2*R_,2)*pow(va,2));

    List<scalar> cJSquare = c*pow(advanceCoeffs_,2);
    List<scalar> diffCurves = thrustCoeffs_ - cJSquare;

    scalar Jts = 0;
    label idx = 1;

    while (idx < thrustCoeffs_.size())
    {
        if (neg(diffCurves[idx]))
        {
            vector2D p1 = vector2D(thrustCoeffs_[idx-1], advanceCoeffs_[idx-1]);
            vector2D p2 = vector2D(thrustCoeffs_[idx], advanceCoeffs_[idx]);
            vector2D p3 = vector2D(cJSquare[idx-1], advanceCoeffs_[idx-1]);
            vector2D p4 = vector2D(cJSquare[idx], advanceCoeffs_[idx]);

            scalar denom = (p1.x() - p2.x())*(p3.y() - p4.y()) - (p1.y() - p2.y())*(p3.x() - p4.x());
            scalar numx = (p1.x()*p2.y() - p1.y()*p2.x())*(p3.x() - p4.x()) - (p1.x() - p2.x())*(p3.x()*p4.y() - p3.y()*p4.x());
            scalar numy = (p1.x()*p2.y() - p1.y()*p2.x())*(p3.y() - p4.y()) - (p1.y() - p2.y())*(p3.x()*p4.y() - p3.y()*p4.x());

            if (denom==0)
            {
                FatalErrorInFunction
                << "Cannot find intersection between KT and cJ^2 curves "
                << "Open-Water Curves are parallel or coincident: check input data! "
                << exit(FatalIOError);
            }
            else
            {
                Kts_ = numx/denom;
                Jts = numy/denom;
                n_ = mag(va/(Jts*2*R_));
                const_cast<radialThrustActuationDiskSource*>(this)
                    ->updateThrust(Kts_*rho_*pow(n_,2)*pow(2*R_,4));

                Info<< "KT from OW curves " << Kts_ << endl;
                Info<< "Thrust from OW curves " << (Kts_*rho_*pow(n_,2)*pow(2*R_,4)) << endl;
                Info<< "J " << Jts << endl;
                Info<< "Rotational speed " << n_ << endl;

                updateWake_ = true;
                return true;
            }
        }
        idx++;
    }

    return false;
}


void Foam::fv::radialThrustActuationDiskSource::sampleDisk() const
{
    coorFramePtr_->updateState();

    switch (diskType_)
    {
        case mdFixed:
        {
            break;
        }
        case mdOWCurves:
        {
            if (!coeffs_.found("wakeFraction") || coeffs_.found("sampleDisk"))
            {
                if (!mesh_.time().writeTime())
                {
                    sampledSurfaces_.update();
                }
                else
                {
                    // Surface updated is done when writing
                    sampledSurfaces_.write();
                }

                // Return sample field on the tri surface
                List<Field<vector>> sampleU;
                sampledSurfaces_.returnField<volVectorField>(sampleU);

                // Get velocity field - hard-coded
                uDisk_ = sampledSurfaces_.returnAverage<volVectorField>(sampleU)[0];
            }
            break;
        }
    }
}


void Foam::fv::radialThrustActuationDiskSource::wakeFraction() const
{
    scalar t = mesh_.time().value();

    if (t >= tr_ && !updateWake_ && !coeffs_.found("wakeFraction"))
    {
        coorFramePtr_->updateState();
        vector dir = coorFramePtr_->coorSys().e1();
        w_ = (uDisk_ & dir)/(uShip_ & dir);
        Info<< "wake fraction " << w_ << endl;

        calcSelfPropulsionPoint();
    }
}


void Foam::fv::radialThrustActuationDiskSource::controlDisk() const
{
    scalar t = mesh_.time().value();

    if (coeffs_.isDict("controller") && t >= tr_)
    {
        scalar setPoint = 0.0;
        scalar control = 0.0;
        scalar target = 0.0;

        switch (propulsionType_)
        {
            case mdSpeed:
            {
                if (coorFramePtr_)
                {
                    coorFramePtr_->updateState();
                    vector dir = coorFramePtr_->coorSys().e1();
                    setPoint = coorFramePtr_->velocity().first() & dir;
                    target = (uShip_ & dir);
                    control = controller_.executePID(setPoint, target);
                    n_ += control;
                }
                else
                {
                    NotImplemented;
                }

                break;
            }
            case mdForce:
            {
                dictionary forceDict = coeffs_.subDict("forces");
                functionObjects::forces f("forces", mesh(), forceDict);
                f.calcForcesMoment();

                setPoint = f.forceEff().x() - T() - SFC_;
                control = controller_.executePID(setPoint, target);
                n_ += control;

                break;
            }
        }

        const_cast<radialThrustActuationDiskSource*>(this)
            ->updateThrust(Kts_*rho_*pow(n_,2)*pow(2*R_,4));

        Info<< "Rotational speed " << n_ << endl;
    }
}


void Foam::fv::radialThrustActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    bool compressible = false;
    if (eqn.dimensions() == dimForce)
    {
        compressible = true;
    }

    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        if (compressible)
        {
            addRadialActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                obr_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addRadialActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                geometricOneField(),
                U
            );
        }
    }
}

void Foam::fv::radialThrustActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    bool compressible = false;
    if (eqn.dimensions() == dimForce)
    {
        compressible = true;
    }
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        if (compressible)
        {
            addRadialActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                rho,
                U
            );
        }
        else
        {
            addRadialActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                geometricOneField(),
                U
            );
        }
    }
}

void Foam::fv::radialThrustActuationDiskSource::writeData(Ostream& os) const
{
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("propPosition", propPosition_);
        propsDict.add("propOrientation", propOrientation_);
        propsDict.regIOobject::write();
    }
    constantThrustActuationDiskSource::writeData(os);
}


bool Foam::fv::radialThrustActuationDiskSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        constantThrustActuationDiskSource::read(dict);
        diskType_ = diskModeTypeNames_.read(coeffs_.lookup("diskType"));
        coeffs_.readIfPresent("PoverD", PoverD_);
        coeffs_.readIfPresent("R", R_);
        coeffs_.readIfPresent("rh", rh_);
        coeffs_.lookup("coeffs") >> radialCoeffs_;
        coeffs_.lookup("propPosition")>> propPosition_;
        coeffs_.lookup("propOrientation")>> propOrientation_;
        shipResistanceThrust_ = coeffs_.lookupOrDefault("shipResistanceThrust", false);
        checkData();

        uShip_ = coeffs_.lookup("vesselSpeed");
        rho_ = readScalar(coeffs_.lookup("rhoInf"));

        switch (diskType_)
        {
            case mdFixed:
            {
                break;
            }
            case mdOWCurves:
            {
                coeffs_.lookup("KT") >> thrustCoeffs_;
                coeffs_.lookup("KQ") >> torqueCoeffs_;
                coeffs_.lookup("J") >> advanceCoeffs_;

                w_ = coeffs_.lookupOrDefault<scalar>("wakeFraction", 0.0);

                if (coeffs_.found("wakeFraction"))
                {
                    calcSelfPropulsionPoint();
                }

                break;
            }
        }

        SFC_ = coeffs_.lookupOrDefault<scalar>("SFC", 0.0);

        if (coeffs_.found("controller"))
        {
            tr_ = readScalar(coeffs_.subDict("controller").lookup("startTime"));
            propulsionType_ = propulsionModeTypeNames_.read(coeffs_.lookup("propulsionTarget"));

            switch (diskType_)
            {
                case mdFixed:
                {
                    n_ = readScalar(coeffs_.subDict("controller").lookup("estimatedRPS"));
                    Kts_ = T()/(rho_*pow(n_, 2)*pow(2*R_, 4));
                    break;
                }
                case mdOWCurves:
                {
                    break;
                }
            }
        }
        else
        {
            tr_ = readScalar(coeffs_.lookup("tr"));
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
