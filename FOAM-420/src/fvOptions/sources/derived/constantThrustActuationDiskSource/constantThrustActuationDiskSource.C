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
    (c) 2013 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "constantThrustActuationDiskSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(constantThrustActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        constantThrustActuationDiskSource,
        dictionary
    );
}

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::constantThrustActuationDiskSource::volumeModeType,
        2
    >::names[] =
    {
        "absolute",
        "specific"
    };
}

const Foam::NamedEnum<Foam::fv::constantThrustActuationDiskSource::volumeModeType, 2>
Foam::fv::constantThrustActuationDiskSource::volumeModeTypeNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::constantThrustActuationDiskSource::checkData() const
{
    if (T_ <= VSMALL)
    {
        FatalErrorInFunction
           << "T must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorInFunction
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }
}


void Foam::fv::constantThrustActuationDiskSource::createCS()
{
    if (coeffs_.found("referenceFrame"))
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, coeffs_);
    }
    else if (coeffs_.found("coordinateSystem"))
    {
        coorSysPtr_ = coordinateSystem::New(mesh_, coeffs_, coordinateSystem::typeName_());
    }
    else
    {
        coorSysPtr_.reset(new coordSystem::cartesian());
    }

    if (debug)
    {
        Info<< "    Rotor gometry:" << nl
            << "    - thrust        = " << T_ << nl
            << "    - disk dir      = " << diskDir_ << nl
            << "    - origin        = " << coorSysPtr_().origin() << nl
            << "    - r-axis        = " << coorSysPtr_().e1() << nl
            << "    - psi-axis      = " << coorSysPtr_().e2() << nl
            << "    - z-axis        = " << coorSysPtr_().e3() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::constantThrustActuationDiskSource::constantThrustActuationDiskSource
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
    diskDir_(vector::zero),
    TF1_(nullptr),
    T_(0),
    coorSysPtr_(),
    coorFramePtr_(nullptr)
{
    read(dict);
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation disk zone: "
        << this->name() << endl;

    IFstream propertiesFile
    (
        mesh_.time().timePath()/"uniform"/(this->name() + "Properties")
    );

    if (propertiesFile.good())
    {
        Info<< "    Reading thrust and diskDir from file" << endl;
        dictionary propertiesDict(dictionary::null, propertiesFile);

        readThrustData(propertiesDict);
    }

    createCS();
    checkData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::constantThrustActuationDiskSource::addSup
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

    if (V() > VSMALL)
    {
        if (compressible)
        {
            addActuationDiskAxialInertialResistance
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
            addActuationDiskAxialInertialResistance
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


void Foam::fv::constantThrustActuationDiskSource::addSup
(
    const volScalarField& rho,
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

    if (V() > VSMALL)
    {
        if (compressible)
        {
            addActuationDiskAxialInertialResistance
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
            addActuationDiskAxialInertialResistance
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


void Foam::fv::constantThrustActuationDiskSource::writeData(Ostream& os) const
{
    Info<< "constantThrustActuationDiskSource::writeData" << endl;
    os  << indent << name_ << endl;
    dict_.write(os);
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

        if (coeffs_.found("diskDir"))
        {
            propsDict.add("diskDir", diskDir_);
        }
        else
        {
            propsDict.add("thrustDir", diskDir_);
        }

        OStringStream out;
        OSstream& out2 = dynamic_cast<OSstream&>(out);
        TF1_->writeData(out2);
        IStringStream in(word(out.str()));
        dictionary newdict(in);
        propsDict.merge(newdict);

        propsDict.regIOobject::write();
    }
}


bool Foam::fv::constantThrustActuationDiskSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        // set a default for interface backward compatibility
        volumeMode_ = vmAbsolute;
        if (coeffs_.found("volumeMode"))
        {
            volumeMode_ = volumeModeTypeNames_.read
            (
                coeffs_.lookup("volumeMode")
            );
        }

        // Set volume normalisation
        if (volumeMode_ == vmAbsolute)
        {
            VDash_ = V_;
        }

        readThrustData(coeffs_);
        checkData();

        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::fv::constantThrustActuationDiskSource::readThrustData(const dictionary& dict)
{
    dict.readIfPresent("diskDir", diskDir_);
    dict.readIfPresent("thrustDir", diskDir_);

    scalar t = mesh_.time().value();
    if (dict.found("T"))
    {
        TF1_ = Function1<scalar>::New("T", dict);
        T_ = TF1_->value(t);
    }
    else if (dict.found("thrust"))
    {
        TF1_ = Function1<scalar>::New("thrust", dict);
        T_ = TF1_->value(t);
    }
    else
    {
        WarningInFunction << "No thrust input specified." << endl;
    }

    return true;
}

bool Foam::fv::constantThrustActuationDiskSource::updateThrust(const scalar& thrust)
{
    T_ = thrust;
    return true;
}

bool Foam::fv::constantThrustActuationDiskSource::updateDiskDirection(const vector& diskDirection)
{
    diskDir_ = diskDirection;
    return true;
}

inline const Foam::coordinateSystem&
Foam::fv::constantThrustActuationDiskSource::csys() const
{
    if (coorFramePtr_)
    {
        return coorFramePtr_->coorSys();
    }
    return *coorSysPtr_;
}


// ************************************************************************* //
