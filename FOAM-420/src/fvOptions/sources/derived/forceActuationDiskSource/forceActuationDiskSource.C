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
    (c) 2022 Esi Ltd.

\*----------------------------------------------------------------------------*/

#include "forceActuationDiskSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(forceActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        forceActuationDiskSource,
        dictionary
    );
}

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::forceActuationDiskSource::volumeModeType,
        2
    >::names[] =
    {
        "absolute",
        "specific"
    };
}

const Foam::NamedEnum<Foam::fv::forceActuationDiskSource::volumeModeType, 2>
Foam::fv::forceActuationDiskSource::volumeModeTypeNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::forceActuationDiskSource::createCS()
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::forceActuationDiskSource::forceActuationDiskSource
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
    force_(nullptr),
    coorSysPtr_(nullptr),
    coorFramePtr_(nullptr),
    radialDistribution_(true)
{
    read(dict);
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating force actuation disk zone: "
        << this->name() << endl;

    createCS();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::forceActuationDiskSource::addSup
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
            addForceActuationDiskResistance
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
            addForceActuationDiskResistance
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


void Foam::fv::forceActuationDiskSource::addSup
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
            addForceActuationDiskResistance
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
            addForceActuationDiskResistance
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


void Foam::fv::forceActuationDiskSource::writeData(Ostream& os) const
{
}


bool Foam::fv::forceActuationDiskSource::read(const dictionary& dict)
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
        radialDistribution_ =
            coeffs_.lookupOrDefault<Switch>("radialDistribution", true);

        force_ = Function1<vector>::New("force", coeffs_);

        return true;
    }
    else
    {
        return false;
    }
}


inline const Foam::coordinateSystem& Foam::fv::forceActuationDiskSource::csys() const
{
    if (coorFramePtr_)
    {
        return coorFramePtr_->coorSys();
    }
    return *coorSysPtr_;
}


// ************************************************************************* //
