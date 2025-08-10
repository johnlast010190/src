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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sensor/patchFieldMotionValue/patchFieldMotionValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{
    defineTypeNameAndDebug(patchFieldMotionValue, 0);

    sensor<scalar>::addmeshConstructorToTable<sensorTypes::patchFieldMotionValue>
        addpatchFieldMotionValueMeshConstructorToTable_;
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sensorTypes::patchFieldMotionValue::patchFieldMotionValue
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    sensor<scalar>(mesh, dict),
    patchName_(),
    point0_(dict.lookup("point0")),
    direction_(dict.lookup("direction")),
    velocity_(Function1<vector>::New("velocity", dict))
{
    Istream& is(dict.lookup("patchName"));
    is  >> patchName_;
}


Foam::sensorTypes::patchFieldMotionValue::patchFieldMotionValue(const patchFieldMotionValue& cv)
:
    sensor<scalar>(cv),
    patchName_(cv.patchName_),
    point0_(cv.point0_),
    direction_(cv.direction_),
    velocity_(cv.velocity_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sensorTypes::patchFieldMotionValue::~patchFieldMotionValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sensorTypes::patchFieldMotionValue::valueField() const
{
    //const volScalarField& field(mesh_.lookupObject<volScalarField>(fieldName()));
    label patchID = mesh_.boundaryMesh().findPatchID(patchName_);

    scalar t0(0.0);
    scalar t(mesh_.time().value());
    // Translation of centre of gravity with constant velocity
    const vector displacement = velocity_->integrate(t0, t);

    if (debug)
    {
        Info<< "displacement : " << displacement << endl;
    }

    tmp<scalarField> mask
    (
        new scalarField(mesh_.boundary()[patchID].size(), 0.0)
    );

    tmp<vectorField> pPts
    (
        new vectorField(mesh_.Cf().boundaryField()[patchID])
    );

    mask.ref() = pos0((pPts() - (point0_ + displacement)) & direction_);

    return mask;
}


// ************************************************************************* //
