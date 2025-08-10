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

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sensor/timeValue/timeValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{
    defineTypeNameAndDebug(timeValue, 0);

    sensor<scalar>::addmeshConstructorToTable<sensorTypes::timeValue>
        addtimeValueMeshConstructorToTable_;
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sensorTypes::timeValue::timeValue
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    sensor<scalar>(mesh, dict),
    patchName_()
{
    Istream& is(dict.lookup("patchName"));
    is  >> patchName_;
}


Foam::sensorTypes::timeValue::timeValue(const timeValue& cv)
:
    sensor<scalar>(cv),
    patchName_(cv.patchName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sensorTypes::timeValue::~timeValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sensorTypes::timeValue::valueField() const
{
    label patchID = mesh_.boundaryMesh().findPatchID(patchName_);
    scalar timeValue = mesh_.time().value();

    return tmp<scalarField>
    (
        new scalarField(mesh_.boundaryMesh()[patchID].size(), timeValue)
    );
}


Foam::scalar Foam::sensorTypes::timeValue::value() const
{
    return mesh_.time().value();
}

// ************************************************************************* //
