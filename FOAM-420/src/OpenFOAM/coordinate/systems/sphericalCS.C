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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "coordinate/systems/sphericalCS.H"

#include "primitives/one/one.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{
    defineTypeName(spherical);
    addToRunTimeSelectionTable(coordinateSystem, spherical, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSystem::spherical::spherical()
:
    coordinateSystem()
{}


Foam::coordSystem::spherical::spherical
(
    const coordinateSystem& csys
)
:
    coordinateSystem(csys)
{}


Foam::coordSystem::spherical::spherical(coordinateSystem&& csys)
:
    coordinateSystem(std::move(csys))
{}


Foam::coordSystem::spherical::spherical(autoPtr<coordinateSystem>&& csys)
:
    coordinateSystem(std::move(csys))
{}


Foam::coordSystem::spherical::spherical
(
    const word& name,
    const coordinateSystem& csys
)
:
    coordinateSystem(name, csys)
{}


Foam::coordSystem::spherical::spherical
(
    const point& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(origin, cr)
{}


Foam::coordSystem::spherical::spherical
(
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(word::null, origin, axis, dirn)
{}


Foam::coordSystem::spherical::spherical
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::coordSystem::spherical::spherical
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


Foam::coordSystem::spherical::spherical
(
    const dictionary& dict
)
:
    coordinateSystem(dict)
{}


Foam::coordSystem::spherical::spherical
(
    const dictionary& dict,
    const word& dictName
)
:
    coordinateSystem(dict, dictName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSystem::spherical::~spherical()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tensor Foam::coordSystem::spherical::R(const point& global) const
{
    tensor rotTensor(rot_);
    const vector ax1 = rotTensor.col<2>();
    vector ax2(global - origin_);

    // Remove colinear component
    ax2 -= ((ax1 & ax2) * ax1);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }
    ax2 /= magAxis2;  // normalise
    vector ax3 =  rotTensor.col<1>();
    rotTensor.col<1>(ax1^ax2);
    rotTensor.col<2>(ax2^ax3);
    return rotTensor;
}


Foam::vector Foam::coordSystem::spherical::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal(toCartesian(local), translate);
}


Foam::tmp<Foam::vectorField> Foam::coordSystem::spherical::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    const label len = local.size();

    auto tresult = tmp<vectorField>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] =
            coordinateSystem::localToGlobal(toCartesian(local[i]), translate);
    }

    return tresult;
}


Foam::vector Foam::coordSystem::spherical::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    return fromCartesian(coordinateSystem::globalToLocal(global, translate));
}


Foam::tmp<Foam::vectorField> Foam::coordSystem::spherical::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    const label len = global.size();

    tmp<vectorField> tresult
    (
        coordinateSystem::globalToLocal(global, translate)
    );
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = fromCartesian(result[i]);
    }

    return tresult;
}


// ************************************************************************* //
