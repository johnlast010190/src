/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coordinate/systems/cylindricalCS.H"
#include "coordinate/rotation/cylindricalRotation.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{
    defineTypeName(cylindrical);
    addToRunTimeSelectionTable(coordinateSystem, cylindrical, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSystem::cylindrical::cylindrical()
:
    coordinateSystem()
{}


Foam::coordSystem::cylindrical::cylindrical(const coordinateSystem& csys)
:
    coordinateSystem(csys)
{}


Foam::coordSystem::cylindrical::cylindrical(coordinateSystem&& csys)
:
    coordinateSystem(std::move(csys))
{}


Foam::coordSystem::cylindrical::cylindrical(autoPtr<coordinateSystem>&& csys)
:
    coordinateSystem(std::move(csys))
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const point& origin,
    const coordinateRotation& crot
)
:
    coordinateSystem(origin, crot)
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const point& origin,
    const vector& axis
)
:
    cylindrical(word::null, origin, axis)
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const word& name,
    const point& origin,
    const vector& axis
)
:
    coordinateSystem
    (
        name,
        origin,
        coordinateRotations::cylindrical(axis)
    )
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    cylindrical(word::null, origin, axis, dirn)
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


Foam::coordSystem::cylindrical::cylindrical(const dictionary& dict)
:
    coordinateSystem(dict)
{}


Foam::coordSystem::cylindrical::cylindrical
(
    const dictionary& dict,
    const word& dictName
)
:
    coordinateSystem(dict, dictName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tensor Foam::coordSystem::cylindrical::R(const point& global) const
{
    // Robuster version of coordinateRotations::axes::rotation()
    // using an E3_E1 order and falling back to the top-level rotation
    // tensor if the directional input is borderline.

    tensor rotTensor(rot_);

    const vector ax1 = rotTensor.col<2>(); // == e3 (already normalized)

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

    // Replace with updated local axes

    rotTensor.col<0>(ax2);
    rotTensor.col<1>(ax1^ax2);

    return rotTensor;
}


Foam::vector Foam::coordSystem::cylindrical::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal(toCartesian(local), translate);
}


Foam::tmp<Foam::vectorField> Foam::coordSystem::cylindrical::localToGlobal
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
            coordinateSystem::localToGlobal
            (
                toCartesian(local[i]),
                translate
            );
    }

    return tresult;
}


Foam::vector Foam::coordSystem::cylindrical::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    return fromCartesian(coordinateSystem::globalToLocal(global, translate));
}


Foam::tmp<Foam::vectorField> Foam::coordSystem::cylindrical::globalToLocal
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
