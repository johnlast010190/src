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
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "tmp.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::tryReregister(regIOobject* newObj, regIOobject* oldObj)
{
    // Swap new and old object registrations in case where old obj would have
    // blocked the new one

    if (!newObj || isNull(*newObj) || !oldObj || isNull(*oldObj))
    {
        // Pointers can be non-null but point to a reinterpreted 'null object'
        // which is not actually an IOobject at all and should not be accessed
        return;
    }

    if
    (
        newObj
     && oldObj
     && newObj->name() == oldObj->name()
     && newObj->registerObject()
     && !newObj->isRegistered()
     && oldObj->isRegistered()
    )
    {
        oldObj->checkOut();
        newObj->checkIn();
    }
}

// ************************************************************************* //
