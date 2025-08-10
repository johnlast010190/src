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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::basicKinematicMPPICParcel

Description
    Definition of basic kinematic MPPIC parcel

SourceFiles
    basicKinematicMPPICParcel.H

\*---------------------------------------------------------------------------*/

#ifndef basicKinematicMPPICParcel_H
#define basicKinematicMPPICParcel_H

#include "primitives/contiguous/contiguous.H"
#include "particle/particle.H"
#include "parcels/Templates/KinematicParcel/KinematicParcel.H"
#include "parcels/Templates/MPPICParcel/MPPICParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef MPPICParcel<KinematicParcel<particle>> basicKinematicMPPICParcel;

    template<>
    constexpr bool contiguous<basicKinematicMPPICParcel>()
    {
        return true;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
