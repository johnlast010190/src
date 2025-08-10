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

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::point& Foam::kinematicParcelInjectionData::x() const
{
    return x_;
}


inline const Foam::vector& Foam::kinematicParcelInjectionData::U() const
{
    return U_;
}


inline Foam::scalar Foam::kinematicParcelInjectionData::d() const
{
    return d_;
}


inline Foam::scalar Foam::kinematicParcelInjectionData::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::kinematicParcelInjectionData::mDot() const
{
    return mDot_;
}


inline Foam::point& Foam::kinematicParcelInjectionData::x()
{
    return x_;
}


inline Foam::vector& Foam::kinematicParcelInjectionData::U()
{
    return U_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionData::d()
{
    return d_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionData::rho()
{
    return rho_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionData::mDot()
{
    return mDot_;
}


// ************************************************************************* //
