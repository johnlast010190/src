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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

inline bool Foam::radiation::absorptionCoeffs::invTemp() const
{
    return  invTemp_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Tcommon() const
{
    return  Tcommon_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Tlow() const
{
    return  Tlow_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Thigh() const
{
    return  Thigh_;
}


inline const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::highACoeffs() const
{
    return  highACoeffs_;
}


inline const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::lowACoeffs() const
{
    return  lowACoeffs_;
}


// ************************************************************************* //
