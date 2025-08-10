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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd

\*---------------------------------------------------------------------------*/


inline const Foam::volScalarField::Boundary&
Foam::radiation::radiativeIntensityRay::qr() const
{
    return qr_();
}


inline Foam::volScalarField::Boundary& Foam::radiation::radiativeIntensityRay::qr()
{
    return qr_();
}

inline const Foam::volScalarField::Boundary& Foam::radiation::
radiativeIntensityRay::qin() const
{
    return qin_();
}


inline Foam::volScalarField::Boundary& Foam::radiation::radiativeIntensityRay::qin()
{
    return qin_();
}


inline const Foam::volScalarField::Boundary& Foam::radiation::
radiativeIntensityRay::qem() const
{
    return qem_();
}


inline Foam::volScalarField::Boundary& Foam::radiation::radiativeIntensityRay::qem()
{
    return qem_();
}

inline const Foam::volScalarField::Boundary& Foam::radiation::
radiativeIntensityRay::qg() const
{
    return qg_();
}


inline Foam::volScalarField::Boundary& Foam::radiation::radiativeIntensityRay::qg()
{
    return qg_();
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRay::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRay::dAve() const
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::omega() const
{
    return omega_;
}

inline Foam::word Foam::radiation::radiativeIntensityRay::rayName() const
{
    return word("I" + name(myRayId_));
}

inline const Foam::volScalarField&
Foam::radiation::radiativeIntensityRay::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
