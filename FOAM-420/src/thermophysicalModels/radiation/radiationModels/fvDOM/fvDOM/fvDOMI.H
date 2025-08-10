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
    (c) 2010-2017, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

inline const Foam::radiation::radiativeIntensityRay&
Foam::radiation::fvDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::fvDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::fvDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::fvDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::fvDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::fvDOM::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::qr() const
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiation::fvDOM::qin() const
{
    return qin_;
}

/*
inline const Foam::volScalarField& Foam::radiation::fvDOM::qem() const
{
    return qem_;
}
*/

inline const Foam::volScalarField& Foam::radiation::fvDOM::qg() const
{
    return qg_;
}


inline const Foam::radiation::blackBodyEmission&
Foam::radiation::fvDOM::blackBody() const
{
    return blackBody_;
}


inline const Foam::waveSolverCoeffs& Foam::radiation::fvDOM::rayWaveCoeffs
(
    const label rayId,
    const label lambdaI
) const
{
    return rayWaveCoeffs_[lambdaI][rayId];
}

inline Foam::scalar Foam::radiation::fvDOM::omegaMax() const
{
    return omegaMax_;
}


inline Foam::vector Foam::radiation::fvDOM::meshOrientation() const
{
    return meshOrientation_;
}

inline const Foam::radiation::boundaryRadiationProperties&
Foam::radiation::fvDOM::boundaryProperties() const
{
    return boundaryRadiationProperties::New(mesh_);
}

// ************************************************************************* //
