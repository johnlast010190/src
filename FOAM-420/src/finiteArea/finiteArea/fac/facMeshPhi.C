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

\*---------------------------------------------------------------------------*/

#include "finiteArea/fac/facMeshPhi.H"
#include "faMesh/faMesh.H"
#include "finiteArea/ddtSchemes/faDdtScheme/faDdtScheme.H"
#include "interpolation/edgeInterpolation/edgeInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::edgeScalarField> Foam::fac::meshPhi
(
    const areaVectorField& vf
)
{
    return fa::faDdtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + vf.name() + ')')
    ).ref().meshPhi(vf);
}


Foam::tmp<Foam::edgeScalarField> Foam::fac::meshPhi
(
    const dimensionedScalar& rho,
    const areaVectorField& vf
)
{
    return fa::faDdtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


Foam::tmp<Foam::edgeScalarField> Foam::fac::meshPhi
(
    const areaScalarField& rho,
    const areaVectorField& vf
)
{
    return fa::faDdtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


void Foam::fac::makeRelative
(
    edgeScalarField& phi,
    const areaVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= fac::meshPhi(U);
    }
}

void Foam::fac::makeRelative
(
    edgeScalarField& phi,
    const dimensionedScalar& rho,
    const areaVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= rho*fac::meshPhi(rho, U);
    }
}

void Foam::fac::makeRelative
(
    edgeScalarField& phi,
    const areaScalarField& rho,
    const areaVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= fac::interpolate(rho)*fac::meshPhi(rho, U);
    }
}


void Foam::fac::makeAbsolute
(
    edgeScalarField& phi,
    const areaVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += fac::meshPhi(U);
    }
}

void Foam::fac::makeAbsolute
(
    edgeScalarField& phi,
    const dimensionedScalar& rho,
    const areaVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += rho*fac::meshPhi(rho, U);
    }
}

void Foam::fac::makeAbsolute
(
    edgeScalarField& phi,
    const areaScalarField& rho,
    const areaVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += fac::interpolate(rho)*fac::meshPhi(rho, U);
    }
}


Foam::tmp<Foam::edgeScalarField> Foam::fac::relative
(
    const tmp<edgeScalarField>& tphi,
    const areaVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fac::meshPhi(U);
    }
    else
    {
        return tmp<edgeScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::edgeScalarField> Foam::fac::relative
(
    const tmp<edgeScalarField>& tphi,
    const areaScalarField& rho,
    const areaVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fac::interpolate(rho)*fac::meshPhi(rho, U);
    }
    else
    {
        return tmp<edgeScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::edgeScalarField> Foam::fac::absolute
(
    const tmp<edgeScalarField>& tphi,
    const areaVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi + fac::meshPhi(U);
    }
    else
    {
        return tmp<edgeScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::edgeScalarField> Foam::fac::absolute
(
    const tmp<edgeScalarField>& tphi,
    const areaScalarField& rho,
    const areaVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi + fac::interpolate(rho)*fac::meshPhi(rho, U);
    }
    else
    {
        return tmp<edgeScalarField>(tphi, true);
    }
}


// ************************************************************************* //
