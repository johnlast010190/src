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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "finiteVolume/fvc/fvc.H"

Foam::scalar Foam::diffusionNumber
(
    const fvMesh& mesh,
    const Time& runTime,
    const volScalarField& rho,
    const volScalarField& alpha
)
{
    scalar DiNum = 0.0;
    scalar meanDiNum = 0.0;

    surfaceScalarField alphabyRhobyDeltaSq
    (
        sqr(mesh.surfaceInterpolation::deltaCoeffs())
      * fvc::interpolate(alpha)
      / fvc::interpolate(rho)
    );

    DiNum = gMax(alphabyRhobyDeltaSq.primitiveField())*runTime.deltaT().value();

    meanDiNum = (average(alphabyRhobyDeltaSq)).value()*runTime.deltaT().value();

    if (mesh.name().empty())
    {
        Info<< "Diffusion Number mean: " << meanDiNum
            << " max: " << DiNum << endl;
    }
    else
    {
        Info<< "Region: " << mesh.name()
            << " Diffusion Number mean: " << meanDiNum
            << " max: " << DiNum << endl;
    }

    return DiNum;
}

// ************************************************************************* //
