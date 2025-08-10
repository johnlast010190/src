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
    (c) 2019-2021 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/CourantNumber/CourantNumber.H"
#include "finiteVolume/fvc/fvc.H"


template <class rhoType>
Foam::scalar Foam::CourantNumber
(
    const fvMesh& mesh,
    const scalar deltaT,
    const rhoType& rho,
    const surfaceScalarField& phi,
    const word description
)
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
      / rho.primitiveField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/mesh.V().field())*deltaT;

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*deltaT;

    if (mesh.name().empty())
    {
        Info<< description << " mean: " << meanCoNum
            << " max: " << CoNum << endl;
    }
    else
    {
        Info<< "Region: " << mesh.name() << " " << description
            << " mean: " << meanCoNum
            << " max: " << CoNum << endl;

    }

    return CoNum;
}


// ************************************************************************* //
