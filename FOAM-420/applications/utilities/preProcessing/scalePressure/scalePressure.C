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
    (c) 2019 Esi Ltd.

Application
    scalePressure

Description
    scale single phase solver pressure to multiphase solver pressure
    (multiply with density and add offset value)

Author
    Daniel Deising

-------------------------------------------------------------------------------
*/

#include "cfdTools/general/include/fvCFD.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "include/setRootCase.H"

#include "include/createTime.H"
#include "include/createMesh.H"

    Info<< "scaling pressure ..." << endl;

    IOobject pHeader
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading p" << endl;
    volScalarField p(pHeader, mesh);

    IOdictionary scalePressureDict
    (
        IOobject
        (
            "scalePressureDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho
    (
        "rho", scalePressureDict.lookup("rho")
    );

    dimensionedScalar pRef
    (
        "pRef", dimensionSet(1,-1,-2,0,0), readScalar(scalePressureDict.lookup("pRef"))
    );

    wordList pBoundaryTypes = p.boundaryField().types();

    Info<< "Creating p_rgh" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p * rho,
        pBoundaryTypes
    );

    p_rgh.forceAssign(p_rgh + pRef);

    p_rgh.correctBoundaryConditions();

    p_rgh.write();

    Info<< endl;

    return(0);
}


// ************************************************************************* //
