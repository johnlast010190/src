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
    (c) 2014 Esi Ltd.

Application
    scaleUprime

Description
    Scale velocity as follows:
        U = UMean + (U-UMean)*(1+alpha);

    For permuting DES flow to study the effect of stochastics on the results

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "scaleUprime: scale perturbation velocity using: "
        "U = UMean + (U-UMean)*(1+alpha)"
    );
    argList::validArgs.append("scaling constant");

#include "include/setRootCase.H"

#include "include/createTime.H"
#include "include/createMesh.H"

    IOobject UHeader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading U" << endl;
    volVectorField U(UHeader, mesh);

    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading UMean" << endl;
    volVectorField UMean(UMeanHeader, mesh);

    const scalar alpha = args.argRead<scalar>(1);

    Info<< "Scaling Uprime." << endl;
    U = UMean + (U - UMean) * (1 + alpha);

    Info<< "Writing modified U to file." << endl;
    U.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
