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

Application
    engineSwirl

Group
    grpPreProcessingUtilities

Description
    Generates a swirling flow for engine calulations.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar Vphi = (constant::mathematical::pi*swirlRPMRatio*rpm/30).value();
    scalar b1 = j1(swirlProfile).value();
    scalar b2 = 2.0*b1/swirlProfile.value() - j0(swirlProfile).value();

    scalar omega = 0.125*(Vphi*bore*swirlProfile/b2).value();

    scalar cylinderRadius = 0.5*bore.value();

    scalar Umax = 0.0;
    forAll(mesh.C(), celli)
    {
        vector c = mesh.C()[celli] - swirlCenter;
        scalar r = ::pow(sqr(c & xT) + sqr(c & yT), 0.5);

        if (r <= cylinderRadius)
        {
            scalar b = j1(swirlProfile*r/cylinderRadius).value();
            scalar vEff = omega*b;
            r = max(r, SMALL);
            U[celli] = ((vEff/r)*(c & yT))*xT + (-(vEff/r)*(c & xT))*yT;
            Umax = max(Umax, mag(U[celli]));
        }
    }

    Info<< "Umax = " << Umax << endl;

    U.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
