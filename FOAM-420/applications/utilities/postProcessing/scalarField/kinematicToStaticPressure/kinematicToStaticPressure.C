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
    (c) 2017 Esi Ltd.

Application
    incompressiblePressure

Description
    Calculates the pressure, P, from the kinematic pressure, p
    Also write constant rho field (for Actran FFT input)

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/turbulentTransportModels/turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("Reference pressure");
    #include "include/addOverwriteOption.H"
    #include "include/addRegionOption.H"
    Foam::argList::addBoolOption
    (
        "volumeFlux",
        "Transform mass to volume flux"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "include/createNamedMesh.H"

    const bool overwrite = args.optionFound("overwrite");
    const bool transformFlux = args.optionFound("volumeFlux");
    const scalar refPressure = args.argRead<scalar>(1);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        Info<< "Reading field p\n" << endl;
        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        #include "cfdTools/incompressible/createPhi.H"

        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
            incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );

        volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            turbulence->rho()()
        );

        volScalarField P
        (
            IOobject
            (
                "pstatic",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            p
        );

        P.dimensions().reset(p.dimensions()*rho.dimensions());

        P.primitiveFieldRef() =
            refPressure + p.primitiveField()*rho.primitiveField();

        forAll(P.boundaryField(), patchi)
        {
            P.boundaryFieldRef()[patchi].fvPatchScalarField::operator=
            (
                P.boundaryField()[patchi]*rho.boundaryField()[patchi]
              + refPressure
            );
        }

        if (transformFlux)
        {
            surfaceScalarField Phi
            (
                IOobject
                (
                    "Phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                phi
            );
            Phi.dimensions().reset(phi.dimensions()*rho.dimensions());
            Phi.primitiveFieldRef() = phi.primitiveField()*rho.primitiveField();
            if (overwrite)
            {
                phi.checkOut();
                Phi.rename("phi");
            }
            Phi.write();
            Info<< "Writing flux to file ("
                << Phi.name() <<")"
                << nl << endl;
        }

        if (overwrite)
        {
            p.checkOut();
            P.rename("p");
        }

        Info<< "Writing pressure and density to file ("
            << P.name() << ", " << rho.name() << ")"
            << nl << endl;

        P.write();
        rho.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
