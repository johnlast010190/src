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
    (c) 2010-2011, Esi Ltd

Application
    heatFlux

Description

    Default behaviour assumes operating in incompressible mode. To apply to
    compressible RAS cases, use the -compressible option.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "psiThermo/psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressibleQ
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    const volScalarField& T,
    volVectorField& Q
)
{
    #include "cfdTools/incompressible/createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turb
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    Q = -turb->Cp() * turb->alphaEff() * fvc::grad(T);
}


void calcCompressibleQ
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& Q
)
{

    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.typeHeaderOk<volScalarField>(true))
    {
        Info<< "    no rho field" << endl;
        return;
    }

    Info<< "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

    #include "cfdTools/compressible/compressibleCreatePhi.H"

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();

    autoPtr<compressible::turbulenceModel> turb
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Q = -turb->alphaEff() * fvc::grad(thermo.he());
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "include/addRegionOption.H"

    argList::validOptions.insert("compressible","");

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "include/createNamedMesh.H"

    bool compressible = args.optionFound("compressible");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volVectorField Q
        (
            IOobject
            (
                "q",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "Q",
                dimensionSet(1, 0 , -3, 0, 0, 0, 0) ,
                vector::zero
            )
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        IOobject THeader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            UHeader.typeHeaderOk<volVectorField>(true)
            && THeader.typeHeaderOk<volScalarField>(true)
        )
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);
            volScalarField T(THeader, mesh);

            if (compressible)
            {
                calcCompressibleQ(mesh, runTime, U, Q);
            }
            else
            {
                calcIncompressibleQ(mesh, runTime, U, T, Q);
            }
        }
        else
        {
            Info<< "    no U or T field" << endl;
        }

        Info<< "Writing heat flux to field " << Q.name() << nl << endl;

        Q.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
