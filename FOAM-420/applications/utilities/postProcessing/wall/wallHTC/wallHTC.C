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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Application
    wallHTC

Group
    grpPostProcessingUtilities

Description
    Calculates and writes the Heat Transfer Coefficient based on a reference
    temperature (Tref) and the heat flux for all patches. Tref can be set
    in a dictionary (e.g.system/default/wallhtcDict). Otherwise a default
    value of 298.15 K will be used.
    Selecting -nusselt it also plots the local Nusselt number. kref and lref
    have to be specified in the dictionary.
    e.g. wallHTC -latestTime -compressible -nusselt -region 'default'
\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "psiThermo/psiThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "rhoThermo/rhoThermo.H"
#include "fluidThermo/fluidThermo.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "regionProperties/regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void zeroToSmall(scalarField& field)
{
    forAll(field, i)
    {
        if (mag(field[i]) < SMALL)
        {
            field[i] = SMALL;
        }
    }
}

void calcIncompressibleQ
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    const volScalarField& Qr,
    volScalarField& Q,
    autoPtr<volScalarField>& totalWallHeatFlux,
    volScalarField& wallHTCconv,
    autoPtr<volScalarField>& totalWallHTC,
    bool printNu,
    volScalarField& wallNu,
    autoPtr<volScalarField>& totalwallNu,
    bool isRegion,
    const word& regionName
)
{
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    //reading the reference temperature from a dict
    IOdictionary wallhtcDict
    (
       IOobject
       (
            "wallhtcDict",
            runTime.system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
       )
    );

    Info<< "Reading reference temperature ..." << endl;

    const scalar Tref = wallhtcDict.lookupOrDefault<scalar>("Tref", 298.15);
    if (printNu)
    {
        Info<< "Reading reference length ..." << nl
            << "Reading fluid thermal conductivity ..." << endl;
    }

    const scalar lref = wallhtcDict.lookupOrDefault<scalar>("lref", 1.0);
    const scalar kref = wallhtcDict.lookupOrDefault<scalar>("kref", 0.02624);

    #include "cfdTools/incompressible/createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turb
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    surfaceScalarField heatFlux
    (
        fvc::interpolate(turb->alphaEff()*turb->rho()*turb->Cp())
       *fvc::snGrad(T)
    );

    const surfaceScalarField::Boundary& patchHeatFlux =
        heatFlux.boundaryField();

    const volScalarField::Boundary& patchRadHeatFlux =
        Qr.boundaryField();


    const surfaceScalarField::Boundary& magSf =
        mesh.magSf().boundaryField();

    const volScalarField::Boundary& patchTempField =
        T.boundaryField();

    Info<< "Reference Temperature      [K]: " << Tref << endl;
    if (printNu)
    {
        Info<< "Reference Length           [m]: " << lref  << nl
            << "Reference conductiviy  [W/m K]: " << kref  << endl;
    }

    Info<< "\n Wall Report" << endl;
    forAll(patchHeatFlux, patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            const scalar convFlux = gSum(magSf[patchi]*patchHeatFlux[patchi]);

            const scalar radFlux =
                -gSum(magSf[patchi]*patchRadHeatFlux[patchi]);

            const scalar wallTemp =
                (gSum(magSf[patchi]*patchTempField[patchi]))
               /gSum(magSf[patchi]);

            const scalar DTemp = (Tref == wallTemp) ? SMALL : (wallTemp-Tref);

            const scalar htc_ref =
                (
                    gSum(magSf[patchi]*patchHeatFlux[patchi])
                   /gSum(magSf[patchi])
                )
               /DTemp;

            Info<< mesh.boundary()[patchi].name() << nl
                << "  convective heat flux [W]: " << convFlux << nl
                << "  radiative  heat flux [W]: " << radFlux << nl
                << "  total heat flux      [W]: " << convFlux + radFlux << nl
                << "  wall temp            [K]: " << wallTemp << nl
                << "  Tw-Tref              [K]: " << DTemp << nl
                << "  conv. htc_ref   [W/m2 K]: " << htc_ref << endl;

            if (printNu)
            {
                const scalar Nusselt = htc_ref*lref/kref;
                Info<< "  Nusselt convection   [-]: " << Nusselt << endl;
            }
        }
    }

    Info<< endl;

    Q.dimensions().reset(heatFlux.dimensions());

    forAll(Q.boundaryField(), patchi)
    {
        Q.boundaryFieldRef()[patchi] = patchHeatFlux[patchi];
    }

    if (totalWallHeatFlux.valid())
    {
        totalWallHeatFlux->dimensions().reset(heatFlux.dimensions());
        forAll(totalWallHeatFlux->boundaryField(), patchi)
        {
            totalWallHeatFlux->boundaryFieldRef()[patchi] =
                patchHeatFlux[patchi] - patchRadHeatFlux[patchi];
        }
    }

    //HTC given by heatflux/(Tw-Tref)
    wallHTCconv.dimensions().reset(heatFlux.dimensions());
    forAll(wallHTCconv.boundaryField(), patchi)
    {
        scalarField dT
        (
            patchTempField[patchi] - Tref
        );

        zeroToSmall(dT);

        wallHTCconv.boundaryFieldRef()[patchi] =
             (patchHeatFlux[patchi]/dT);
    }

    if (totalWallHTC.valid())
    {
        totalWallHTC->dimensions().reset(heatFlux.dimensions());
        forAll(totalWallHTC->boundaryField(), patchi)
        {
            scalarField dT
            (
                patchTempField[patchi] - Tref
            );

            zeroToSmall(dT);

            totalWallHTC->boundaryFieldRef()[patchi] =
                (patchHeatFlux[patchi] - patchRadHeatFlux[patchi])
               /dT;
        }
    }

    //local Nusselt number based on Tref
    if (printNu)
    {
        wallNu.dimensions().reset(heatFlux.dimensions());
        forAll(wallNu.boundaryField(), patchi)
        {
            scalarField dT
            (
                patchTempField[patchi] - Tref
            );

            zeroToSmall(dT);

            wallNu.boundaryFieldRef()[patchi] =
                (patchHeatFlux[patchi]/dT)*lref/kref;
        }
        if (totalwallNu.valid())
        {
            totalwallNu->dimensions().reset(heatFlux.dimensions());
            forAll(totalwallNu->boundaryField(), patchi)
            {
                scalarField dT
                (
                    patchTempField[patchi] - Tref
                );

                zeroToSmall(dT);

                totalwallNu->boundaryFieldRef()[patchi] =
                    (
                        (patchHeatFlux[patchi] - patchRadHeatFlux[patchi])
                       /dT
                    )
                   *lref/kref;
            }
        }
    }
}


void calcCompressibleQ
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    const volScalarField& Qr,
    volScalarField& Q,
    autoPtr<volScalarField>& totalWallHeatFlux,
    volScalarField& wallHTCconv,
    autoPtr<volScalarField>& totalWallHTC,
    bool printNu,
    volScalarField& wallNu,
    autoPtr<volScalarField>& totalwallNu,
    bool isRegion,
    const word& regionName

)
{
    word thermoType("none");
    {
        IOobject thermoHeader
        (
            basicThermo::dictName,
            isRegion
          ? runTime.caseConstant()/regionName
          : runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        IOobject materialsHeader
        (
            basicThermo::matDictName,
            isRegion
          ? runTime.caseSystem()/regionName
          : runTime.caseSystem(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        const bool isMaterials(materialsHeader.typeHeaderOk<IOdictionary>(true));

        IOdictionary thermoDict(isMaterials ? materialsHeader : thermoHeader);

        word tt = "none";
        if (isMaterials)
        {
            tt = word(thermoDict.lookup("materialType"));
            if (tt == "fluid")
            {
                thermoType = "rho";
                Info<< "Detected 'rho'-type material models." << endl;
            }
            else if (tt == "psiFluid")
            {
                thermoType = "psi";
                Info<< "Detected 'psi'-type material models." << endl;
            }
        }
        else
        {
            // Now it will crush in case it doesn't find it as a materialType or subdict
            tt = word(thermoDict.subDict("thermoType").lookup("type"));
            if (tt.find("Rho", 0) != string::npos)
            {
                thermoType = "rho";
                Info<< "Detected 'rho'-type thermo." << endl;
            }
            else if (tt.find("Psi", 0) != string::npos)
            {
                thermoType = "psi";
                Info<< "Detected 'psi'-type thermo." << endl;
            }
        }
        if (thermoType == "none")
        {
            thermoType = tt;
        }
    }

    // The temperature field has to be loaded before thermo
    // otherwise the thermo will use the old temperature field
    // from object registry
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Clean-up material models sub-registry for new timestep
    if (mesh.thisDb().found("materialModels"))
    {
        mesh.thisDb().checkOut("materialModels");
    }

    autoPtr<fluidThermo> thermo;

    if (thermoType == "psi")
    {
        thermo.reset(psiThermo::New(mesh).ptr());
    }
    else if (thermoType == "rho")
    {
        thermo.reset(rhoThermo::New(mesh).ptr());
    }
    else
    {
        FatalError << "Unsupported thermo type: " << thermoType
                   << exit(FatalError);
    }

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        thermo->rho()
    );

    #include "cfdTools/compressible/compressibleCreatePhi.H"

    autoPtr<compressible::turbulenceModel> turb
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo()
        )
    );

    //    reading the reference temperature from a dict
    IOdictionary wallhtcDict
    (
       IOobject
       (
            "wallhtcDict",
            runTime.system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
       )
    );

    Info<< "Reading reference temperature ..." << endl;

    const scalar Tref = wallhtcDict.lookupOrDefault<scalar>("Tref", 298.15);
    if (printNu)
    {
        Info<< "Reading reference length ..." << nl
            << "Reading fluid thermal conductivity ..." << endl;
    }
    const scalar lref = wallhtcDict.lookupOrDefault<scalar>("lref", 1.0);
    const scalar kref = wallhtcDict.lookupOrDefault<scalar>("kref", 0.02624);

    surfaceScalarField heatFlux
    (
        fvc::interpolate(turb->alphaEff())*fvc::snGrad(thermo->he())
    );

    const surfaceScalarField::Boundary& patchHeatFlux =
        heatFlux.boundaryField();

    const volScalarField::Boundary& patchRadHeatFlux =
        Qr.boundaryField();

    const surfaceScalarField::Boundary& magSf =
        mesh.magSf().boundaryField();

    const volScalarField::Boundary& patchTempField =
        T.boundaryField();

    Info<< "Reference Temperature      [K]: " << Tref  << endl;
    if (printNu)
    {
        Info<< "Reference Length           [m]: " << lref  << nl
            << "Reference Conductiviy  [W/m K]: " << kref  << endl;
    }

    Info<< "\n Wall Report" << endl;
    forAll(patchHeatFlux, patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            const scalar convFlux = gSum(magSf[patchi]*patchHeatFlux[patchi]);

            const scalar radFlux =
                -gSum(magSf[patchi]*patchRadHeatFlux[patchi]);

            const scalar wallTemp =
                (gSum(magSf[patchi]*patchTempField[patchi]))
               /gSum(magSf[patchi]);

            const scalar DTemp = (Tref == wallTemp) ? SMALL : (wallTemp-Tref);

            const scalar htc_ref =
                (gSum(magSf[patchi]*patchHeatFlux[patchi])
               /gSum(magSf[patchi]))/DTemp;

            Info<< mesh.boundary()[patchi].name() << nl
                << "  convective heat flux [W]: " << convFlux << nl
                << "  radiative  heat flux [W]: " << radFlux << nl
                << "  total heat flux      [W]: " << convFlux + radFlux << nl
                << "  wall temp            [K]: " << wallTemp << nl
                << "  Tw-Tref              [K]: " << DTemp << nl
                << "  conv.  htc_ref  [W/m2 K]: " << htc_ref << endl;

            if (printNu)
            {
                const scalar Nusselt = htc_ref*lref/kref;
                Info<< "  Nusselt convection   [-]: " << Nusselt << endl;
            }
        }
    }

    Info<< endl;

    Q.dimensions().reset(heatFlux.dimensions());

    forAll(Q.boundaryField(), patchi)
    {
        Q.boundaryFieldRef()[patchi] = patchHeatFlux[patchi];
    }
    if (totalWallHeatFlux.valid())
    {
        totalWallHeatFlux->dimensions().reset(heatFlux.dimensions());
        forAll(totalWallHeatFlux->boundaryField(), patchi)
        {
            totalWallHeatFlux->boundaryFieldRef()[patchi] =
                patchHeatFlux[patchi] - patchRadHeatFlux[patchi];
        }
    }

    //HTC given by (heatflux/(Tw-Tref)
    wallHTCconv.dimensions().reset(heatFlux.dimensions());
    forAll(wallHTCconv.boundaryField(), patchi)
    {
        scalarField dT
        (
            patchTempField[patchi] - Tref
        );

        zeroToSmall(dT);

        wallHTCconv.boundaryFieldRef()[patchi] =
            (patchHeatFlux[patchi]/dT);
    }
    if (totalWallHTC.valid())
    {
        totalWallHTC->dimensions().reset(heatFlux.dimensions());
        forAll(totalWallHTC->boundaryField(), patchi)
        {
            scalarField dT
            (
                patchTempField[patchi] - Tref
            );

            zeroToSmall(dT);

            totalWallHTC->boundaryFieldRef()[patchi] =
                (patchHeatFlux[patchi] - patchRadHeatFlux[patchi])
                /dT;
        }
    }

    //local Nusselt number based on Tref
    if (printNu)
    {
        wallNu.dimensions().reset(heatFlux.dimensions());
        forAll(wallNu.boundaryField(), patchi)
        {
            scalarField dT
            (
                patchTempField[patchi] - Tref
            );

            zeroToSmall(dT);

            wallNu.boundaryFieldRef()[patchi] =
                (patchHeatFlux[patchi]/dT)*lref/kref;
        }
        if (totalwallNu.valid())
        {
            totalwallNu->dimensions().reset(heatFlux.dimensions());
            forAll(totalwallNu->boundaryField(), patchi)
            {
                scalarField dT
                (
                    patchTempField[patchi] - Tref
                );

                zeroToSmall(dT);

                totalwallNu->boundaryFieldRef()[patchi] =
                    (
                        (patchHeatFlux[patchi] - patchRadHeatFlux[patchi])
                       /dT
                    )
                   *lref/kref;
            }
        }
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validOptions.insert("compressible","");
    argList::validOptions.insert("nusselt","");
    #include "include/addRegionOption.H"
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "include/createNamedMesh.H"

    // Load all the region meshes and make them available through object
    // registry. This is needed for new region coupled boundary conditions
    PtrList<fvMesh> regionMeshes;
    if
    (
        IOobject
        (
            "regionProperties",
            runTime.caseConstant(),
            runTime.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ).typeHeaderOk<dictionary>(true)
    )
    {
        const wordList regionNames(regionProperties(runTime).allRegionNames());
        forAll(regionNames, regioni)
        {
            const word& regionNameI = regionNames[regioni];
            if (regionNameI != regionName)
            {
                regionMeshes.append
                (
                    new fvMesh
                    (
                        IOobject
                        (
                            regionNameI,
                            runTime.timeName(),
                            runTime,
                            IOobject::MUST_READ
                        )
                    )
                );
            }
        }
    }

    //force compressible
    bool compressible = args.optionFound("compressible");
    //force Nusselt
    bool printNu = args.optionFound("nusselt");
    bool isRegion = args.optionFound("region");

    if (printNu)
    {
        Info<< " Additional output: Nusselt Number " << endl;
    }

    // if incompressible
    if (!compressible)
    {
        //if compressible not forced, check if transportProperties exists
        IOobject thermoHeader
        (
            basicThermo::dictName,
            isRegion
          ? runTime.caseConstant()/regionName
          : runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        IOobject materialsHeader
        (
            basicThermo::matDictName,
            isRegion
          ? runTime.caseSystem()/regionName
          : runTime.caseSystem(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        IOobject transportHeader
        (
            "transportProperties",
            isRegion
          ? runTime.caseConstant()/regionName
          : runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if
        (
            !transportHeader.typeHeaderOk<IOdictionary>(true)
         && (
                thermoHeader.typeHeaderOk<IOdictionary>(true)
             || materialsHeader.typeHeaderOk<IOdictionary>(true)
            )
        )
        {
            //thermo exists, but transportProperties doesnt
            //so almost definitely compressible
            compressible = true;
        }
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        forAll(regionMeshes, regioni)
        {
            regionMeshes[regioni].readUpdate();
        }

        #include "createFields.H"

        // here we are init wallHeatFlux
        // which will be printed in a Qconv file later
        volScalarField wallHeatFlux
        (
            IOobject
            (
                "Qconv",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallHeatFlux", dimless, 0.0)
        );

        // Write the total heat-flux including the radiative contribution
        // if available
        autoPtr<volScalarField> totalWallHeatFlux;
        if (Qr.typeHeaderOk<volScalarField>(true))
        {
            totalWallHeatFlux.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "totalWallHeatFlux",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "totalWallHeatFlux",
                        wallHeatFlux.dimensions(),
                        0.0
                    )
                )
            );
        }

        //HTC init
        volScalarField wallHTCconv
        (
            IOobject
            (
                "HTCconv",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallHTCconv", dimless, 0.0)
        );

        // Write the HTC including the radiative contribution
        // if available
        autoPtr<volScalarField> totalWallHTC;
        if (Qr.typeHeaderOk<volScalarField>(true))
        {
            totalWallHTC.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "totalWallHTC",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "totalWallHTC",
                        wallHTCconv.dimensions(),
                        0.0
                    )
                )
            );
        }

        volScalarField wallNusselt
        (
            IOobject
            (
                "wallNusselt",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallNusselt", dimless, 0.0)
        );

        autoPtr<volScalarField> totalwallNusselt;
        if (Qr.typeHeaderOk<volScalarField>(true))
        {
            totalwallNusselt.set
            (
                new volScalarField
                (
                    IOobject
                    (
                       "totalwallNusselt",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar
                    ("totalwallNusselt",wallNusselt.dimensions(), 0.0)
                )
            );
        }

        if (compressible)
        {
            calcCompressibleQ
            (
                mesh,
                runTime,
                U,
                Qr,
                wallHeatFlux,
                totalWallHeatFlux,
                wallHTCconv,
                totalWallHTC,
                printNu,
                wallNusselt,
                totalwallNusselt,
                isRegion,
                regionName
            );
        }
        else
        {
            calcIncompressibleQ
            (
                mesh,
                runTime,
                U,
                Qr,
                wallHeatFlux,
                totalWallHeatFlux,
                wallHTCconv,
                totalWallHTC,
                printNu,
                wallNusselt,
                totalwallNusselt,
                isRegion,
                regionName
            );
        }

        // print wallHeatFlux and HTC
        wallHeatFlux.write();
        wallHTCconv.write();
        if (totalWallHeatFlux.valid())
        {
            totalWallHeatFlux->write();
            totalWallHTC->write();
        }

        // print Nusselt
        if (printNu)
        {
            wallNusselt.write();
            if (totalwallNusselt.valid())
            {
                //print this version including rad
                totalwallNusselt->write();
            }
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
