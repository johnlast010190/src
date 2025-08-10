    Info<< "Reading Temperature field \n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
         ),
        mesh
    );

    Info<< "Reading Velocity field \n" << endl;
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

    volScalarField alphaConv
    (
        IOobject
        (
            "wallheattransfercoefficient",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
            ("alphaConv", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
    );

    if (args.optionFound("compressible"))
    {
        IOobject rhoHeader
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

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

        forAll(alphaConv.boundaryField(), i)
        {
            if (!alphaConv.boundaryField()[i].coupled())
            {
                alphaConv.boundaryFieldRef()[i]
                    = thermo.Cp()->boundaryField()[i]
                    * turb->alphaEff()().boundaryField()[i]
                    * mesh.boundary()[i].deltaCoeffs();
            }
        }
    }
    else
    {
        #include "cfdTools/incompressible/createPhi.H"
        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turb
        (
            incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );

        forAll(alphaConv.boundaryField(), i)
        {
            if (!alphaConv.boundaryField()[i].coupled())
            {
                alphaConv.boundaryFieldRef()[i]
                    = turb->Cp()->boundaryField()[i]
                    * turb->alphaEff()().boundaryField()[i]
                    * mesh.boundary()[i].deltaCoeffs();
            }
        }
    }
