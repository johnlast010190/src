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

    scalar TRef;
    bool useTRef = false;
    if (args.optionReadIfPresent("TRef", TRef))
    {
        useTRef = true;
        Info<< "Reference Temperature used to evaluate HTC: " << TRef <<" K"<<endl;
    }

    if (args.optionFound("compressible"))
    {
        Info<< "Reading thermophysical properties\n" << endl;

        autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
        rhoThermo& thermo = pThermo();
        thermo.validate(args.executable(), "h", "e");

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
            thermo.rho()
        );

        #include "cfdTools/compressible/compressibleCreatePhi.H"

        Info<< "Creating turbulence model\n" << endl;
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

        /*
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
        */

        forAll(alphaConv.boundaryField(), i)
        {
            if (!alphaConv.boundaryField()[i].coupled())
            {
                if (useTRef)
                {
                    alphaConv.boundaryFieldRef()[i] =
                    (
                        thermo.Cp()->boundaryField()[i]
                        * turb->alphaEff()().boundaryField()[i]
                        * fvc::snGrad(T)->boundaryField()[i]
                    )
                    / stabilise(T.boundaryField()[i] - TRef, VSMALL);
                }
                else
                {
                    alphaConv.boundaryFieldRef()[i]
                        = thermo.Cp()->boundaryField()[i]
                        * turb->alphaEff()().boundaryField()[i]
                        * mesh.boundary()[i].deltaCoeffs();
                }
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
                if (useTRef)
                {
                    alphaConv.boundaryFieldRef()[i] =
                    (
                        turb->Cp()->boundaryField()[i]
                        * turb->rho()->boundaryField()[i]
                        * turb->alphaEff()().boundaryField()[i]
                        * fvc::snGrad(T)->boundaryField()[i]
                    )
                    / stabilise(T.boundaryField()[i] - TRef, VSMALL);
                }
                else
                {
                    alphaConv.boundaryFieldRef()[i]
                        = turb->Cp()->boundaryField()[i]
                        * turb->rho()->boundaryField()[i]
                        * turb->alphaEff()().boundaryField()[i]
                        * mesh.boundary()[i].deltaCoeffs();
                }
            }
        }
    }
