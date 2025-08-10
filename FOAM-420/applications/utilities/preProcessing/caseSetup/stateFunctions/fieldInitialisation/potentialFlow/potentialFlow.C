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
    (c) 2016-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "potentialFlow.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "stateFunction/stateFunction.H"
#include "cfdTools/general/include/fvCFD.H"
#include "cfdTools/general/adjustU/adjustU.H"
#include "sources/derived/MRFSource/IOMRFSourceList.H"
#include "cfdTools/general/adjustSplitBCs/adjustSplitBCs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(potentialFlow, 0);
    addToRunTimeSelectionTable(fieldInit, potentialFlow, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

bool Foam::fieldInitializations::potentialFlow::solvePotentialFlow
(
    const dictionary& potDict
) const
{
    word fieldName = fieldInit::name();

    if (word(potDict.lookup("type")) != word("potentialFlow"))
    {
        //- Not potential flow initialization
        return false;
    }

    if (localDb().foundObject<volVectorField>(word(fieldName+"pot")))
    {
        //- If U*pot is there, no need to solve this field since it is already
        //  done
        return true;
    }
    else if (localDb().foundObject<volScalarField>(word(fieldName+"Pot")))
    {
        //- If pressure p*Pot is there, no need to solve this field since it
        // is already done
        return true;
    }
    else if
    (
        fieldName == "p_rgh" &&
        localDb().foundObject<volScalarField>(word("pPot"))
    )
    {
        //exception
        return true;
    }


    //potential flow variables
    autoPtr<volVectorField> UPtr;
    if (localDb().foundObject<volVectorField>(name()))
    {
        UPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    word(fieldName+"pot"),
                    mesh().time().timeName(),
                    localDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                localDb().lookupObject<volVectorField>(name())
            )
        );
    }
    else if (localDb().foundObject<volVectorField>("U"))
    {
        UPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    word("Upot"),
                    mesh().time().timeName(),
                    localDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                localDb().lookupObject<volVectorField>("U")
            )
        );
    }
    else
    {
        FatalError << "No velocity field defined, required for "
            << "potential flow initialisation." << exit(FatalError);
    }

    UPtr->primitiveFieldRef() = vector::zero;

    autoPtr<volScalarField> pPtr;

    if (localDb().foundObject<volScalarField>("p_rgh"))
    {
        pPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pPot",
                    mesh().time().timeName(),
                    localDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                localDb().lookupObject<volScalarField>("p_rgh")
            )
        );
    }
    else if (localDb().foundObject<volScalarField>("p"))
    {
        pPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pPot",
                    mesh().time().timeName(),
                    localDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                localDb().lookupObject<volScalarField>("p")
            )
        );
    }
    else
    {
        FatalError << "No pressure field defined, required for "
            << "potential flow initialisation." << exit(FatalError);
    }

    //- Some of the states create a rho (rhoCalculated) as a field with values
    //  of 1.0 kg/m^3. Here we need to create a local rho with values
    //  of either read from disc or set using rhoRef.
    //  This rho is not registered to the registry.
    //  If rho exists on memory then the one of memory is synced with the local
    //  one because BCs will look the one in the registry.
    //  This is not the best for now; Requires caseSetup refactoring to be
    //  cleaner

    dimensionedScalar rhoRef
    (
        "rhoRef",
        dimensionSet(1, -3, 0, 0, 0),
        1.0
    );
    bool pdiv = false;

    const word rhoName(potDict.lookupOrDefault<word>("rhoName", "rho"));

    bool rhoInRegistry = localDb().foundObject<volScalarField>(rhoName);

    IOobject rhoIO
    (
        rhoName,
        mesh().time().timeName(),
        localDb(),
        mesh().conformal() ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
        IOobject::NO_WRITE,
        !rhoInRegistry
    );

    if (pPtr->dimensions() == dimensionSet(1, -1, -2, 0, 0, 0, 0))
    {
        if (potDict.found("rhoRef"))
        {
            rhoRef.value() = readScalar(potDict.lookup("rhoRef"));
            Info<< "rhoRef is found and set to " << rhoRef.value() << endl;
        }
        else
        {
            if (rhoIO.headerOk())
            {
                Info<< "Density field " << rhoIO.name() << " found on disc."
                     << " Selecting it for initialiation."
                     << endl;
            }
            else
            {
                rhoRef.value() = 1.205;
                Info<< "rhoRef is not found. Setting it to air (default): "
                     << rhoRef.value() << endl;
            }
        }
        pdiv = true;
    }

    autoPtr<volScalarField> initRho;

    if (pdiv)
    {
        initRho.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh(),
                rhoRef,
                zeroGradientFvPatchScalarField::typeName
            )
        );
        initRho->correctBoundaryConditions();
        if (rhoInRegistry)
        {
            //- sync local rho with the one in registry if exists
            const volScalarField& crhoMat =
                localDb().lookupObject<volScalarField>(rhoName);
            volScalarField& rhoMat = const_cast<volScalarField&>(crhoMat);
            rhoMat = initRho();
        }
    }

    // special treatment for pTot-p systems
    bool pTotPSystem(potDict.lookupOrDefault<bool>("pTotPSystem", false));

    if (pTotPSystem)
    {
        if (fieldName == "p")
        {
            FatalErrorInFunction
                << "Unsupported initialisation potentialFlow "
                << "for p in ptot-p systems "
                << exit(FatalError);
        }

        // optional input parameter
        word inletPatchName
        (
            potDict.lookupOrDefault<word>("inletPatchName", "inlet")
        );
        word outletPatchName
        (
            potDict.lookupOrDefault<word>("outletPatchName", "outlet")
        );
        scalar deltaP(potDict.lookupOrDefault<scalar>("deltaP", 0.0));

        // Initialize BCs list from p
        wordList pcorrTypes
        (
            pPtr().boundaryField().types()
        );

        // average inlet/outlet pressure
        scalar pInlet = 0;
        scalar pOutlet = 0;

        // Set inlet BC for p to zeroGradient and
        // compute inlet/outlet pressure
        forAll(pPtr().boundaryField(), patchi)
        {
            if (pPtr().boundaryField()[patchi].patch().name() == inletPatchName)
            {
                pcorrTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
                pInlet = gAverage(pPtr().boundaryField()[patchi]);
            }
            if (pPtr().boundaryField()[patchi].patch().name() == outletPatchName)
            {
                pOutlet = gAverage(pPtr().boundaryField()[patchi]);
            }
        }

        // compute and reset U inlet
        scalar Uin = sqrt(2*(pInlet-pOutlet-deltaP));

        forAll(UPtr().boundaryField(), patchi)
        {
            if (UPtr().boundaryField()[patchi].patch().name() == inletPatchName)
            {
                UPtr().boundaryFieldRef()[patchi].forceAssign(
                    Uin*mesh().Sf().boundaryField()[patchi]
                  / mesh().magSf().boundaryField()[patchi]
                );
            }
        }

        if (debug)
        {
            Info<< "pInlet = " << pInlet << endl;
            Info<< "pOutlet = " << pOutlet << endl;
            Info<< "deltaP = " << deltaP << endl;
            Info<< "U inlet = " << Uin << endl;
        }

        // reset p with new bc types
        pPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pPot",
                    mesh().time().timeName(),
                    localDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                localDb().lookupObject<volScalarField>("p"),
                pcorrTypes
            )
        );
    }

    //generate unit 1/AU field for boundary condition compatibility
    volScalarField rUA
    (
        IOobject
        (
            "(1|A(U))",
            mesh().time().timeName(),
            localDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("(1|A(U))", dimTime, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );
    if (pdiv)
    {
        rUA.ref() /= initRho();
        rUA.correctBoundaryConditions();
    }


    fvMesh& refMesh = const_cast<fvMesh&>(mesh());
    dictionary fvSchemes = refMesh.schemes().localSchemeDict();
    {
        if (!fvSchemes.found("gradSchemes"))
        {
            fvSchemes.add
            (
                word("gradSchemes"), dictionary(), false
            );
        }
        char pPotGrad[] = "Gauss linear";

        fvSchemes.subDict("gradSchemes").add
        (
            word("grad(pPot)"), pPotGrad, false
        );

        fvSchemes.subDict("gradSchemes").add
        (
            word("snGradCorr(pPot)"), pPotGrad, false
        );


        if (!fvSchemes.found("laplacianSchemes"))
        {
            fvSchemes.add(word("laplacianSchemes"), dictionary(), false);
        }

        char pPotLapl[]
            = "Gauss linear limited 0.333";

        fvSchemes.subDict("laplacianSchemes").add
        (
            word("laplacian(1,pPot)"), pPotLapl, false
        );
        fvSchemes.subDict("laplacianSchemes").add
        (
            word("laplacian((1|A(U)),pPot)"), pPotLapl, false
        );

        if (!fvSchemes.found("fluxRequired"))
        {
            fvSchemes.add(word("fluxRequired"), dictionary(), false);
        }
        fvSchemes.subDict("fluxRequired").add(word("pPot"), "");


        if (!fvSchemes.found("interpolationSchemes"))
        {
            fvSchemes.add(word("interpolationSchemes"), dictionary(), false);
        }
        fvSchemes.subDict("interpolationSchemes").add
        (
            word("interpolate(Upot)"),
            "linear",
            false
        );

        refMesh.schemes().setLocalSchemeDict(fvSchemes);
    }

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh().time().timeName(),
            localDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(UPtr(), "interpolate(Upot)") & mesh().Sf()
    );

    dictionary refDict;
    word refPointName(pPtr().name()+"RefPoint");
    word refCellName(pPtr().name()+"RefCell");
    word refValueName(pPtr().name()+"RefValue");
    if (!potDict.found(refPointName))
    {
        refDict.add
        (
            refCellName,
            potDict.lookupOrDefault<label>(refCellName, 0)
        );
    }
    else
    {
        refDict.add
        (
            refPointName, potDict.lookupOrDefault<vector>
            (
                refPointName,
                vector::zero
            )
        );
    }
    refDict.add
    (
        refValueName,
        potDict.lookupOrDefault<scalar>(refValueName, 0.0)
    );
    label pRefCell = 0;
    scalar pRefValue = 0.0;

    setRefCell
    (
        pPtr(),
        pPtr(),
        refDict,
        pRefCell,
        pRefValue
    );

    label nNonOrthCorr
        = potDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 10);

    //run potFlow
    setOutFlowU(rhoIO, phi, UPtr());
    adjustPhi(phi, UPtr(), pPtr());
    fv::IOMRFSourceList MRF(localDb());

    if (!MRF.PtrList::size())
    {
        Info<< "No finite volume options present" << endl;
    }

    Switch initialiseUBCs
        = potDict.lookupOrDefault<Switch>("initialiseUBCs", false);

    label nCorrectors
        = potDict.lookupOrDefault<label>("nCorrectors", initialiseUBCs ? 3 : 1);

    if (initialiseUBCs)
    {
        Info<< "   " << "Initialising velocity boundary conditions" <<endl;
        UPtr->correctBoundaryConditions();

        //set phi again since there is a dependence on U
        phi = (fvc::interpolate(UPtr()) & mesh().Sf());
    }

    for (label corr = 0; corr < nCorrectors; corr++)
    {
        adjustSplitBCs(phi, UPtr(), pPtr());
        bool adjustedU = adjustSplitBCs(phi, UPtr(), UPtr());

        if (adjustedU || (corr && initialiseUBCs))
        {
            // In the case of split velocity BCs, we need an initial
            // correctBoundaryConditions before adjustSplitBCs as well as this
            // one after, in order to get the correct velocity into both inlet
            // and outlet BCs
            Info<< "   " << "Updating velocity boundary conditions" <<endl;
            UPtr->correctBoundaryConditions();

            //set phi again since there is a dependence on U
            phi = (fvc::interpolate(UPtr()) & mesh().Sf());
        }

        Info<< "   " << "Switching the absolute flux to relative" <<endl;
        MRF.makeRelative(phi);
        adjustPhi(phi, UPtr(), pPtr());

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    rUA,
                    pPtr()
                )
            ==
                fvc::div(phi)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve(potDict.subDict("pSolver"));

            if (nonOrth == nNonOrthCorr)
            {
                phi -= pEqn.flux();
            }
        }

        Info<< "   " << "Switching the relative flux to absolute" <<endl;
        MRF.makeAbsolute(phi);

        Info<< "continuity error = "
            << mag(fvc::div(phi))().weightedAverage(mesh().V()).value()
            << endl;

        UPtr() = fvc::reconstruct(phi);
        UPtr->correctBoundaryConditions();

        Info<< "   " << "Interpolated Velocity error = "
            << (sqrt(sum(sqr((fvc::interpolate(UPtr())
                & mesh().Sf()) - phi)))
                /sum(mesh().magSf())).value()
            << endl;

    }

    //store potentialflow variables for later assignment
    UPtr->store(UPtr);
    pPtr->store(pPtr);

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::potentialFlow::potentialFlow
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::potentialFlow::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    dictionary potDict
    (
        stateFunction::etcDictionary
        (
            "caseDicts/preProcessing/caseSetup/potentialFlowInitialisation.cfg",
            true
        )().subDict("potentialFlow")
    );

    potDict.merge(initDict());

    solvePotentialFlow(potDict);

    if (localDb().foundObject<volVectorField>(name()))
    {
        volVectorField& f
            = const_cast<volVectorField&>
            (localDb().lookupObject<volVectorField>(name()));

        f.primitiveFieldRef()
            = localDb().lookupObject<volVectorField>
            (word(name()) + "pot").primitiveField();

        forAll(f.boundaryField(), pI)
        {
            f.boundaryFieldRef()[pI].fvPatchVectorField::operator=
            (
                fvPatchVectorField
                (
                    localDb().lookupObject<volVectorField>(word(name()) + "pot")
                    .boundaryField()[pI]
                )
            );
        }

        scalar sigmaU
            = initDict().lookupOrDefault<scalar>("sigmaU", 10.0);

        scalar totalVolume = gSum(mesh().V());
        scalar meanF = gSum(mag(f.primitiveField()) * mesh().V());
        meanF /= totalVolume;

        scalar stdDevF
            = gSum(mesh().V() * magSqr(mag(f.primitiveField()) - meanF));
        stdDevF /= totalVolume;
        stdDevF = Foam::sqrt(stdDevF);
        scalar maxF = meanF + sigmaU * stdDevF;

        Info<< "   " << "Velocity magnitude limited to: " << maxF << endl;

        forAll(f, i)
        {
            scalar magfi = mag(f[i]);

            if (magfi > maxF)
            {
                f[i] = f[i] * (maxF/magfi);
            }
        }

        volVectorField::Boundary& fbf = f.boundaryFieldRef();
        forAll(fbf, patchi)
        {
            if (!fbf[patchi].fixesValue())
            {
                forAll(fbf[patchi], facei)
                {
                    vector& fi = fbf[patchi][facei];
                    scalar magfi = mag(fi);

                    if (magfi > maxF)
                    {
                        fi = fi * maxF/magfi;
                    }
                }
            }
        }

        //calculate reasonable boundary layer thickness
        //mean wall dist * constant
        volScalarField d = wallDist(mesh()).y();

        scalar cBL = 0.5;

        dimensionedScalar blThickness
        (
            "boundaryLayerThickness",
            dimLength,
            gAverage(d)*cBL
        );

        if (initDict().found("boundaryLayerThickness"))
        {
            blThickness.value() = readScalar
            (
                initDict().lookup("boundaryLayerThickness")
            );
        }

        Info<< "   " << "Assumed boundary layer thickness: "
             << blThickness.value() << " m" << endl;


        // Modify velocity by applying a 1/7th power law boundary-layer
        // u/U0 = (d/D)^(1/7)
        // assuming U0 is the same as the current cell velocity

        forAll(f, cellI)
        {
            if (d[cellI] <= blThickness.value())
            {
                f[cellI] = f[cellI]
                    *::pow(d[cellI]/blThickness.value(), (1.0/7.0));
            }
        }

    }
    else
    {
        volScalarField& f
            = const_cast<volScalarField&>
            (localDb().lookupObject<volScalarField>(name()));

        f.primitiveFieldRef()
            = localDb().lookupObject<volScalarField>("pPot").primitiveField();

        forAll(f.boundaryField(), pI)
        {
            f.boundaryFieldRef()[pI].fvPatchScalarField::operator=
            (
                fvPatchScalarField
                (
                    localDb().lookupObject<volScalarField>("pPot")
                    .boundaryField()[pI]
                )
            );
        }
    }

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
