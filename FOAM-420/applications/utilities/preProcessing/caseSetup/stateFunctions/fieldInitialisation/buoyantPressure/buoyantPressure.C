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
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "buoyantPressure.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"
#include "rhoThermo/rhoThermo.H"
#include "stateFunction/stateFunction.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "planarAverage/planarAverage.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(buoyantPressure, 0);
    addToRunTimeSelectionTable(fieldInit, buoyantPressure, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::buoyantPressure::buoyantPressure
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::buoyantPressure::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    label nNonOrthCorr =
        initDict().lookupOrDefault<label>("nNonOrthogonalCorrectors", 10);

    // Number of piecewise-linear height segments
    // NOTE: Can default to -1 to have it automatically determined based on min
    // edge length, but this may be overkill
    label nAveragingPlanes =
        initDict().lookupOrDefault<label>("nAveragingPlanes", 1000);

    bool resetT = initDict().lookupOrDefault<bool>("resetT", false);
    // This should be unnecessary due to the planar averaging performed on rho
    if (initDict().found("resetT"))
    {
        DeprecationIOWarningInFunction(initDict(), "resetT", "option", 40100)
            << nl << endl;
    }

    fvMesh& refMesh = const_cast<fvMesh&>(solReg().mesh());

    // Can be p or p_rgh
    volScalarField& p_rgh =
        const_cast<volScalarField&>
        (
            localDb().lookupObject<volScalarField>(name())
        );

    word phaseName = initDict().lookupOrDefault<word>("phaseName", "");

    Info<< "Reading thermophysical properties\n" << endl;
    rhoThermo* thermoPtr =
        &refCast<rhoThermo>
        (
            multiphaseThermo::lookupOrCreate(localDb(), phaseName)
        );
    const bool distinctBuoyancy
    (
        basicThermo::dictName == basicThermo::matDictName
     && thermoPtr->distinctBuoyancy()
    );

    if (initDict().found("initialValue"))
    {
        // Pre-initialise with supplied value in order to bootstrap the density
        // calculation
        p_rgh.primitiveFieldRef() =
            scalarField(word("initialValue"), initDict(), p_rgh.size());
        initCoupledBoundaries();
    }
    else if (!thermoPtr->isochoric())
    {
        IOWarningInFunction(initDict())
            << "No 'initialValue' entry specified for starting iteration of "
            << "buoyant pressure calculation - "
            << p_rgh.name() << " field's default value will be used."
            << endl;
    }

    volScalarField& p = thermoPtr->p();

    const Time& runTime = refMesh.time();
    const fvMesh& mesh = refMesh;
    Info<< "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            localDb(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    Info<< "\nReading hRef" << endl;
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            localDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("hRef", dimLength, 0)
    );
    #include "cfdTools/general/include/gh.H"

    // Non-orthogonal correction for boundary gradient
    forAll(gh.boundaryField(), patchi)
    {
        const vectorField pd(mesh.boundary()[patchi].delta());
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        if
        (
            !isA<emptyFvPatchField<scalar>>(gh.boundaryField()[patchi])
         && !gh.boundaryField()[patchi].coupled()
        )
        {
            fvPatchScalarField& pgh = gh.boundaryFieldRef()[patchi];
            forAll(pgh, bfi)
            {
                pgh[bfi] =
                    (g.value() & (mesh.C()[fc[bfi]] + pd[bfi]))
                  - ghRef.value();
            }
            ghf.boundaryFieldRef()[patchi].forceAssign
            (
                gh.boundaryField()[patchi]
            );
        }
    }

    // Create dummy zero U field to pass to constrainPressure
    volVectorField Uzero
    (
        IOobject
        (
            "Uzero",
            runTime.constant(),
            localDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("Zero", dimVelocity, vector::zero)
    );

    // phi may be needed for BC's
    const volVectorField& U =
        localDb().lookupObject<volVectorField>
        (
            initDict().lookupOrDefault<word>("UName", "U")
        );

    volScalarField& T = thermoPtr->T();
    volScalarField Tcpy(T);

    if (resetT)
    {
        scalar Tave = gSum(T.primitiveField()*mesh.V())/gSum(mesh.V());
        T.forceAssign(dimensionedScalar("T", T.dimensions(), Tave));
        thermoPtr->he() = thermoPtr->he(p, T);
        thermoPtr->correct();
    }

    volScalarField rho
    (
        IOobject
        (
            initDict().lookupOrDefault<word>("rhoName", "rho"),
            refMesh.time().timeName(),
            localDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        distinctBuoyancy ? thermoPtr->buoyantRho()() : thermoPtr->rho()
    );
    rho =
        planarAverage::planarAverage
        (
            rho, g.value()/stabilise(mag(g.value()), SMALL), nAveragingPlanes
        );

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
            (
                "phi",
                runTime.timeName(),
                localDb(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
        fvc::interpolate(rho)*fvc::flux(U)
    );

    //solve steady buoyantPressure equation
    Info<< "Solving the Steady buoyantPressure Equation" << endl;

    bool isp_rgh = (p_rgh.member() == "p_rgh");
    if (isp_rgh)
    {
        p = p_rgh + rho*gh;
    }
    thermoPtr->correct();

    for (label i=0; i<nNonOrthCorr; i++)
    {
        surfaceScalarField phig
        (
            "phig",
           -ghf
           *fvc::snGrad(rho, "snGrad("+p().name()+')', "grad("+p().name()+')')
           *refMesh.magSf()
        );
        if (!isp_rgh)
        {
            phig +=
                fvc::snGrad
                (
                    rho*gh, "snGrad("+p().name()+')', "grad("+p().name()+')'
                )
               *refMesh.magSf();
        }

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, rho, Uzero, phig, geometricOneField());

        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(p_rgh) == fvc::div(phig)
        );

        if (p_rgh.needReference())
        {
            label pRefCell(0);
            scalar pRefValue(0);
            setRefCell
            (
                p_rgh,
                mesh.solution().subDictPtr("PIMPLE")
              ? mesh.solution().subDict("PIMPLE")
              : mesh.solution().subDict("SIMPLE"),
                pRefCell,
                pRefValue
            );

            p_rghEqn.setReference
            (
                pRefCell,
                getRefCellValue(p_rgh, pRefCell)
            );
        }

        if (refMesh.solution().solverDictExists(name()))
        {
            p_rghEqn.solve(refMesh.solution().solver(name()));
        }
        else
        {
            dictionary pSolverDict
            (
                stateFunction::etcDictionary
                (
                    "caseDicts/preProcessing/caseSetup/potentialFlowInitialisation.cfg",
                    true
                )().subDict("potentialFlow")
            );

            p_rghEqn.solve(pSolverDict.subDict("pSolver"));
        }

        if (isp_rgh)
        {
            p = p_rgh + rho*gh;
        }
        // Prevent T from changing due to change in p
        thermoPtr->he().forceAssign(thermoPtr->he(p, Tcpy));
        thermoPtr->correct();
        rho =
            planarAverage::planarAverage
            (
                distinctBuoyancy
              ? thermoPtr->buoyantRho()()
              : thermoPtr->rho(),
                g.value()/stabilise(mag(g.value()), SMALL),
                nAveragingPlanes
            );
    }
    T.forceAssign(Tcpy);

    // the field has been initialised
    initialised() = true;
}



// ************************************************************************* //
