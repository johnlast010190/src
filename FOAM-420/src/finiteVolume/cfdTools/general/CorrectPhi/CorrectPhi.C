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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2018-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvScalarMatrix/fvScalarMatrix.H"
#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchFields.H"
#include "cfdTools/general/adjustPhi/adjustPhi.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RAUfType, class DivUType>
void Foam::CorrectPhi
(
    volVectorField& U,
    surfaceScalarField& phi,
    const volScalarField& p,
    const RAUfType& rAUf,
    const DivUType& divU,
    solutionControl& ximple,
    const word& solver
)
{
    const fvMesh& mesh = U.mesh();
    const Time& runTime = mesh.time();

    correctUphiBCs(U, phi);

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (isA<fixedValueZoneFvPatchField<vector>>(U.boundaryField()[patchi]))
        {
            pcorrTypes[patchi] = fixedValueZoneFvPatchVectorField::typeName;
        }
        else if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    if (pcorr.needReference())
    {
        fvc::makeRelative(phi, U);
        adjustPhi(phi, U, pcorr);
        fvc::makeAbsolute(phi, U);
    }

    mesh.schemes().setFluxRequired(pcorr.name());

    while (ximple.correctNonOrthogonal())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi) - divU
        );

        pcorrEqn.setReference(0, 0);

        pcorrEqn.solve
        (
            mesh.solution().solver
            (
                solver == word::null
              ? pcorr.select(ximple.finalNonOrthogonalIter())
              : solver
            )
        );

        if (ximple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

    // Early update of old-time face velocities, for consistency with
    // Rhie-Chow interpolation when phi is corrected.
    // Applied if indirect patches are present.
    bool hasIndirectPatches = false;
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isA<indirectPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            hasIndirectPatches = true;
            break;
        }
    }
    if (hasIndirectPatches)
    {
        const word UfName("Uf");
        if (mesh.foundObject<surfaceVectorField>(UfName))
        {
            surfaceVectorField& Uf =
                mesh.lookupObjectRef<surfaceVectorField>(UfName);
            Uf.oldTime() = fvc::interpolate(U.oldTime());
            surfaceVectorField n(mesh.Sf()/mesh.magSf());
            Uf.oldTime() += n*(phi/mesh.magSf() - (n & Uf.oldTime()));
        }
        else if (mesh.changing())
        {
            WarningInFunction
                << "Could not find surfaceVectorField " << UfName
                << ": Inconsistency in Rhie-Chow" << endl;
        }
    }
}


template<class RAUfType, class DivRhoUType>
void Foam::CorrectPhi
(
    volVectorField& U,
    surfaceScalarField& phi,
    const volScalarField& p,
    const volScalarField& rho,
    const volScalarField& psi,
    const RAUfType& rAUf,
    const DivRhoUType& divRhoU,
    solutionControl& ximple,
    const word& solver
)
{
    const fvMesh& mesh = U.mesh();
    const Time& runTime = mesh.time();

    correctUphiBCs(rho, U, phi);

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    mesh.schemes().setFluxRequired(pcorr.name());

    while (ximple.correctNonOrthogonal())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rAUf, pcorr)
         ==
            divRhoU
        );

        pcorrEqn.solve
        (
            mesh.solution().solver
            (
                solver == word::null
              ? pcorr.select(ximple.finalNonOrthogonalIter())
              : solver
            )
        );

        if (ximple.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux();
        }
    }

    //- consistency with Rhie-Chow if phi is corrected
    {
        word rhoUfName("rhoUf");
        if (mesh.foundObject<surfaceVectorField>(rhoUfName))
        {
            surfaceVectorField& rhoUf = const_cast<surfaceVectorField&>
            (
                mesh.lookupObject<surfaceVectorField>(rhoUfName)
            );
            rhoUf.oldTime() = fvc::interpolate(rho*U.oldTime());
            surfaceVectorField n(mesh.Sf()/mesh.magSf());
            rhoUf.oldTime() += n*(phi/mesh.magSf() - (n & rhoUf.oldTime()));
        }
        else if (mesh.changing())
        {
            WarningInFunction
                << "Could not find surfaceVectorField " << rhoUfName
                << ": Inconsistency in Rhie-Chow" << endl;
        }
    }
}


// ************************************************************************* //
