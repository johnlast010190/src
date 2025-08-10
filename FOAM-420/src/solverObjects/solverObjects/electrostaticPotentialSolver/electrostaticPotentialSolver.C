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
    (c) 2019-2021 Esi Ltd.
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "electrostaticPotentialSolver.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(electrostaticPotentialSolver, 0);
}
}

makeFvSolverOption(electrostaticPotentialSolver);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::coordinateSystem&
Foam::fv::electrostaticPotentialSolver::coordSys() const
{
    if (!coordinates_.valid() && !coorFramePtr_)
    {
        FatalErrorInFunction
            << "Co-ordinate system invalid"
            << abort(FatalError);
    }
    else if (coorFramePtr_)
    {
        return (*coorFramePtr_).coorSys();
    }

    return coordinates_();
}


void Foam::fv::electrostaticPotentialSolver::transformSigma
(
    const volVectorField& sigmaLocal,
    volSymmTensorField& aniSigma
)
{
    aniSigma_().primitiveFieldRef() =
        coordSys().transformPrincipal(mesh_.C(), aniLocalSigma_());
    forAll(aniSigma_().boundaryField(), patchi)
    {
        fvPatchSymmTensorField& pf = aniSigma_().boundaryFieldRef()[patchi];
        if (!isA<emptyFvPatch>(pf.patch()))
        {
            pf.forceAssign
            (
                coordSys().transformPrincipal
                (
                    mesh_.C().boundaryField()[patchi],
                    aniLocalSigma_().boundaryField()[patchi]
                )
            );
        }
    }
}


void Foam::fv::electrostaticPotentialSolver::updateSigma()
{
    if (anisotropic_)
    {
        updateSigma(aniSigmaFn_, aniLocalSigma_());
        transformSigma(aniLocalSigma_(), aniSigma_());
    }
    else
    {
        updateSigma(sigmaFn_, sigma_());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::electrostaticPotentialSolver::electrostaticPotentialSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    solnControlPtr_(nullptr),
    V_
    (
        volScalarField
        (
            IOobject
            (
                "electrical_V",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    anisotropic_(false),
    coorFramePtr_(nullptr)
{
    // TODO: This should probably be moved to thermophysical properties as a
    // new electricalProperties section
    dict.lookup("anisotropicConductivity") >> anisotropic_;
    if (anisotropic_)
    {
        aniSigmaFn_ = Function1<vector>::New("sigma", dict);
    }
    else
    {
        sigmaFn_ = Function1<scalar>::New("sigma", dict);
    }

    if (anisotropic_)
    {
        Info<< "    Using anisotropic electrical conductivity" << endl;
        aniLocalSigma_.set(createSigma(aniSigmaFn_).ptr());
        aniLocalSigma_().rename("electical_localSigma");
        if (dict.found("referenceFrame"))
        {
            coorFramePtr_ = coordinateFrame::lookupNew(mesh_, dict);
        }
        else
        {
            WarningInFunction
                << "coordinateSystem will be deprecated. "
                << "Please use the referenceFrame instead." << nl << endl;
            coordinates_ = coordinateSystem::New(mesh_, dict);
        }
        aniSigma_.set
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "electrical_sigma",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedSymmTensor("0", aniLocalSigma_().dimensions(), Zero)
            )
        );
    }
    else
    {
        Info<< "    Using scalar electrical conductivity" << endl;
        sigma_.set(createSigma(sigmaFn_).ptr());
    }

    // Only update the sigma field in initialise() (needs T)

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::electrostaticPotentialSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("electrical_V");

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
    correctorMembers.insert
    (
        "nonOrthogonalCorrector:"+solveNames[0], solveNames
    );
}


bool Foam::fv::electrostaticPotentialSolver::initialise()
{
    // Read controls
    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    updateSigma();
    return true;
}


bool Foam::fv::electrostaticPotentialSolver::isFinalCorrector
(
    const label corrector,
    const word& correctorName
)
{
    if (correctorName == "nonOrthogonalCorrector:"+V_.name())
    {
        // Non-orthogonal correctors
        return (corrector >= solnControlPtr_->nNonOrthCorr());
    }
    else if (correctorName == solverObject::outerCorrectorName)
    {
        return true;
    }
    else
    {
        return false;
    }
}


tmp<fvScalarMatrix>
Foam::fv::electrostaticPotentialSolver::assembleScalarMatrix
(
    const word& fieldName,
    bool& finalSolve,
    word& dictName
)
{
    fv::options& fvOptions = this->fvOptions();

    updateSigma();

    autoPtr<volVectorField> gradV;
    // Pre-store the grad for use in non-orthog correction and BCs
    // Prevents duplicate calculation for anisotropic regions; we call it
    // regardless to allow consolidated group to contain both isotropic and
    // anisotropic regions
    gradV.set(fvc::grad(V_).ptr());

    // Assemble electrical potential equation
    tmp<fvScalarMatrix> tVEqn
    (
        new fvScalarMatrix
        (
            anisotropic_ ?
           -fvm::laplacian(aniSigma_(), V_) :
           -fvm::laplacian(sigma_(), V_)
        )
    );

    tVEqn->relax();
    fvOptions.constrain(tVEqn.ref());

    return tVEqn;
}


void Foam::fv::electrostaticPotentialSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    this->fvOptions().correct();
}


// ************************************************************************* //
