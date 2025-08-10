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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description
     Finite-Area scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "faMatrices/faScalarMatrix/faScalarMatrix.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<>
void faMatrix<scalar>::setComponentReference
(
    const label patchI,
    const label edgeI,
    const direction,
    const scalar value
)
{
    const labelUList& faceLabels =
        psi_.mesh().boundary()[patchI].edgeFaces();

    internalCoeffs_[patchI][edgeI] +=
        diag()[faceLabels[edgeI]];

    boundaryCoeffs_[patchI][edgeI] = value;
}


template<>
autoPtr<faMatrix<scalar>::faSolver> faMatrix<scalar>::solver
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "faMatrix<scalar>::solver(const dictionary&) : "
               "solver for faMatrix<scalar>"
            << endl;
    }

    scalarField saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    autoPtr<faMatrix<scalar>::faSolver> solverPtr
    (
        new faMatrix<scalar>::faSolver
        (
            *this,
            lduMatrix::solver::New
            (
                psi_.name(),
                *this,
                boundaryCoeffs_,
                internalCoeffs_,
                psi_.boundaryField().scalarInterfaces(),
                solverControls
            )
        )
    );

    diag() = saveDiag;

    return solverPtr;
}


template<>
solverPerformance faMatrix<scalar>::faSolver::solve
(
    const dictionary& solverControls
)
{
    GeometricField<scalar, faPatchField, areaMesh>& psi =
        const_cast<GeometricField<scalar, faPatchField, areaMesh>&>
        (faMat_.psi());

    scalarField saveDiag = faMat_.diag();
    faMat_.addBoundaryDiag(faMat_.diag(), 0);

    scalarField totalSource = faMat_.source();
    faMat_.addBoundarySource(totalSource, faMat_.psi(), false);

    // assign new solver controls
    solver_->read(solverControls);

    solverPerformance solverPerf =
        solver_->solve(psi.primitiveFieldRef(), totalSource);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info);
    }

    faMat_.diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
solverPerformance faMatrix<scalar>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "faMatrix<scalar>::solve(const dictionary&) : "
               "solving faMatrix<scalar>"
            << endl;
    }
    GeometricField<scalar, faPatchField, areaMesh>& psi =
        const_cast<GeometricField<scalar, faPatchField, areaMesh>&>
        (psi_);

    this->boundaryManipulate(psi.boundaryFieldRef());

    scalarField saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    scalarField totalSource = source_;
    addBoundarySource(totalSource, psi_, false);

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psi_.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        psi.boundaryField().scalarInterfaces(),
        solverControls
    )->solve(psi.primitiveFieldRef(), totalSource);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info);
    }

    diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


// Return the matrix residual
template<>
tmp<scalarField> faMatrix<scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.primitiveField(),
            source_ - boundaryDiag*psi_.primitiveField(),
            boundaryCoeffs_,
            psi_.boundaryField().scalarInterfaces(),
            0
        )
    );

    addBoundarySource(tres.ref(), psi_);

    return tres;
}


// H operator
template<>
tmp<areaScalarField> faMatrix<scalar>::H() const
{
    tmp<areaScalarField> tHphi
    (
        new areaScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField Hphi = tHphi();

    Hphi.primitiveFieldRef() = (lduMatrix::H(psi_.primitiveField()) + source_);
    addBoundarySource(Hphi.primitiveFieldRef(), psi_);

    Hphi.primitiveFieldRef() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
