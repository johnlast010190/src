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
    Finite-Area matrix basic solvers.

\*---------------------------------------------------------------------------*/

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"
#include "fields/Fields/diagTensorField/diagTensorField.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<class Type>
void faMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    internalCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

    boundaryCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]*value;
}

template<class Type>
Foam::solverPerformance Foam::faMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "faMatrix<Type>::solve(const dictionary& solverControls) : "
               "solving faMatrix<Type>"
            << endl;
    }

    label maxIter = -1;
    if (solverControls.readIfPresent("maxIter", maxIter))
    {
        if (maxIter == 0)
        {
            return solverPerformance();
        }
    }

    word type(solverControls.lookupOrDefault<word>("type", "segregated"));

    if (type == "segregated")
    {
        return solveSegregated(solverControls);
    }
    else
    {
        FatalIOErrorIn
        (
            "faMatrix<Type>::solve(const dictionary& solverControls)",
            solverControls
        )   << "Unknown type " << type
            << "; currently supported solver types are segregated"// and coupled"
            << exit(FatalIOError);

        return solverPerformance();
    }
}

template<class Type>
Foam::solverPerformance Foam::faMatrix<Type>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "faMatrix<Type>::solve(const dictionary&) : "
               "solving faMatrix<Type>"
            << endl;
    }

    GeometricField<Type, faPatchField, areaMesh>& psi =
       const_cast<GeometricField<Type, faPatchField, areaMesh>&>(psi_);

    solverPerformance solverPerfVec
    (
        "faMatrix<Type>::solve",
        psi_.name()
    );

    this->boundaryManipulate(psi.boundaryFieldRef());

    scalarField saveDiag(diag());

    Field<Type> source(source_);
    addBoundarySource(source, psi_);

    typename Type::labelType validComponents
    (
        pow
        (
            psi.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        // copy field and source

        scalarField psiCmpt( psi_.primitiveField().component(cmpt) );
        addBoundaryDiag(diag(), solvingComponent);

        scalarField sourceCmpt(source.component(cmpt));

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        lduInterfaceFieldPtrsList interfaces( psi.boundaryField().scalarInterfaces() );

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        initMatrixInterfaces
        (
            true,
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        updateMatrixInterfaces
        (
            true,
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        if (solverPerformance::debug)
        {
            solverPerf.print(Info);
        }

        solverPerfVec = max(solverPerfVec, solverPerf);
        solverPerfVec.solverName() = solverPerf.solverName();

        psi.primitiveFieldRef().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerfVec);

    return solverPerfVec;
}

template<class Type>
autoPtr<typename faMatrix<Type>::faSolver> faMatrix<Type>::solver()
{
    return solver
    (
        psi_.mesh().solution().solverDict
        (
            psi_.select
            (
                psi_.mesh().data::template lookupOrDefault<bool>
                ("finalIteration", false)
            )
        )
    );
}

template<class Type>
solverPerformance faMatrix<Type>::faSolver::solve()
{
     return solve
     (
         faMat_.psi_.mesh().solution().solverDict
         (
             faMat_.psi_.select
             (
                 faMat_.psi_.mesh().data::template lookupOrDefault<bool>
                 ("finalIteration", false)
             )
         )
     );
}


template<class Type>
solverPerformance faMatrix<Type>::solve()
{
     return solve
     (
         psi_.mesh().solution().solverDict
         (
             psi_.select
             (
                 psi_.mesh().data::template lookupOrDefault<bool>
                 ("finalIteration", false)
             )
         )
     );
}


// Return the matrix residual
template<class Type>
tmp<Field<Type>> faMatrix<Type>::residual() const
{
    tmp<Field<Type>> tres(new Field<Type>(source_));
    Field<Type>& res = tres.ref();

    addBoundarySource(res, psi_);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
//    lduInterfaceFieldPtrsList interfaces = psi_.boundaryFieldRef().interfaces();

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.primitiveField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                psi_.boundaryField().scalarInterfaces(),
                cmpt
            )
        );
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
