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
    (c) 2016 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/solvers/smoothSolver/smoothSolver.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smoothSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<smoothSolver>
        addsmoothSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<smoothSolver>
        addsmoothSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothSolver::smoothSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothSolver::readControls()
{
    lduMatrix::solver::readControls();
    nSweeps_ = controlDict_.lookupOrDefault<label>("nSweeps", 1);
}


Foam::solverPerformance Foam::smoothSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // If the nSweeps_ is negative do a fixed number of sweeps
    if (nSweeps_ < 0)
    {
        #if !defined( WIN32 ) && !defined( WIN64 )
        addProfiling(solve, "lduMatrix::smoother." + fieldName_);
        #endif

        autoPtr<lduMatrix::smoother> smootherPtr = lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        );

        smootherPtr->smooth
        (
            psi,
            source,
            cmpt,
            -nSweeps_
        );

        solverPerf.nIterations() -= nSweeps_;
    }
    else
    {
        scalar normFactor = 0;

        {
            scalarField Apsi(psi.size());
            scalarField temp(psi.size());

            // Calculate A.psi
            matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

            // Calculate normalisation factor
            normFactor = this->normFactor(psi, source, Apsi, temp);

            scalarField residual(source - Apsi);

            matrix().setResidualField(residual, fieldName_);

            // Calculate residual magnitude
            solverPerf.initialResidual() = gSumMag
            (
                residual,
                matrix().mesh().comm()
            )/normFactor;
            solverPerf.finalResidual() = solverPerf.initialResidual();
        }

        if (lduMatrix::debug >= 2)
        {
            Info.masterStream(matrix().mesh().comm())
                << "   " << fieldName_ << "Eqn normalisation factor = "
                << normFactor << endl;
        }


        // Check convergence, solve if not converged or minIter > 0
        if
        (
            minIter_ > 0
         || !solverPerf.checkConvergence(tolerance_, relTol_)
        )
        {
            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(solve, "lduMatrix::smoother." + fieldName_);
            #endif

            autoPtr<lduMatrix::smoother> smootherPtr = lduMatrix::smoother::New
            (
                fieldName_,
                matrix_,
                interfaceBouCoeffs_,
                interfaceIntCoeffs_,
                interfaces_,
                controlDict_
            );

            // Smoothing loop
            do
            {
                smootherPtr->smooth
                (
                    psi,
                    source,
                    cmpt,
                    nSweeps_
                );

                // Calculate the residual to check convergence
                solverPerf.finalResidual() = gSumMag
                (
                    matrix_.residual
                    (
                        psi,
                        source,
                        interfaceBouCoeffs_,
                        interfaces_,
                        cmpt
                    )(),
                    matrix().mesh().comm()
                )/normFactor;
            } while
            (
                (
                    (solverPerf.nIterations() += nSweeps_) < maxIter_
                && !solverPerf.checkConvergence(tolerance_, relTol_)
                )
             || solverPerf.nIterations() < minIter_
            );
        }
    }

    return solverPerf;
}


// ************************************************************************* //
