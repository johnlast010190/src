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
    (c) 2004-6 H. Jasak All rights reserved
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/BlockSolverPerformance.H"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance
)
{
    if (BlockLduMatrix<Type>::debug >= 2)
    {
        Info<< solverName_
            << ":  Iteration " << nIterations_
            << " residual = " << finalResidual_
            << endl;
    }

    if
    (
        mag(finalResidual_) < Tolerance
     || (
            RelTolerance > SMALL
         && mag(finalResidual_) <= RelTolerance*mag(initialResidual_)
        )
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance,
    const label comp
)
{
    scalar initResComp = initialResidual_.component(comp);
    scalar finalResComp = finalResidual_.component(comp);

    return checkConvergence(Tolerance, RelTolerance, initResComp, finalResComp);
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance,
    const scalar initResComp,
    const scalar finalResComp
)
{
    if (BlockLduMatrix<Type>::debug >= 2)
    {
        Info<< solverName_
            << ":  Iteration " << nIterations_
            << " residual = " << finalResidual_
            << endl;
    }

    if
    (
        finalResComp < Tolerance
     || (
            RelTolerance > SMALL
         &&
            (
                finalResComp <=
                RelTolerance*initResComp
            )
        )
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkSingularity
(
    const scalar& residual
)
{
    if (mag(residual) > VSMALL)
    {
        singular_ = false;
    }
    else
    {
        singular_ = true;
    }

    return singular_;
}


template<class Type>
void Foam::BlockSolverPerformance<Type>::print() const
{
    /*
    if (BlockLduMatrix<Type>::debug)
    {
        Info<< solverName_ << ":  Solving for " << fieldName_;

        if (singular())
        {
            Info<< ":  solution singularity" << endl;
        }
        else
        {
            Info<< ", Initial residual = " << initialResidual_
                << ", Final residual = " << finalResidual_
                << ", No Iterations " << nIterations_
                << endl;
        }
    }
    */
    Info<< solverName_ << ":  Solving for " << fieldName_;
    Info<< ", Initial residual = " << initialResidual_
        << ", Final residual = " << finalResidual_
        << ", No Iterations " << nIterations_
        << endl;
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::operator!=
(
    const BlockSolverPerformance<Type>& sp
) const
{
    return
    (
        this->solverName()      != sp.solverName()
     || this->fieldName()       != sp.fieldName()
     || this->initialResidual() != sp.initialResidual()
     || this->finalResidual()   != sp.finalResidual()
     || this->nIterations()     != sp.nIterations()
     || this->converged()       != sp.converged()
     || this->singular()        != sp.singular()
    );
}


template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    typename Foam::BlockSolverPerformance<Type>& sp
)
{
    is.readBeginList("SolverPerformance<Type>");
    is  >> sp.solverName_
        >> sp.fieldName_
        >> sp.initialResidual_
        >> sp.finalResidual_
        >> sp.nIterations_
        >> sp.converged_
        >> sp.singular_;
    is.readEndList("SolverPerformance<Type>");

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const typename Foam::BlockSolverPerformance<Type>& sp
)
{
    os  << token::BEGIN_LIST
        << sp.solverName_ << token::SPACE
        << sp.fieldName_ << token::SPACE
        << sp.initialResidual_ << token::SPACE
        << sp.finalResidual_ << token::SPACE
        << sp.nIterations_ << token::SPACE
        << sp.converged_ << token::SPACE
        << sp.singular_ << token::SPACE
        << token::END_LIST;

    return os;
}





// ************************************************************************* //
