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
    (c) 2020 Esi Ltd.

Description
    Virtual base class for block iterative solvers

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduSolvers/BlockIterativeSolver/BlockIterativeSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
template<class Type>
Foam::BlockIterativeSolver<Type>::BlockIterativeSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<Type>(fieldName, matrix, dict),
    tolerance_(dict.lookupOrDefault<Foam::scalar>("tolerance", 1.e-30)),
    relTolerance_(dict.lookupOrDefault<Foam::scalar>("relTol", 1.e-1)),
    convCriteria_(dict.lookupOrDefault<Foam::word>
    (
        "convergenceCriteria", "magnitude")
    ),
    component_(0),
    minIter_(dict.lookupOrDefault<Foam::label>("minIter", 3)),
    maxIter_(dict.lookupOrDefault<Foam::label>("maxIter",30)),
    maxFirstIter_(dict.lookupOrDefault<Foam::label>("maxFirstIter", maxIter_))
{
    if (convCriteria_ == "component")
    {
        //- Currently for Up systems, check p equation (continuity)
        component_ = dict.lookupOrDefault<Foam::label>("component", 3);
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::BlockIterativeSolver<Type>::normFactor
(
    Field<Type>& x,
    const Field<Type>& b,
    bool onMaster
) const
{
    const BlockLduMatrix<Type>& matrix = this->matrix_;

    // Calculate the normalisation factor
    const label nRows = x.size();

    Field<Type> pA(nRows);
    Field<Type> wA(nRows);

    // Calculate reference value of x
    Type xRef = pTraits<Type>::zero;
    if (onMaster)
    {
        xRef = average(x);
    }
    else
    {
        xRef = gAverage(x);
    }

    // Calculate A.x
    matrix.Amul(wA, x);

    // Calculate A.xRef, temporarily using pA for storage
    matrix.Amul
    (
        pA,
        Field<Type>(nRows, xRef)
    );

    scalar normFactor = 0;
    if (onMaster)
    {
        normFactor = sum(mag(wA - pA) + mag(b - pA)) + this->small_;
    }
    else
    {
        normFactor = gSum(mag(wA - pA) + mag(b - pA)) + this->small_;
    }

    if (BlockLduMatrix<Type>::debug >= 2)
    {
        Info<< "Iterative solver normalisation factor = " << normFactor << endl;
    }

    return normFactor;
}


template<class Type>
Type Foam::BlockIterativeSolver<Type>::typeNormFactor
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    const BlockLduMatrix<Type>& matrix = this->matrix_;

    // Calculate the normalisation factor
    const label nRows = x.size();

    Field<Type> pA(nRows);
    Field<Type> wA(nRows);

    // Calculate reference value of x
    Type xRef = gAverage(x);

    // Calculate A.x
    matrix.Amul(wA, x);

    // Calculate A.xRef, temporarily using pA for storage
    matrix.Amul
    (
        pA,
        Field<Type>(nRows, xRef)
    );

    Type smallType = pTraits<Type>::one*this->small_;
    Type normFactor = gSum(cmptMag(wA - pA) + cmptMag(b - pA)) + smallType;

    if (BlockLduMatrix<Type>::debug >= 2)
    {
        Info<< "Iterative solver normalisation factor = " << normFactor << endl;
    }

    return normFactor;
}


template<class Type>
bool Foam::BlockIterativeSolver<Type>::stop
(
    BlockSolverPerformance<Type>& solverPerf
) const
{
    if (solverPerf.nIterations() < minIter_)
    {
        return false;
    }
    if (this->timeIndex() == 1)
    {
        if (solverPerf.nIterations() >= maxFirstIter_)
        {
            return true;
        }
    }
    else
    {
        if
        (
            solverPerf.nIterations() >= maxIter_
        )
        {
            return true;
        }
    }

    if (convCriteria_ == "magnitude")
    {
        return solverPerf.checkConvergence(tolerance_, relTolerance_);
    }
    else if (convCriteria_ == "component")
    {
        return solverPerf.checkConvergence
            (
                tolerance_, relTolerance_, component_
            );
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
