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

Description
    Vector-matrix multiplication operations for a block matrix

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::BlockLduMatrix<Type>::decoupledH(const Field<Type>& x) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;

    // Create result
    tmp<Field<Type>> tresult
    (
        new Field<Type>(lduAddr().size(), pTraits<Type>::zero)
    );
    Field<Type>& result = tresult.ref();

    const labelUList& u = lduAddr().upperAddr();
    const labelUList& l = lduAddr().lowerAddr();

    const TypeCoeffField& Upper = this->upper();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Lower multiplication

    if (symmetric())
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[u[coeffI]] -= mult(activeUpper[coeffI], x[l[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeUpper = Upper.asLinear();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[u[coeffI]] -= mult(activeUpper[coeffI], x[l[coeffI]]);
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& Lower = this->lower();

        if (Lower.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeLower = Lower.asScalar();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[u[coeffI]] -= mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
        else if (Lower.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeLower = Lower.asLinear();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[u[coeffI]] -= mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
    }


    // Upper multiplication

    if (Upper.activeType() == blockCoeffBase::SCALAR)
    {
        const scalarTypeField& activeUpper = Upper.asScalar();

        for (label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            result[l[coeffI]] -= mult(activeUpper[coeffI], x[u[coeffI]]);
        }
    }
    else if (Upper.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& activeUpper = Upper.asLinear();

        for (label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            result[l[coeffI]] -= mult(activeUpper[coeffI], x[u[coeffI]]);
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::BlockLduMatrix<Type>::decoupledFaceH(const Field<Type>& x) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;

    const labelUList& u = lduAddr().upperAddr();
    const labelUList& l = lduAddr().lowerAddr();

    // Create result
    tmp<Field<Type>> tresult(new Field<Type>(u.size(), pTraits<Type>::zero));
    Field<Type>& result = tresult.ref();

    const TypeCoeffField& Upper = this->upper();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Lower multiplication

    if (symmetric())
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                // This can be optimised with a subtraction.
                // Currently not done for clarity.  HJ, 31/Oct/2007
                result[coeffI] =
                    mult(activeUpper[coeffI], x[u[coeffI]])
                  - mult(activeUpper[coeffI], x[l[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeUpper = Upper.asLinear();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                // This can be optimised with a subtraction.
                // Currently not done for clarity.  HJ, 31/Oct/2007
                result[coeffI] =
                    mult(activeUpper[coeffI], x[u[coeffI]])
                  - mult(activeUpper[coeffI], x[l[coeffI]]);
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& Lower = this->lower();

        if (Lower.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeUpper = Upper.asScalar();
            const scalarTypeField& activeLower = Lower.asScalar();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[coeffI] =
                    mult(activeUpper[coeffI], x[u[coeffI]])
                  - mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
        else if (Lower.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeUpper = Upper.asLinear();
            const linearTypeField& activeLower = Lower.asLinear();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[coeffI] =
                    mult(activeUpper[coeffI], x[u[coeffI]])
                  - mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
    }

    return tresult;
}


// ************************************************************************* //
