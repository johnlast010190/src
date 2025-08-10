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

Description
    Vector-matrix multiplication operations for a block matrix

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::BlockLduMatrix<Type>::H(const Field<Type>& x) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    // Create result
    tmp<Field<Type>> tresult
    (
        new Field<Type>(lduAddr().size(), pTraits<Type>::zero)
    );
    Field<Type>& result = tresult.ref();

    if (this->thereIsUpper())
    {
        const labelUList& l = lduAddr().lowerAddr();
        const labelUList& u = lduAddr().upperAddr();
        const TypeCoeffField& Upper = this->upper();

        // Create multiplication function object
        typename BlockCoeff<Type>::multiply mult;

        // Lower multiplication

        if (symmetric())
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& cU = Upper.asScalar();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(cU[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& cU = Upper.asLinear();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(cU[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& cU = Upper.asSquare();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // Use transpose upper coefficient
                    result[u[coeffI]] -=
                        mult(cU[coeffI].T(), x[l[coeffI]]);
                }
            }
        }
        else // Asymmetric matrix
        {
            const TypeCoeffField& Lower = this->lower();

            if (Lower.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& cL = Lower.asScalar();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(cL[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& cL = Lower.asLinear();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(cL[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& cL = Lower.asSquare();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[u[coeffI]] -= mult(cL[coeffI], x[l[coeffI]]);
                }
            }
        }


        // Upper multiplication

        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& cU = Upper.asScalar();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[l[coeffI]] -= mult(cU[coeffI], x[u[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& cU = Upper.asLinear();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[l[coeffI]] -= mult(cU[coeffI], x[u[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& cU = Upper.asSquare();

            for (label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                result[l[coeffI]] -= mult(cU[coeffI], x[u[coeffI]]);
            }
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::BlockLduMatrix<Type>::faceH(const Field<Type>& x) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const labelUList& u = lduAddr().upperAddr();
    const labelUList& l = lduAddr().lowerAddr();

    // Create result
    tmp<Field<Type>> tresult(new Field<Type>(u.size(), pTraits<Type>::zero));
    Field<Type>& result = tresult.ref();

    if (this->thereIsUpper())
    {
        const TypeCoeffField& Upper = this->upper();

        // Create multiplication function object
        typename BlockCoeff<Type>::multiply mult;

        // Lower multiplication

        if (symmetric())
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& cU = Upper.asScalar();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // This can be optimised with a subtraction.
                    // Currently not done for clarity.  HJ, 31/Oct/2007
                    result[coeffI] =
                        mult(cU[coeffI], x[u[coeffI]])
                    - mult(cU[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& cU = Upper.asLinear();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // This can be optimised with a subtraction.
                    // Currently not done for clarity.  HJ, 31/Oct/2007
                    result[coeffI] =
                        mult(cU[coeffI], x[u[coeffI]])
                    - mult(cU[coeffI], x[l[coeffI]]);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& cU = Upper.asSquare();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // Use transpose upper coefficient
                    result[coeffI] =
                        mult(cU[coeffI], x[u[coeffI]])
                    - mult(cU[coeffI].T(), x[l[coeffI]]);
                }
            }
        }
        else // Asymmetric matrix
        {
            const TypeCoeffField& Lower = this->lower();

            if (Lower.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& cU = Upper.asScalar();
                const scalarTypeField& cL = Lower.asScalar();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[coeffI] =
                        mult(cU[coeffI], x[u[coeffI]])
                    - mult(cL[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& cU = Upper.asLinear();
                const linearTypeField& cL = Lower.asLinear();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[coeffI] =
                        mult(cU[coeffI], x[u[coeffI]])
                    - mult(cL[coeffI], x[l[coeffI]]);
                }
            }
            else if (Lower.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& cU = Upper.asSquare();
                const squareTypeField& cL = Lower.asSquare();

                for (label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    result[coeffI] =
                        mult(cU[coeffI], x[u[coeffI]])
                    - mult(cL[coeffI], x[l[coeffI]]);
                }
            }
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::BlockLduMatrix<Type>::faceH
(
    const direction sI,
    const Field<scalar>& x,
    const Field<Type>& psi
) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const labelUList& u = lduAddr().upperAddr();
    const labelUList& l = lduAddr().lowerAddr();

    tmp<Field<scalar>> tresult(new Field<scalar>(u.size(), 0));
    Field<scalar>& result = tresult.ref();

    const direction nCmpts = pTraits<Type>::nComponents;

    if (this->thereIsUpper())
    {
        const TypeCoeffField& Upper = this->upper();

        if (symmetric())
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& cU = Upper.asScalar();

                for (label fI = 0; fI < u.size(); fI++)
                {
                    const scalar& xU = psi[u[fI]][sI];
                    const scalar& xL = psi[u[fI]][sI];
                    result[fI] = cU[fI]*(xU-xL);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& cU = Upper.asLinear();

                for (label fI = 0; fI < u.size(); fI++)
                {
                    const scalar& xU = psi[u[fI]][sI];
                    const scalar& xL = psi[u[fI]][sI];

                    result[fI] = cU[fI][sI]*(xU-xL);
                }
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& cU = Upper.asSquare();

                for (label fI = 0; fI < u.size(); fI++)
                {
                    const scalar& xU = psi[u[fI]][sI];
                    const scalar& xL = psi[u[fI]][sI];
                    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                    {
                        result[fI] +=
                            cU[fI](sI, cmptI)*xU
                          - cU[fI](cmptI, sI)*xL;
                    }
                }
            }
        }
        else
        {
            const TypeCoeffField& Lower = this->lower();

            if (Lower.activeType() == blockCoeffBase::SCALAR)
            {
                const scalarTypeField& cU = Upper.asScalar();
                const scalarTypeField& cL = Lower.asScalar();

                for (label fI = 0; fI < u.size(); fI++)
                {
                    const scalar& xU = psi[u[fI]][sI];
                    const scalar& xL = psi[l[fI]][sI];

                    result[fI] = cU[fI]*xU - cL[fI]*xL;
                }
            }
            else if (Lower.activeType() == blockCoeffBase::LINEAR)
            {
                const linearTypeField& cU = Upper.asLinear();
                const linearTypeField& cL = Lower.asLinear();

                for (label fI = 0; fI < u.size(); fI++)
                {
                    const scalar& xU = psi[u[fI]][sI];
                    const scalar& xL = psi[l[fI]][sI];
                    result[fI] = cU[fI][sI]*xU - cL[fI][sI]*xL;
                }
            }
            else if (Lower.activeType() == blockCoeffBase::SQUARE)
            {
                const squareTypeField& cU = Upper.asSquare();
                const squareTypeField& cL = Lower.asSquare();

                for (label fI = 0; fI < u.size(); fI++)
                {
                    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
                    {
                        const scalar& xU = psi[u[fI]][cmptI];
                        const scalar& xL = psi[l[fI]][cmptI];
                            result[fI] +=
                                cU[fI](sI, cmptI)*xU
                              - cL[fI](sI, cmptI)*xL;
                    }
                }
            }
        }
    }

    return tresult;
}


// ************************************************************************* //
