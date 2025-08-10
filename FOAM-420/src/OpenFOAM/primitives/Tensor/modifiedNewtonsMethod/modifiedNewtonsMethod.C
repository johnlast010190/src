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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "modifiedNewtonsMethod.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::modifiedNewtonsMethod::modifiedCholesky(tensor& tA)
{
    //Modified Cholesky factorization which produces a new matrix A' = A + E
    //Where E is a (small enough) diagonal contribution so that A' is positive
    //definite. If A is already positive definite then E = 0

    //Make tensor as a List of List for easier implementation of the algorithm
    List<List<scalar>> A(3);
    forAll(A, i)
    {
        A[i].setSize(3, 0.);
    }
    for (int i=0; i<9; i++)
    {
        A[i/3][i%3] = tA.v_[i];
    }

    //Create Auxiliary Identity Matrix P
    List<List<scalar>> P(3);
    forAll(P, i)
    {
        P[i].setSize(3, 0.);
        P[i][i] = 1.;
    }

    //Machine Accuracy
    scalar macheps = 1e-19;

    //Max off-diagonal element
    scalar ksi = max(max(mag(A[1][0]),mag(A[2][0])),mag(A[2][1]));

    //Max diagonal element
    scalar gamma = max(max(mag(A[0][0]),mag(A[1][1])),mag(A[2][2]));

    scalar delta = macheps*max( (gamma + ksi), 1);

    //Diagonal Hessian Modification
    vector e=vector::zero;

    //GMW81 Algorithm
    for (int j=0; j<3; j++)
    {
        scalar th = 0;
        for (int i=0; i<3; i++)
        {
            scalar sum = 0;
            for (int k=0; k<i; k++)
            {
                sum += P[k][i]*P[k][j];
            }
            P[i][j] = (A[i][j]-sum)/P[i][i];

            if ((A[i][j] - sum) > th)
            {
                th = A[i][j] -sum;
            }
            if (i>j)
            {
                P[i][j] = 0;
            }
        }
        scalar sum = 0;
        for (int i=0; i<j; i++)
        {
            sum += sqr(P[i][j]);
        }

        scalar phij = A[j][j] - sum;

        scalar di = 0;

        if (j<=1)
        {
            if (j==0)
            {
                di = max(mag(A[1][0]), mag(A[2][0]));
            }
            else
            {
                di = mag(A[2][1]);
            }
        }
        else
        {
            di = mag(A[2][2]);
        }

        scalar b = Foam::sqrt(max(macheps, max(gamma , di/3.)));

        if (delta >= max(mag(phij), sqr(th/b)))
        {
            e[j] = delta - phij;
        }
        else if (mag(phij) >= max(sqr(delta/b),delta))
        {
            e[j] = mag(phij) - phij;
        }
        else if (sqr(th/b) >= max(delta, mag(phij)))
        {
            e[j] = sqr(th/b) - phij;
        }

        P[j][j] = Foam::sqrt(mag(A[j][j] - sum + e[j]));
    }

    //Add diagonal modification
    tA.v_[0] += e[0];
    tA.v_[4] += e[1];
    tA.v_[8] += e[2];
}

Foam::vector Foam::modifiedNewtonsMethod::newtonsStep
(
    const vector& G,
    const tensor& T
)
{
    vector displacement = vector::zero;
    tensor hessian = T;

    vector eigen = eigenValues(hessian);
    if
    (
         eigen.x()<0
     || eigen.y()<0
     || eigen.z()<0
    )
    {
        modifiedCholesky(hessian);
    }

    scalar dett(det(hessian));

    if (mag(dett)>VSMALL)
    {
        hessian = inv(hessian);
    }
    else
    {
        hessian = tensor::zero;
    }

    displacement = hessian&G;

    return displacement;
}
// ************************************************************************* //
