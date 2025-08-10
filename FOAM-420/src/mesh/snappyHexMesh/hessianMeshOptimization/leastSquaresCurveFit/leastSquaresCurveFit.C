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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/leastSquaresCurveFit/leastSquaresCurveFit.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::leastSquaresCurveFit::calculateParametricValues()
{
    scalar totalDistance = 0;
    for (int i=0; i< points_.size()-1; i++)
    {
        totalDistance += mag(points_[i]-points_[i+1]);
    }

    parametricValues_[0] = 0;

    scalar localDistance = 0;
    for (int i=1; i<points_.size(); i++)
    {
        localDistance += mag(points_[i]-points_[i-1]);
        parametricValues_[i] = localDistance/totalDistance;
    }

}

void Foam::leastSquaresCurveFit::solveLeastSquaresSystem()
{
    //Construct A matrix
    scalar A_11=0, A_12=0, A_22 =0 ;
    forAll(points_, pI)
    {
        const scalar& ti = parametricValues_[pI];
        A_11 += 9*ti*ti*pow((1-ti), 4);
        A_12 += 9*ti*ti*ti*pow((1-ti), 3);
        A_22 += 9*ti*ti*ti*ti*(1-ti)*(1-ti);
    }
    //Construct B matrix
    vector B_1 = vector::zero, B_2 = vector::zero;
    label last = points_.size()-1;
    forAll(points_, pI)
    {
        const scalar& ti = parametricValues_[pI];
        B_1 += 3*(points_[pI]- pow((1-ti), 3)*points_[0] - ti*ti*ti*points_[last])
                *(1-ti)*(1-ti)*ti;
        B_2 += 3*(points_[pI]- pow((1-ti), 3)*points_[0] - ti*ti*ti*points_[last])
                        *(1-ti)*ti*ti;
    }

    cP_[0] = points_[0];
    cP_[3] = points_[last];

    //Inverse of A
    scalar invA_11=0, invA_12=0, invA_22 =0 ;
    scalar det = A_11*A_22-A_12*A_12;

    invA_11 = A_22/det;
    invA_22 = A_11/det;
    invA_12 = -A_12/det;

    cP_[1] = invA_11*B_1 + invA_12*B_2;
    cP_[2] = invA_12*B_1 + invA_22*B_2;
}

void Foam::leastSquaresCurveFit::solveQuadraticLeastSquaresSystem()
{
    label last = points_.size()-1;
    cP_[0] = points_[0];
    cP_[2] = points_[last];
    vector coef1 = 0;
    scalar coef2 = 0;

    forAll(points_, pI)
    {
        const scalar& ti = parametricValues_[pI];
        coef1 += points_[pI]-(1-ti)*(1-ti)*cP_[0]-ti*ti*cP_[2];
        coef2 += 2*(1-ti)*ti;
    }

    cP_[1] = coef1/coef2;
}

void Foam::leastSquaresCurveFit::calculateNewPoints()
{
    for (int i=1; i<points_.size()-1; i++)
    {
        const scalar& ti = parametricValues_[i];
        newPoints_[i] = pow((1-ti), 3)*cP_[0] + 3*(1-ti)*(1-ti)*ti*cP_[1] +
                3*(1-ti)*ti*ti*cP_[2] + ti*ti*ti*cP_[3];
    }
}

void Foam::leastSquaresCurveFit::calculateNewPointsQuadratic()
{
    for (int i=1; i<points_.size()-1; i++)
    {
        const scalar& ti = parametricValues_[i];
        newPoints_[i] = sqr(1-ti)*cP_[0] + 2*(1-ti)*ti*cP_[1] + ti*ti*cP_[2];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::leastSquaresCurveFit::leastSquaresCurveFit
(
    const pointField& points
)
:
    points_(points),
    parametricValues_(points.size(), 0),
    newPoints_(points),
    cP_(min(points.size()-1,4), vector::zero)
{
    if (points_.size()<4)
    {
        cP_.resize(0);
    }
    else
//    else if (points_.size()==4)
    {
        calculateParametricValues();
        solveQuadraticLeastSquaresSystem();
        calculateNewPointsQuadratic();
    }
//    else
//    {
//        calculateParametricValues();
//        solveLeastSquaresSystem();
//        calculateNewPoints();
//    }
}

Foam::leastSquaresCurveFit::~leastSquaresCurveFit()
{}

// * * * * * * * * * * * * *  Public Member Functions* * * * * * * * * * * //

Foam::pointField Foam::leastSquaresCurveFit::getNewPoints()
{
    return newPoints_;
}
