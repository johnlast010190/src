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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    splineBasis

Description
    B-spline basis function

\*---------------------------------------------------------------------------*/

#include "vNurbsDeformation/splineBasis/splineBasis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// private member functions

void splineBasis::createKnots()
{
    for (label i = 0; i <= nDegree; i++)
    {
        knots[i]= 0.0;
    }

    label fI = nDegree + 1;
    label lI  = nKnots - nDegree - 1;
    label size  = nKnots - 2*nDegree - 2;
    for (label i = 0; i < size; i++)
    {
        knots[i + fI] = scalar(i+1)/scalar(size+1);
    }

    for (label i = 0; i < nDegree + 1; i++)
    {
        knots[i + lI] = 1.0;
    }
}

// Constructor

splineBasis::splineBasis
(
    const label nPoints,
    const label deg
)
{
    nCPs = nPoints;

    nDegree = deg;

    nKnots = nPoints+nDegree+1;

    knots = scalarField(nKnots, 0.0);

    this->createKnots();
}

// copy constructor

splineBasis::splineBasis
(
    const splineBasis& basis
) :
    nCPs(basis.nCPs),
    nDegree(basis.nDegree),
    nKnots(basis.nKnots),
    knots(basis.knots)
{}

// public member functions

scalar splineBasis::value
(
    const label& iCP,
    const label& deg,
    const scalar& u
) const
{

    if (deg == 0)
    {
        if ((u >= knots[iCP]) && (u < knots[iCP+1]))
        {
            return 1.0;
        }
        else if ((u == 1.0) && (knots[iCP+1] == 1))
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    else
    {
        scalar denom1 = (knots[iCP+deg] - knots[iCP]);

        scalar denom2 = (knots[iCP + deg+1] - knots[iCP+1]);

        scalar value(0.0);

        if (denom1 != 0.0)
        {
            label d = deg-1;
            value += (u - knots[iCP])*this->value(iCP, d, u) / denom1;
        }
        if (denom2 != 0.0)
        {
            label d = deg-1;
            label i = iCP+1;
            value += (knots[iCP+deg+1] - u)*this->value(i, d, u) / denom2;
        }

        return value;
    }

}

scalar splineBasis::derivative
(
    const label& iCP,
    const label& deg,
    const scalar& u
) const
{

    if (deg == 0)
    {
        return 0.0;
    }
    else
    {
        scalar denom1 = (knots[iCP+deg] - knots[iCP]);

        scalar denom2 = (knots[iCP + deg + 1] - knots[iCP + 1]);

        scalar derivative(0.0);

        if (denom1 != 0.0)
        {
            label d = deg-1;
            derivative += ( (u - knots[iCP])*this->derivative(iCP, d, u) + this->value(iCP, d, u) ) /denom1;
        }
        if (denom2 != 0.0)
        {
            label d = deg-1;
            label i = iCP+1;
            derivative +=( (knots[iCP+deg+1] - u)*this->derivative(i, d, u) - this->value(i, d, u) ) / denom2;
        }
        return derivative;
    }
}

label splineBasis::span
(
    const scalar& u
) const
{
    if (u == knots[nCPs]) return nCPs-1;

    label low, high, mid;

    low = nDegree;
    high = nCPs;
    mid = (high+low)/2;

    while (u < knots[mid] || u >= knots[mid+1])
    {
        if (u < knots[mid]) high = mid;
    else low = mid;

    mid = (low+high)/2;
    }

    return mid;
}

const scalarField& splineBasis::getKnots() const
{
    return knots;
}

void splineBasis::changeKnots
(
    const scalarField& k
)
{
    if (nKnots != k.size())
    {
        FatalErrorInFunction
            << "The length of the knot array is not equal to nCP + degree + 1. "
            << exit(FatalError);
    }

    forAll(knots, i)
    {
        knots[i] = k[i];
    }

    return;
}

const label& splineBasis::nbPoints() const
{
    return nCPs;
}

const label& splineBasis::nbKnots() const
{
    return nKnots;
}

const label& splineBasis::degree() const
{
    return nDegree;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
