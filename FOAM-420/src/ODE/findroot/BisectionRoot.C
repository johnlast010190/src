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
    (c) held by original author

Class
    BisectionRoot

Description
    Bisection root Based on Numerical Recipes in C++, Section 9.1, page 358.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "findroot/BisectionRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Func>
const Foam::label Foam::BisectionRoot<Func>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::BisectionRoot<Func>::BisectionRoot(const Func& f, const scalar eps)
:
    f_(f),
    eps_(eps)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Func>
Foam::scalar Foam::BisectionRoot<Func>::root
(
    const scalar x0,
    const scalar x1
) const
{
    scalar f, fMid, dx, rtb, xMid;

    f = f_(x0);
    fMid = f_(x1);

    if (f*fMid >= 0)
    {
        FatalErrorInFunction
            << "Root is not bracketed.  f(x0) = " << f << " f(x1) = " << fMid
            << abort(FatalError);
    }

    // Orient the search such that f > 0 lies at x + dx
    if (f < 0)
    {
        dx = x1 - x0;
        rtb = x0;
    }
    else
    {
        dx = x0 - x1;
        rtb = x1;
    }

    for (label nIter = 0; nIter < maxIter; nIter++)
    {
        dx *= 0.5;
        xMid = rtb + dx;

        fMid = f_(xMid);

        if (fMid <= 0)
        {
            rtb = xMid;
        }

        if (mag(dx) < eps_ || mag(fMid) < SMALL)
        {
            return rtb;
        }
    }

    FatalErrorInFunction
        << "Maximum number of iterations exceeded" << abort(FatalError);

    // Dummy return to keep compiler happy
    return x0;
}


// ************************************************************************* //
