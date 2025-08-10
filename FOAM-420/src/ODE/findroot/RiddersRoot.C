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
    RiddersRoot

Description
    Ridder's method of root finiding given a function, bracketed root
    and accuracy.  Based on Numerical Recipes in C++, Section 9.2,
    page 362.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "findroot/RiddersRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Func>
const Foam::label Foam::RiddersRoot<Func>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::RiddersRoot<Func>::RiddersRoot(const Func& f, const scalar eps)
:
    f_(f),
    eps_(eps)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Func>
Foam::scalar Foam::RiddersRoot<Func>::root
(
    const scalar x0,
    const scalar x1
) const
{
    scalar fl = f_(x0);
    scalar fh = f_(x1);

    // Check bracketing of the root
    if ((fl > 0 && fh < 0) || (fl < 0 && fh > 0))
    {
        scalar xl = x0;
        scalar xh = x1;

        // Bad guess for answer, to simplify logic below
        scalar ans = -1.11e-30;

        scalar xm, fm, s, xNew, fNew;

        for (label nIter = 0; nIter < maxIter; nIter++)
        {
            // Create a mean value and evaluate function
            xm = 0.5*(xl + xh);
            fm = f_(xm);

            // Solve quadratic equation if well-posed
            s = sqrt(sqr(fm) - fl*fh);

            if (s < SMALL)
            {
                return ans;
            }

            // Updating formula
            if (fl >= fh)
            {
                xNew = xm + (xm - xl)*fm/s;
            }
            else
            {
                xNew = xm - (xm - xl)*fm/s;
            }

            // Check answer
            if (mag(xNew - ans) <= eps_)
            {
                return ans;
            }

            ans = xNew;

            fNew = f_(ans);

            if (mag(fNew) < SMALL)
            {
                return ans;
            }

            if (checkSign(fm, fNew))
            {
                xl = xm;
                fl = fm;
                xh = ans;
                fh = fNew;
            }
            else if (checkSign(fl, fNew))
            {
                xh = ans;
                fh = fNew;
            }
            else if (checkSign(fh, fNew))
            {
                xl = ans;
                fl = fNew;
            }
            else
            {
                // Never get here
                FatalErrorInFunction
                    << "Error in search logic" << abort(FatalError);
            }

            if (mag(xh - xl) <= eps_)
            {
                return ans;
            }
        }

        FatalErrorInFunction
            << "Maximum number of iterations exceeded" << abort(FatalError);
    }
    else if (mag(fl) < SMALL)
    {
        return x0;
    }
    else if (mag(fh) < SMALL)
    {
        return x1;
    }
    else
    {
        FatalErrorInFunction
            << "Root is not bracketed.  f(x0) = " << fl << " f(x1) = " << fh
            << abort(FatalError);

    }

    // Dummy return to keep compiler happy
    return x0;
}


// ************************************************************************* //
