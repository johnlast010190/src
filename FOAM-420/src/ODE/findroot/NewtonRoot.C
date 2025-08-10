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
    NewtonRoot

Description
    Newton root Based on Numerical Recipes in C++, Section 9.1, page 358.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)


\*----------------------------------------------------------------------------*/

#include "findroot/NewtonRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Func>
const Foam::label Foam::NewtonRoot<Func>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::NewtonRoot<Func>::NewtonRoot(const Func& f, const scalar eps)
:
    f_(f),
    eps_(eps)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Func>
Foam::scalar Foam::NewtonRoot<Func>::root
(
    const scalar x0

) const
{
    scalar f, df, resid=1, Iter=1,xstart=x0;

    while (Iter<maxIter)
    {
     f = f_(xstart);
     df = f_.d(xstart);
     resid=f_(xstart);
     xstart=xstart-f/df;
     if (mag(resid)<eps_)
        { return xstart;}

     Iter=Iter+1;

    }
    FatalErrorInFunction
        << "Maximum number of iterations exceeded" << abort(FatalError);

    // Dummy return to keep compiler happy
    return xstart;
}




// ************************************************************************* //
