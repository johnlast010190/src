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
    (c) 2011 OpenFOAM Foundation
    (c) 2019-2020 Esi Ltd.

Class
    Foam::MHRICLimiter

Description
    Class with limiter function which returns the limiter for the
    MHRIC differencing scheme based on phict obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedICScheme.

SourceFiles
    MHRIC.C

\*---------------------------------------------------------------------------*/

#ifndef MHRIC_H
#define MHRIC_H

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HyperCLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class MHRICLimiter
:
    public LimiterFunc
{

public:

    MHRICLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const scalar& phiP,
        const scalar& phiN,
        const vector& gradcP,
        const vector& gradcN,
        const vector& d,
        const scalar& Cof,
        const volScalarField& phi,
        const point& p,
        const label  faceI,
        const scalar& blending
    ) const
    {
        scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        scalar MHRIClimiter = max(min(2.0*r, 2.0), 0.0);

        if (Cof <= 0.3)
        {
            MHRIClimiter = max(min(2.0*blending*r, 2.0), 0.0);
        }
        else if (Cof > 0.3 && Cof < 0.7)
        {
            MHRIClimiter = max(min(2.0*blending*r*(0.7-Cof)/(0.7-0.3), 2.0), 0.0);
        }
        else
        {
            MHRIClimiter = 0;
        }

        return MHRIClimiter;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
