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
    (c) 2019 Esi Ltd.

Class
    Foam::v2kLagAnalyticalProfile

Description
    Analytical profile of v2k for Elliptic Blending Lag kEpsilon.

SourceFiles
    v2kLagAnalyticalProfile.C

\*---------------------------------------------------------------------------*/
#ifndef v2kEBLagAnalyticalProfile_H
#define v2kEBLagAnalyticalProfile_H

#include "cfdTools/general/include/fvCFD.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class v2kLagAnalyticalProfile Declaration
\*---------------------------------------------------------------------------*/
class v2kEBLagAnalyticalProfile
{

    //- Hermite Polynomial Weights
    RectangularMatrix <scalar> v2k_;

public:

    //- Constructor

        //- Empty Constructor
        v2kEBLagAnalyticalProfile();

    //- Destructor
    ~v2kEBLagAnalyticalProfile(){}

    // Member Functions

    scalar getv2k(const scalar& yPlus);

};

}//End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
