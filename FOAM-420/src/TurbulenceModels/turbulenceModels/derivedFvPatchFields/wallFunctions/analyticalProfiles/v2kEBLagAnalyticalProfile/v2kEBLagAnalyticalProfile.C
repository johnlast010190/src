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

\*---------------------------------------------------------------------------*/

#include "v2kEBLagAnalyticalProfile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::v2kEBLagAnalyticalProfile::v2kEBLagAnalyticalProfile()
{
    v2k_.setSize(5,3);

    v2k_[0][0]  = 1;  v2k_[0][1]  = 0.00147;  v2k_[0][2]  = 0.002299;
    v2k_[1][0]  = 3;  v2k_[1][1]  = 0.008;    v2k_[1][2]  = 0.0123;
    v2k_[2][0]  = 5;  v2k_[2][1]  = 0.0165;   v2k_[2][2]  = 0.0264;
    v2k_[3][0]  = 11; v2k_[3][1]  = 0.0596;   v2k_[3][2]  = 0.088;
    v2k_[4][0]  = 30; v2k_[4][1]  = 0.177;    v2k_[4][2]  = 0.138;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::v2kEBLagAnalyticalProfile::getv2k(const scalar& yPlus)
{
    scalar v2k = 0;

    if (yPlus < 1)//v2kPlus[0][0])
    {
        v2k = 0.0015*yPlus*yPlus;
    }
    else if (yPlus >= v2k_[4][0])
    {
        v2k = min(2./3., 0.3077*log(yPlus)/log(10.0) - 0.2775);
    }
    else
    {
        int i = 0;
        while (i != 4)
        {
        if (yPlus >= v2k_[i][0] && yPlus < v2k_[i+1][0])
            {
            scalar dy = log(v2k_[i+1][0]/v2k_[i][0]);
            scalar t = 1./dy*log(yPlus/v2k_[i][0]);
            v2k = (2*t*t*t - 3*t*t + 1)*v2k_[i][1] +
                  (3*t*t - 2*t*t*t)*v2k_[i+1][1] +
                  (t*t*t - 2*t*t + t)*dy*v2k_[i][2] +
                  (t*t*t - t*t)*dy*v2k_[i+1][2];
            break;
            }
        ++i;
        }
    }

    return v2k;
}
