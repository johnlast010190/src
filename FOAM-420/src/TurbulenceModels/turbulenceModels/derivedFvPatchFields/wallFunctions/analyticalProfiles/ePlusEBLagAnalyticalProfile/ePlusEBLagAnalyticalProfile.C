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

#include "ePlusEBLagAnalyticalProfile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::ePlusEBLagAnalyticalProfile::ePlusEBLagAnalyticalProfile(const bool& lag)
:
    lag_(lag)
{
    ePlus_.setSize(5,3);

    if (lag_)
    {
        ePlus_[0][0]  = 1;  ePlus_[0][1]  = 0.0803; ePlus_[0][2]  = -0.045;
        ePlus_[1][0]  = 3;  ePlus_[1][1]  = 0.0774; ePlus_[1][2]  = 0.045;
        ePlus_[2][0]  = 5;  ePlus_[2][1]  = 0.127;  ePlus_[2][2]  = 0.122;
        ePlus_[3][0]  = 11; ePlus_[3][1]  = 0.206;  ePlus_[3][2]  = -0.07;
        ePlus_[4][0]  = 30; ePlus_[4][1]  = 0.0863; ePlus_[4][2]  = -0.074;
    }
    else
    {
        ePlus_[0][0]  = 1;  ePlus_[0][1]  = 0.0732;  ePlus_[0][2]  = -0.021;
        ePlus_[1][0]  = 3;  ePlus_[1][1]  = 0.063;  ePlus_[1][2]  = 0.025;
        ePlus_[2][0]  = 5;  ePlus_[2][1]  = 0.091;  ePlus_[2][2]  = 0.082;
        ePlus_[3][0]  = 11; ePlus_[3][1]  = 0.149;  ePlus_[3][2]  = 0.0;
        ePlus_[4][0]  = 30; ePlus_[4][1]  = 0.0793; ePlus_[4][2]  = -0.074;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::ePlusEBLagAnalyticalProfile::getePlus(const scalar& yPlus)
{
    scalar ePlus = 0;
    scalar Ceps2 = 1.83;

    if (yPlus < 1)
    {
        ePlus = 5*1600/Ceps2/Ceps2/pow((yPlus+12.43), 4.);
    }
    else if (yPlus >= 30)
    {
        ePlus = 1./0.41/yPlus;
    }
    else
    {
        int i = 0;
        while (i != 4)
        {
            if (yPlus >= ePlus_[i][0] && yPlus < ePlus_[i+1][0])
            {
                scalar Dy = log(ePlus_[i+1][0]/ePlus_[i][0]);
                scalar t = 1/Dy*log(yPlus/ePlus_[i][0]);
                ePlus =
                    (2*t*t*t - 3*t*t + 1)*ePlus_[i][1] +
                    (3*t*t - 2*t*t*t)*ePlus_[i+1][1] +
                    (t*t*t - 2*t*t + t)*Dy*ePlus_[i][2] +
                    (t*t*t - t*t)*Dy*ePlus_[i+1][2];
            break;
            }
            i++;
        }
    }

    return ePlus;
}
