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

#include "FkEBLagAnalyticalProfile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::FkEBLagAnalyticalProfile::FkEBLagAnalyticalProfile(const bool& lag)
:
lag_(lag)
{
    if (lag_)
    {
        Fk_.setSize(9,3);
        Fk_[0][0] = 0;       Fk_[0][1] = 0.0;     Fk_[0][2] = 10;
//        Fk_[1][0] = 0.5; Fk_[1][1] = 1.53;    Fk_[1][2] = 0.77;
//        Fk_[2][0] = 5;   Fk_[2][1] = 2.846;   Fk_[2][2] = 0.143;
        Fk_[1][0] = 0.5;     Fk_[1][1] = 1.2;     Fk_[1][2] = 0.53;
        Fk_[2][0] = 5;       Fk_[2][1] = 2.337;   Fk_[2][2] = 0.123;
        Fk_[3][0] = 40;      Fk_[3][1] = 4.343;   Fk_[3][2] = 0.0316;
        Fk_[4][0] = 300;     Fk_[4][1] = 8.561;   Fk_[4][2] = 0.0123;
        Fk_[5][0] = 1000;    Fk_[5][1] = 14.12;   Fk_[5][2] = 0.006295;
        Fk_[6][0] = 10000;   Fk_[6][1] = 48.25;   Fk_[6][2] = 0.002325;
        Fk_[7][0] = 100000;  Fk_[7][1] = 169.451; Fk_[7][2] = 0.000859;
        Fk_[8][0] = 1000000; Fk_[8][1] = 584.451; Fk_[8][2] = 0.000359;
    }
    else
    {
        Fk_.setSize(8,3);
        Fk_[0][0] = 0;   Fk_[0][1] = 0.0;       Fk_[0][2] = 10;
        Fk_[1][0] = 0.5; Fk_[1][1] = 1.53;    Fk_[1][2] = 0.77;
        Fk_[2][0] = 5;   Fk_[2][1] = 2.846;   Fk_[2][2] = 0.143;
        Fk_[3][0] = 40;  Fk_[3][1] = 5.058;   Fk_[3][2] = 0.0386;
        Fk_[4][0] = 300; Fk_[4][1] = 9.926;   Fk_[4][2] = 0.0123;
        Fk_[5][0] = 1000; Fk_[5][1] = 16.21;   Fk_[5][2] = 0.007295;
        Fk_[6][0] = 10000; Fk_[6][1] = 53.64;   Fk_[6][2] = 0.003025;
        Fk_[7][0] = 1000000; Fk_[7][1] = 217.451; Fk_[7][2] = 0.00159;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::FkEBLagAnalyticalProfile::getYPlus(const scalar& x)
{
    label i=-1;
    for (int j=0; j<Fk_.m()-1; ++j)
    {
        if (x>=Fk_[j][0] && x<=Fk_[j+1][0])
        {
            i = j;
            break;
        }
    }
    scalar Dx = Fk_[i+1][0]-Fk_[i][0];
    scalar t = (x-Fk_[i][0])/Dx;
    scalar yPlus =
        (2*t*t*t -3*t*t + 1)*Fk_[i][1]
      + (t*t*t -2*t*t + t)*Fk_[i][2]*Dx
      + (3*t*t - 2*t*t*t)*Fk_[i+1][1]
      + (t*t*t - t*t)*Fk_[i+1][2]*Dx;

    if (i == 0 && lag_)
    {
        yPlus = pow(t, 0.2) + pow(t, 1.5);// - 0.44;
    }

    return yPlus;
}

Foam::scalar Foam::FkEBLagAnalyticalProfile::maxFkValue()
{
    return 1e06;
}
