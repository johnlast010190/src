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
    (c) 2022 Esi Ltd.

Application
    velocityProfile


\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    //Constants
    scalar argDen = 50;
    scalar argDen2 = 27;//27;
    scalar k = 0.41;
    int N = 3;

    Info<<"yPlus"<<"  " <<"uF"<<nl;
    for (int i=1; i<400000; ++i)
    {
        scalar yPlus = 0.004*i;

        scalar F_Rei = Foam::log(1.+0.4*yPlus)/k +
                7.8*(1-Foam::exp(-yPlus/11.)-yPlus/11.0*Foam::exp(-yPlus/3.0));

        scalar F_log = Foam::log(yPlus)/k + 5.1;

        scalar arg = yPlus/argDen2;
        scalar phi =  Foam::tanh(arg*arg*arg*arg);
        scalar F_Rei_Blend = (1-phi)*F_Rei + phi*F_log;

        bool converged = false;
        scalar uPlus_Old;
        scalar uPlus_New;
        scalar Fsp, dFsp;
        scalar c1 = Foam::exp(-k*5.2);

        uPlus_Old = F_Rei_Blend;
        while (!converged)
        {
            scalar Sum = 1;
            for (int n=1; n<N; ++n)
            {
                Sum += pow(k*uPlus_Old, n)/scalar(factorial(n));
            }
            Fsp  = uPlus_Old + c1*(Foam::exp(k*uPlus_Old)-Sum) - yPlus;

            scalar c2 = Foam::exp(k*uPlus_Old);

            dFsp = 1 + c1*(k*c2 - k - k*uPlus_Old - 0.5*(k*uPlus_Old)*(k*uPlus_Old));

            uPlus_New = uPlus_Old - Fsp/dFsp;

            scalar deltaU = mag(uPlus_New - uPlus_Old);

//            Info<< "Old U "<<uPlus_Old<<endl;
//            Info<< "New U "<<uPlus_New<<endl;
//            Info<<"F is: "<< Fsp<<endl;
//            Info<<"dF is: "<< dFsp<<endl;

            uPlus_Old = uPlus_New;
            if (deltaU < 0.000001)
            {
                converged = true;
            }
        }

        scalar F_Spald = uPlus_New;

        arg = yPlus/argDen;
        phi = Foam::tanh(arg*arg);
        //scalar uF = (1-phi)*F_Spald + phi*F_Rei_Blend;
        scalar uF = F_Spald;
        Info<<yPlus<<" " <<uF<<nl;

        scalar uTau = 0.03917126174;
        scalar Cmu = 0.09;
        scalar kappa = 0.41;
        scalar E = 9.81;
        scalar nuw = 1.51e-05;
        scalar beta1  = 0.075;
        scalar Rey = yPlus/pow025(Cmu);
        scalar lamFrac = Foam::exp(-Rey/11);
        scalar turbFrac = 1 - lamFrac;
        scalar omegaVis = 6*uTau*uTau/beta1/nuw/yPlus/yPlus;

        scalar uPlus = (1/kappa)*Foam::log(E*yPlus);
        scalar omegaLog = uTau*uTau/(Foam::sqrt(Cmu)*kappa*nuw*yPlus);
        scalar omegab1 = omegaVis + omegaLog;

        scalar omegab2 = Foam::pow(Foam::pow(omegaVis,1.2) + Foam::pow(omegaLog,1.2),1/1.2);

        scalar argNew = 0.1*yPlus;

        scalar phiNew = Foam::tanh(pow(argNew,4));

//        scalar omega =  phiNew*omegab1 + (1-phiNew)*omegab2;
 //       scalar omega = lamFrac*omegaVis + turbFrac*omegaLog;
        label n_ = 2;
        scalar omega = Foam::pow
        (
             Foam::pow(omegaVis, n_) + Foam::pow(omegaLog, n_),
             scalar(1)/n_
         );

//        Info<<yPlus<<" " <<omega<<nl;

    }


    return 0;
}

// ************************************************************************* //
