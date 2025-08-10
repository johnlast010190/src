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

Application
    foamVOF

Group
    grpMultiphaseSolvers

Description
    Solver for multi-phase flows with two or more fluids. This file contains
    implementation details of the member functions of the Cavitate class

\*---------------------------------------------------------------------------*/

#include "cavitate.h"
using namespace Foam;

Cavitate::Cavitate
(
    const MatProp &prop,
    const FoamVOF &vof
)
:
VOFModels
(
    prop,
    vof,
    "cavitation"
)
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    nseed_=scalar(transportProperties.subDict("cavitation").lookupOrDefault<scalar>("n", 1.0e+12));
    R0_=scalar(transportProperties.subDict("cavitation").lookupOrDefault<scalar>("R0", 1.0e-6));
    vapName_=word(transportProperties.subDict("cavitation").lookupOrDefault<word>("vapour_phase_name", "vapour"));
    liqName_=word(transportProperties.subDict("cavitation").lookupOrDefault<word>("liquid_phase_name", "water"));
    rhol_=scalar(transportProperties.subDict("cavitation").lookupOrDefault<scalar>("rhol", 1000));
    rhov_=scalar(transportProperties.subDict("cavitation").lookupOrDefault<scalar>("rhov", 0.1));

    scalar fv=4.0*M_PI*nseed_*pow(R0_,3)/3.0;
    fvmin_=fv/(1.0+fv);
    fvmax_=0.99;


}


const volScalarField& Cavitate::vapour()
{
    if
    (
      vof().alpha1().name()==vapName_
      || vof().alpha1().name()=="alpha."+vapName_
    )
    {
        return  vof().alpha1();
    }
    else
    {
        return vof().alpha2();
    }
}


void Cavitate::specSource
(
    const word& scname,
    volScalarField &Su,
    volScalarField &Sp
)
{
   volTransRate(vof().volTrRate());
   Su=vof().volTrRate();
   return;
}


void Cavitate::volTransRate
(
    volScalarField &vtr
)
{

    vtr.primitiveFieldRef()=0.0;
    scalar cutoff=0.98;
    scalar beta=fvconst();
    scalar coef=2.0/(3.0*rhol_);
    scalar smallP=0.01;
    forAll(vtr.primitiveField(),i)
    {
        scalar alphav(vapour().primitiveField()[i]);
        alphav=max(alphav,fvmin_);
        alphav=min(alphav,fvmax_);
        scalar alphav1(alphav/(beta*(1-alphav)));
        scalar R(pow(alphav1,1.0/3.0));
        scalar dp(prop().psat()-vof().p().primitiveField()[i]);
        scalar dRdt=0.0;
        if (dp>0.0)
        {
          dRdt=sqrt(coef*dp);
        }
        else if (alphav>1.01*fvmin_)
        {
          dRdt=-sqrt(coef*fabs(dp));
        }

        scalar vtri=3.0*alphav*(1.0-alphav)*dRdt/R;
        scalar alpnew=vof().alpha1().oldTime().primitiveField()[i]+vtri*vof().dt();
        if (alpnew>cutoff && vtri>0)
        {
          vtri=(cutoff-vof().alpha1().primitiveField()[i])/vof().dt();
        }
        else if (alpnew<fvmin_)
        {
          vtri=0.0;
        }
        vtr.primitiveFieldRef()[i]=vtri;
        vof().vdotP().primitiveFieldRef()[i]=fabs(vtri)/(fabs(dp)+smallP);
    }

}
