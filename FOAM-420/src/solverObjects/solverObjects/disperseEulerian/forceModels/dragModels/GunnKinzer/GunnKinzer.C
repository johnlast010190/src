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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2020-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "GunnKinzer.H"
#include "disperseEulerian/phase/phase.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decoupledEulerian
{
    defineTypeNameAndDebug(GunnKinzer, 0);
    addToRunTimeSelectionTable(dragModel, GunnKinzer, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decoupledEulerian::GunnKinzer::GunnKinzer
(
    const dictionary& dict,
    const phase& phase,
    const bool registerObject
)
:
    dragModel(dict, phase, registerObject),
    autoMaterialProp_(true)
{
    if (dict.found("GunnKinzerCoeffs"))
    {
        autoMaterialProp_ =
            dict.subDict("GunnKinzerCoeffs").lookupOrDefault<Switch>("autoMaterialProperties", true);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decoupledEulerian::GunnKinzer::~GunnKinzer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::decoupledEulerian::GunnKinzer::CdRe() const
{
    volScalarField Re(phase_.Re());

    if (autoMaterialProp_)
    {
       Re = rhoa(phase_.Td())*phase_.magUr()*phase_.diam()/mua(phase_.Td());
    }

    tmp<volScalarField> tCdRe
    (
        new volScalarField
        (
            IOobject
            (
                "tCdRe",
                Re.instance(),
                Re.mesh()
            ),
            Re
        )
    );
    volScalarField& CdRe = tCdRe.ref();

    //--------------------------------------------------------------------------------------
    // Matrices containing Reynolds numbers and the corresponding drag coefficient values
    // according to Gunn & Kinzer
    double Re_M[]={1.80, 9.61, 23.4, 43.2, 68.7, 98.9, 134.0, 175.0, 220.0, 269.0,
              372.0, 483.0, 603.0, 731.0, 866.0, 1013.0, 1164.0, 1313.0, 1461.0, 1613.0,
              1764.0, 1915.0, 2066.0, 2211.0, 2357.0, 2500.0, 2636.0, 2772.0, 2905.0, 3033.0,
              3164.0, 3293.0, 3423.0, 3549.0};
    double Cd_M[]={15.0, 4.2, 2.4, 1.66, 1.28, 1.07, 0.926, 0.815, 0.729, 0.671,
              0.607, 0.570, 0.545, 0.528, 0.517, 0.504, 0.495, 0.494, 0.498, 0.503,
              0.511, 0.520, 0.529, 0.544, 0.559, 0.575, 0.594, 0.615, 0.635, 0.660,
              0.681, 0.700, 0.727, 0.751};
    //--------------------------------------------------------------------------------------

    forAll(CdRe, celli)
    {
        int i;

        if (Re[celli] < scalar(1.8))
        {
            i = 0;
            CdRe[celli] = ( Cd_M[i] + (((Cd_M[i+1] - Cd_M[i])/(Re_M[i+1] - Re_M[i]))*(Re[celli] - Re_M[i])) ) * Re[celli];
        }
        else if (Re[celli] >= scalar(3549))
        {
            i = 32;
            CdRe[celli] = ( Cd_M[i] + (((Cd_M[i+1] - Cd_M[i])/(Re_M[i+1] - Re_M[i]))*(Re[celli] - Re_M[i])) ) * Re[celli];
        }
        else
        {
            for (i=0; i<=32; ++i)
            {
                if ((Re_M[i] <= Re[celli]) && (Re[celli] < Re_M[i+1]))
                {
                    CdRe[celli] = ( Cd_M[i] + (((Cd_M[i+1] - Cd_M[i])/(Re_M[i+1] - Re_M[i]))*(Re[celli] - Re_M[i])) ) * Re[celli];
                    break;
                }
            }
        }
    }

    if (debug)
    {
        Info<<"min/max(Re) : " << min(Re) << " , " << max(Re) << endl;
        Info<<"min/max(CdRe) : " << min(CdRe) << " , " << max(CdRe) << endl;
    }

    CdRe.correctBoundaryConditions();
    return tCdRe;
}


Foam::tmp<Foam::volScalarField> Foam::decoupledEulerian::GunnKinzer::K() const
{
    if (autoMaterialProp_)
    {
        return
            0.75
           *CdRe()
           *mua(phase_.Td())
           /rhop(phase_.Td())
           /sqr(phase_.diam());
    }

    return
        0.75
       *CdRe()
       *phase_.muc()
       /phase_.rhod()
       /sqr(phase_.diam());
}


const Foam::dimensionedScalar Foam::decoupledEulerian::GunnKinzer::mua(const dimensionedScalar& T) const
{
    int i;
    dimensionedScalar mua("mua", phase_.muc().dimensions(), 0.0);

    double T_1_M[]={-18.0, -6.7, 0.0, 4.4, 15.6, 26.7, 37.8};
    double MU_A_M[] = {0.0157E-3, 0.0168E-3, 0.0171E-3, 0.0173E-3, 0.0179E-3, 0.0184E-3, 0.0190E-3};

    if (T.value()-scalar(273.15) < scalar(-18))
    {
        Info<< "Air temperature is too low!";
    }
    else if (T.value()-scalar(273.15) >= scalar(37.8))
    {
        Info<< "Air temperature is too high!";
    }
    else
    {
        for (i=0; i<=5; ++i)
        {
            if ((T_1_M[i] <= T.value()-scalar(273.15)) && (T.value()-scalar(273.15) < T_1_M[i+1]))
            {
                mua.value() = MU_A_M[i] + ( ((MU_A_M[i+1] - MU_A_M[i])/(T_1_M[i+1] - T_1_M[i])) * (T.value()-scalar(273.15) - T_1_M[i]));
                break;
            }
        }
    }

    return mua;
}


const Foam::dimensionedScalar Foam::decoupledEulerian::GunnKinzer::rhoa(const dimensionedScalar& T) const
{
    int i;
    dimensionedScalar rhoa("rhoa", dimDensity, 0.0);

    double T_1_M[]={-18.0, -6.7, 0.0, 4.4, 15.6, 26.7, 37.8};
    double RHO_A_M[] = {1.38, 1.32, 1.293, 1.27, 1.22, 1.18, 1.13};

    if (T.value()-scalar(273.15) < scalar(-18))
    {
        Info<< "Air temperature is too low!";
    }
    else if (T.value()-scalar(273.15) >= scalar(37.8))
    {
        Info<< "Air temperature is too high!";
    }
    else
    {
        for (i=0; i<=5; ++i)
        {
            if ((T_1_M[i] <= T.value()-scalar(273.15)) && (T.value()-scalar(273.15) < T_1_M[i+1]))
            {
                rhoa.value() = RHO_A_M[i] + ( ((RHO_A_M[i+1] - RHO_A_M[i])/(T_1_M[i+1] - T_1_M[i])) * (T.value()-scalar(273.15) - T_1_M[i]));
                break;
            }
        }
    }

    return rhoa;
}


const Foam::dimensionedScalar Foam::decoupledEulerian::GunnKinzer::rhop(const dimensionedScalar& T) const
{
    int i;
    dimensionedScalar rhop("rhop", dimDensity, 0.0);

    double T_2_M[] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0,
                20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0};
    double RHO_P_M[] = {999.87, 999.97, 1000.00, 999.97, 999.88, 999.73, 999.52, 999.27, 998.97, 998.62,
                998.23, 997.80, 997.33, 996.81, 996.26, 995.68, 994.00, 992.00};

    if (T.value()-scalar(273.15) <= scalar(0))
    {
        Info<< "Air temperature is too low!";
    }
    else if (T.value()-scalar(273.15) >= scalar(40))
    {
        Info<< "Air temperature is too high!";
    }
    else
    {
        for (i=0; i<=16; ++i)
        {
            if ((T_2_M[i] <= T.value()-scalar(273.15)) && (T.value()-scalar(273.15) < T_2_M[i+1]))
            {
                rhop.value() = RHO_P_M[i] + ( ((RHO_P_M[i+1] - RHO_P_M[i])/(T_2_M[i+1] - T_2_M[i])) * (T.value()-scalar(273.15) - T_2_M[i]));
                break;
            }
        }
    }

    return rhop;
}


// ************************************************************************* //
