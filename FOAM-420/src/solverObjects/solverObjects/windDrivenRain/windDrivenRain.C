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
    (c) 2017-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "windDrivenRain.H"
#include "solverOption/SolverOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvc.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/Fields/oneField/oneField.H"
#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "cfdTools/general/bound/bound.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(windDrivenRain, 0);
}
}

makeFvSolverOption(windDrivenRain);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::windDrivenRain::computeMaterialParam()
{
    int i;

    double T_1_M[] = {-18.0, -6.7, 0.0, 4.4, 15.6, 26.7, 37.8};
    double RHO_A_M[] = {1.38, 1.32, 1.293, 1.27, 1.22, 1.18, 1.13};
    double MU_A_M[] = {0.0157E-3, 0.0168E-3, 0.0171E-3, 0.0173E-3, 0.0179E-3, 0.0184E-3, 0.0190E-3};

    double T_2_M[] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0,
                      20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0};
    double RHO_P_M[] = {999.87, 999.97, 1000.00, 999.97, 999.88, 999.73, 999.52, 999.27, 998.97, 998.62,
                        998.23, 997.80, 997.33, 996.81, 996.26, 995.68, 994.00, 992.00};

    if (temp_.value()-scalar(273.15) < scalar(-18))
    {
        Info<< "Air temperature is too low!";
    }
    else if (temp_.value()-scalar(273.15) >= scalar(37.8))
    {
        Info<< "Air temperature is too high!";
    }
    else
    {
        for (i=0; i<=5; ++i)
        {
            if ((T_1_M[i] <= temp_.value()-scalar(273.15)) && (temp_.value()-scalar(273.15) < T_1_M[i+1]))
            {
                mua_.value() = MU_A_M[i] + ( ((MU_A_M[i+1] - MU_A_M[i])/(T_1_M[i+1] - T_1_M[i])) * (temp_.value()-scalar(273.15) - T_1_M[i]));
                rhoa_.value() = RHO_A_M[i] + ( ((RHO_A_M[i+1] - RHO_A_M[i])/(T_1_M[i+1] - T_1_M[i])) * (temp_.value()-scalar(273.15) - T_1_M[i]));
                break;
            }
        }
    }

    if (temp_.value()-scalar(273.15) <= scalar(0))
    {
        Info<< "Air temperature is too low!";
    }
    else if (temp_.value()-scalar(273.15) >= scalar(40))
    {
        Info<< "Air temperature is too high!";
    }
    else
    {
        for (i=0; i<=16; ++i)
        {
            if ((T_2_M[i] <= temp_.value()-scalar(273.15)) && (temp_.value()-scalar(273.15) < T_2_M[i+1]))
            {
                rhop_.value() = RHO_P_M[i] + ( ((RHO_P_M[i+1] - RHO_P_M[i])/(T_2_M[i+1] - T_2_M[i])) * (temp_.value()-scalar(273.15) - T_2_M[i]));
                break;
            }
        }
    }
}

void Foam::fv::windDrivenRain::computeCdRe()
{
    forAll(CdRe(), celli)
    {
        int i;

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

        if (Re()[celli] < scalar(1.8))
        {
            i = 0;
            CdRe()[celli] = ( Cd_M[i] + (((Cd_M[i+1] - Cd_M[i])/(Re_M[i+1] - Re_M[i]))*(Re()[celli] - Re_M[i])) ) * Re()[celli];
        }
        else if (Re()[celli] >= scalar(3549))
        {
            i = 32;
            CdRe()[celli] = ( Cd_M[i] + (((Cd_M[i+1] - Cd_M[i])/(Re_M[i+1] - Re_M[i]))*(Re()[celli] - Re_M[i])) ) * Re()[celli];
        }
        else
        {
            for (i=0; i<=32; ++i)
            {
                if ((Re_M[i] <= Re()[celli]) && (Re()[celli] < Re_M[i+1]))
                {
                    CdRe()[celli] = ( Cd_M[i] + (((Cd_M[i+1] - Cd_M[i])/(Re_M[i+1] - Re_M[i]))*(Re()[celli] - Re_M[i])) ) * Re()[celli];
                    break;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::windDrivenRain::windDrivenRain
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),

    Rh_(dimensionedScalar("Rh", dimVelocity, 0.)),
    temp_(dimensionedScalar("temp", dimTemperature, 0.)),
    rhoa_(dimensionedScalar("rhoa", dimMass/dimVolume, 0.)),
    rhop_(dimensionedScalar("rhob", dimMass/dimVolume, 0.)),
    mua_(dimensionedScalar("mua", dimMass/dimLength/dimTime, 0.)),
    scalingFactor_(),
    stabDivAlpha_(dict.lookupOrDefault<scalar>("stabDivAlpha", 0.001)),

    solveTD_(),
    phases_(),

    UPtr_(nullptr),
    //phiwindPtr_(),
    RePtr_(),
    CdRePtr_(),

    CtrainPtrL_(),
    kPtr_(),
    epsilonPtr_(),
    nutPtr_(),
    nutrainPtr_(),

    UrainPtrL_(),
    phirainPtrL_(),
    alpharainPtrL_(),
    scrPtrL_(),

    nNonOrthCorr_(0),

    fieldDependency_(dict.lookupOrDefault<word>("fieldDependency", "none"))
{
    read(dict);
    #include "createFieldsWdr.H"
    #include "createRainFieldsWdr.H"
    #include "createTDFieldsWdr.H"
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::windDrivenRain::~windDrivenRain()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::windDrivenRain::read(const dictionary& dict)
{
    nNonOrthCorr_ = dict.lookupOrDefault("nNonOrthogonalCorrectors", 0);
    dict.readIfPresent("fieldDependency", fieldDependency_);
}


bool Foam::fv::windDrivenRain::initialise()
{
    return true;
}


void Foam::fv::windDrivenRain::correct
(
    const word& solveName,
    const word& regionName
)
{
    if (debug)
    {
        Info<< "    " << "Solving for decoupled Eulerian particles. " << endl;
    }

    const Time& runTime(mesh_.time());
    const fvMesh& mesh(mesh_);

    #include "cfdTools/general/include/readGravitationalAcceleration.H"

    // compute density and viscosity based on temperature
    computeMaterialParam();

    // for all non-orth corrcetions
    for (int nonOrth=0; nonOrth <= nNonOrthCorr_; nonOrth++)
    {
        // for all rain fields
        for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
        {
            // need a field named "phi" for boundary conditions and Courant number calculation
            surfaceScalarField phi(phirain(phase_no));

            #include "cfdTools/incompressible/CourantNo.H"
            #include "alphaEqnsWdr.H"
            #include "computePhysicalQuantitiesWdr.H"
            #include "UEqnsWdr.H"
        }
    }

    // write result
    if (runTime.write())
    {
        for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
        {
            Urain(phase_no).write();
            alpharain(phase_no).write();
        }
    }

    //#include "calculateCatchRatioWdr.H"
}


// ************************************************************************* //
