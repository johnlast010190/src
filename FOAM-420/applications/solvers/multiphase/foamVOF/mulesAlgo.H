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
    (c) ESI Ltd, 2019.

Description :

  This file defines the MulesAlgo class in FoamVOF model implementation.
  FoamVOF is a general-purpose multi-fluid module developed under the
  framework of volume-of-fluid methodology.
\*---------------------------------------------------------------------------*/

#include "mulesAlgo.h"

using namespace Foam;

MulesAlgo::MulesAlgo
(
    const Time& time,
    const dynamicFvMesh& mesh,
    const MatProp &prop
)
:
FoamVOF
(
    time,
    mesh,
    prop
)
{
    name_("Muels");
}

void MulesAlgo::updateFields
(
    pimpleControl &pimple,
    incompressible::turbulenceModel &turbulence,
    fv::options& fvOptions
)
{
    readTimeControls();
    updateAcceleration(frameAcceleration);

    if (LTS)
    {
        setRDeltaT(pimple);
    }
    else
    {
        CourantNo();
        alphaCourantNoIF();
        setDeltaT();
    }

    const_cast<Time&> (runTime_)++;

    Info<< "Time = " << runTime().timeName() << nl << endl;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
      if (!isA<staticFvMesh>(mesh()))
      {
         if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
           dynamicMeshUpdate(pimple);
        }
      }

      alphaControls();

      alphaEqnSubCycle();
      mixture->correct();

      if (pimple.frozenFlow())
      {
        continue;
      }

      fvVectorMatrix UEqn(U(), dimensionSet(1, 1, -2, 0, 0, 0, 0));
      UEqnSolve(UEqn,turbulence,fvOptions,mixture(),pimple);

      // --- Pressure corrector loop
      while (pimple.correct())
      {
           pEqn(UEqn, fvOptions,mixture(),pimple);
      }

     if (pimple.turbCorr())
     {
           turbulence.correct();
     }
   }


   return;

}


void MulesAlgo::createFields
(
    pimpleControl &pimple
)
{
    FoamVOF::createFields(pimple);
    return;
}

