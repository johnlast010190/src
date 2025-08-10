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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "materialModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    declareMaterialModelScalarDims(W, dimMass/dimMoles);
    declareMaterialModelScalarDims(R, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Y, dimless);
    declareMaterialModelScalarDims(limit, dimless);
    declareMaterialModelScalarDims(rho, dimMass/dimVolume);
    declareMaterialModelScalarDims(SDeparture, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(CpDeparture, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(HDeparture, dimEnergy/dimMass);
    declareMaterialModelScalarDims(psi, sqr(dimTime)/dimArea);
    declareMaterialModelScalarDims(Z, dimless);
    declareMaterialModelScalarDims(CpMCv, dimEnergy/dimMass/dimTemperature);

    declareMaterialModelScalarDims(CpThermo, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(HaThermo, dimEnergy/dimMass);
    declareMaterialModelScalarDims(HsThermo, dimEnergy/dimMass);
    declareMaterialModelScalarDims(Hc, dimEnergy/dimMass);
    declareMaterialModelScalarDims(SThermo, dimEnergy/dimMass/dimTemperature);

    declareMaterialModelScalarDims(A, dimEnergy/dimMass);
    declareMaterialModelScalarDims(CpByCpv, dimless);
    declareMaterialModelScalarDims(Cp, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Cpv, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Cv, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Ea, dimEnergy/dimMass);
    declareMaterialModelScalarDims(Es, dimEnergy/dimMass);
    declareMaterialModelScalarDims(gamma, dimless);
    declareMaterialModelScalarDims(G, dimEnergy/dimMass);
    declareMaterialModelScalarDims(Ha,dimEnergy/dimMass);
    declareMaterialModelScalarDims(HE, dimEnergy/dimMass);
    declareMaterialModelScalarDims(Hs, dimEnergy/dimMass);
    declareMaterialModelScalarDims(S, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(strainRate, dimless/dimTime);
    declareMaterialModelScalarDims(mu, dimDynamicViscosity);
    declareMaterialModelScalarDims(kappa, dimPower/(dimLength*dimTemperature));
    declareMaterialModelVectorDims(vKappa, dimPower/(dimLength*dimTemperature));
    declareMaterialModelScalarDims(alphah, dimDynamicViscosity);
    declareMaterialModelVectorDims(vAlphah, dimDynamicViscosity);
    declareMaterialModelScalarDims(D, dimArea/dimTime);

    // Mixture models
    declareMaterialModelScalarMixture(weightedAverageMixture, dimless);
    declareMaterialModelVectorMixture(vectorWeightedAverageMixture, dimless);
    declareMaterialModelTensorMixture(tensorWeightedAverageMixture, dimless);
    declareMaterialModelScalarMixture(clippedWeightedAverageMixture, dimless);
    declareMaterialModelVectorMixture(vectorClippedWeightedAverageMixture, dimless);
    declareMaterialModelTensorMixture(tensorClippedWeightedAverageMixture, dimless);
    declareMaterialModelScalarMixture(harmonicMixture, dimless);

    declareMaterialModelScalarDims(T, dimTemperature);
    declareMaterialModelScalarDims(THa, dimTemperature);
    declareMaterialModelScalarDims(THs, dimTemperature);
    declareMaterialModelScalarDims(TEa, dimTemperature);
    declareMaterialModelScalarDims(TEs, dimTemperature);
    declareMaterialModelScalarDims(p, dimPressure);

    // Dummy vector model
    declareMaterialModelVectorDims(testVector, dimless);
    declareMaterialModelTensorDims(testTensor, dimless);

}


// ************************************************************************* //
