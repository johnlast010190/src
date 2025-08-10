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

Class
    Foam::materialModels

Description
    Interface classes.

SourceFiles
    materialModels.C

\*---------------------------------------------------------------------------*/


#include "fields/volFields/volFields.H"
#include "baseModels.H"
#include "materialModels/materialMacros.H"

#ifndef WModel_H
#define WModel_H
    defineScalarMaterialModel(W);
#endif

#ifndef YModel_H
#define YModel_H
    defineScalarMaterialModel(Y);
#endif

#ifndef RModel_H
#define RModel_H
    defineScalarMaterialModel(R);
#endif

#ifndef limitModel_H
#define limitModel_H
    defineScalarMaterialModel(limit);
#endif

#ifndef rhoModel_H
#define rhoModel_H
    defineScalarMaterialModel(rho);
#endif

#ifndef CpDepartureModel_H
#define CpDepartureModel_H
    defineScalarMaterialModel(CpDeparture);
#endif

#ifndef HDepartureModel_H
#define HDepartureModel_H
    defineScalarMaterialModel(HDeparture);
#endif

#ifndef psiModel_H
#define psiModel_H
    defineScalarMaterialModel(psi);
#endif

#ifndef ZModel_H
#define ZModel_H
    defineScalarMaterialModel(Z);
#endif

#ifndef CpMCvModel_H
#define CpMCvModel_H
    defineScalarMaterialModel(CpMCv);
#endif

#ifndef SDepartureModel_H
#define SDepartureModel_H
    defineScalarMaterialModel(SDeparture);
#endif

#ifndef CpThermoModel_H
#define CpThermoModel_H
    defineScalarMaterialModel(CpThermo);
#endif

#ifndef HaThermoModel_H
#define HaThermoModel_H
    defineScalarMaterialModel(HaThermo);
#endif

#ifndef HsThermoModel_H
#define HsThermoModel_H
    defineScalarMaterialModel(HsThermo);
#endif

#ifndef HcModel_H
#define HcModel_H
    defineScalarMaterialModel(Hc);
#endif

#ifndef SThermoModel_H
#define SThermoModel_H
    defineScalarMaterialModel(SThermo);
#endif

#ifndef AModel_H
#define AModel_H
    defineScalarMaterialModel(A);
#endif

#ifndef CpByCpvModel_H
#define CpByCpvModel_H
    defineScalarMaterialModel(CpByCpv);
#endif

#ifndef CpModel_H
#define CpModel_H
    defineScalarMaterialModel(Cp);
#endif

#ifndef CpvModel_H
#define CpvModel_H
    defineScalarMaterialModel(Cpv);
#endif

#ifndef CvModel_H
#define CvModel_H
    defineScalarMaterialModel(Cv);
#endif

#ifndef EaModel_H
#define EaModel_H
    defineScalarMaterialModel(Ea);
#endif

#ifndef EsModel_H
#define EsModel_H
    defineScalarMaterialModel(Es);
#endif

#ifndef gammaModel_H
#define gammaModel_H
    defineScalarMaterialModel(gamma);
#endif

#ifndef GModel_H
#define GModel_H
    defineScalarMaterialModel(G);
#endif

#ifndef HaModel_H
#define HaModel_H
    defineScalarMaterialModel(Ha);
#endif

#ifndef HEModel_H
#define HEModel_H
    defineScalarMaterialModel(HE);
#endif

#ifndef HsModel_H
#define HsModel_H
    defineScalarMaterialModel(Hs);
#endif

#ifndef SModel_H
#define SModel_H
    defineScalarMaterialModel(S);
#endif

#ifndef strainRateModel_H
#define strainRateModel_H
    defineScalarMaterialModel(strainRate);
#endif

#ifndef muModel_H
#define muModel_H
    defineScalarMaterialModel(mu);
#endif

#ifndef kappaModel_H
#define kappaModel_H
    defineScalarMaterialModel(kappa);
#endif

#ifndef vKappaModel_H
#define vKappaModel_H
    defineVectorMaterialModel(vKappa);
#endif

#ifndef alphahModel_H
#define alphahModel_H
    defineScalarMaterialModel(alphah);
#endif

#ifndef vAlphahModel_H
#define vAlphahModel_H
    defineVectorMaterialModel(vAlphah);
#endif

#ifndef DModel_H
#define DModel_H
    defineScalarMaterialModel(D);
#endif

#ifndef weightedAverageMixtureModel_H
#define weightedAverageMixtureModel_H
    defineScalarMaterialModel(weightedAverageMixture);
#endif

#ifndef vectorWeightedAverageMixtureModel_H
#define vectorWeightedAverageMixtureModel_H
    defineVectorMaterialModel(vectorWeightedAverageMixture);
#endif

#ifndef tensorWeightedAverageMixtureModel_H
#define tensorWeightedAverageMixtureModel_H
    defineTensorMaterialModel(tensorWeightedAverageMixture);
#endif

#ifndef clippedWeightedAverageMixtureModel_H
#define clippedWeightedAverageMixtureModel_H
    defineScalarMaterialModel(clippedWeightedAverageMixture);
#endif

#ifndef vectorClippedWeightedAverageMixtureModel_H
#define vectorClippedWeightedAverageMixtureModel_H
    defineVectorMaterialModel(vectorClippedWeightedAverageMixture);
#endif

#ifndef tensorClippedWeightedAverageMixtureModel_H
#define tensorClippedWeightedAverageMixtureModel_H
    defineTensorMaterialModel(tensorClippedWeightedAverageMixture);
#endif

#ifndef harmonicMixtureModel_H
#define harmonicMixtureModel_H
    defineScalarMaterialModel(harmonicMixture);
#endif

#ifndef TModel_H
#define TModel_H
    defineScalarMaterialModel(T);
#endif

#ifndef THaModel_H
#define THaModel_H
    defineScalarMaterialModel(THa);
#endif

#ifndef THsModel_H
#define THsModel_H
    defineScalarMaterialModel(THs);
#endif

#ifndef TEaModel_H
#define TEaModel_H
    defineScalarMaterialModel(TEa);
#endif

#ifndef TEsModel_H
#define TEsModel_H
    defineScalarMaterialModel(TEs);
#endif

#ifndef pModel_H
#define pModel_H
    defineScalarMaterialModel(p);
#endif

#ifndef testVectorModel_H
#define testVectorModel_H
    defineVectorMaterialModel(testVector);
#endif

#ifndef testTensorModel_H
#define testTensorModel_H
    defineTensorMaterialModel(testTensor);
#endif

// ************************************************************************* //
