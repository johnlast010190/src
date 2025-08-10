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

\*---------------------------------------------------------------------------*/

#include "IAPWSIF97EOSModels.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IAPWSIF97EOSModels, 0);
    addToRunTimeSelectionTable(materialModel, IAPWSIF97EOSModels, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IAPWSIF97EOSModels::IAPWSIF97EOSModels
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    waterSteamModel_(obr)
{
    IAPWSIF97EOSModels::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::IAPWSIF97EOSModels::castScalarModel
(
    const word& modelName
)
{
    if (modelName == rhoModel::typeName)
    {
        return dynamic_cast<rhoModel*>(this);
    }
    else if (modelName == psiModel::typeName)
    {
        return dynamic_cast<psiModel*>(this);
    }
    else if (modelName == ZModel::typeName)
    {
        return dynamic_cast<ZModel*>(this);
    }
    else if (modelName == CpMCvModel::typeName)
    {
        return dynamic_cast<CpMCvModel*>(this);
    }
    else if (modelName == CpModel::typeName)
    {
        return dynamic_cast<CpModel*>(this);
    }
    else if (modelName == CvModel::typeName)
    {
        return dynamic_cast<CvModel*>(this);
    }
    else if (modelName == HaModel::typeName)
    {
        return dynamic_cast<HaModel*>(this);
    }
    else if (modelName == HsModel::typeName)
    {
        return dynamic_cast<HsModel*>(this);
    }
    else if (modelName == HcModel::typeName)
    {
        return dynamic_cast<HcModel*>(this);
    }
    else if (modelName == SModel::typeName)
    {
        return dynamic_cast<SModel*>(this);
    }
    else if (modelName == muModel::typeName)
    {
        return dynamic_cast<muModel*>(this);
    }
    else if (modelName == kappaModel::typeName)
    {
        return dynamic_cast<kappaModel*>(this);
    }
    return nullptr;
}


bool Foam::IAPWSIF97EOSModels::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    p_ = constructOrReturnRefFieldPtr<scalar>("p");
    Hf_ = dict_->lookupOrDefault("Hf", 0.0);
    Sf_ = dict_->lookupOrDefault("Sf", 0.0);

    // This can cause issues with unstable runs
    // (artificial jumps in temperature in other region) or in not suitable
    // initial condions. Hence for now it will be by default allowed to change
    // phase regions. However, if the simulation actually suppose to simulate
    // phase change in the domain this will not be correct for solvers not
    // supporting phase change.
    supportPhaseChange_ = dict_->lookupOrDefault("supportPhaseChange", true);

    // We assume that the whole domain should be in one region when the phase
    // change isn't supported.
    regionNumber_ =
        p_->primitiveField()().size()
      ? waterSteamModel_.whichRegion
        (
            p_->primitiveField()().first(),
            T_->primitiveField()().first()
        )
      : 0;
    reduce(regionNumber_, maxOp<label>());
    return true;
}


// ************************************************************************* //
