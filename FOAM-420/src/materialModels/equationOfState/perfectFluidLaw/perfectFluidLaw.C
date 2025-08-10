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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "perfectFluidLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(perfectFluidLaw, 0);
    addToRunTimeSelectionTable(materialModel, perfectFluidLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::perfectFluidLaw::perfectFluidLaw
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(4);
    perfectFluidLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::perfectFluidLaw::incompressible() const
{
    if (p_->isConst())
    {
        return true;
    }
    return false;
}


bool Foam::perfectFluidLaw::isochoric() const
{
    if (T_->isConst() && p_->isConst())
    {
        return true;
    }
    return false;
}


void Foam::perfectFluidLaw::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    const word tName = materialTables_.tableName(phaseName_, specieName_);

    // Model names
    const word& RName = RModel::typeName;
    const word& SDepartureName = SDepartureModel::typeName;
    const word& psiName = psiModel::typeName;
    const word& rhoName = rhoModel::typeName;
    const word& CpMCvName = CpMCvModel::typeName;

    if (modelName == rhoName)
    {
        sMod_.set(R, models[RName]);
        dep_[0].model = models[rhoName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, models[RName]);
        // Faster access R
        models[RName]->updateTable(RName);
        R_ = models[RName]->operator[](0);
    }
    else if (modelName == SDepartureName)
    {
        sMod_.set(R, models[RName]);
        dep_[1].model = models[SDepartureName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[RName]);
        // Faster access R
        models[RName]->updateTable(RName);
        R_ = models[RName]->operator[](0);
    }
    else if (modelName == psiName)
    {
        sMod_.set(R, models[RName]);
        dep_[2].model = models[psiName];
        dep_[2].dependencies.setSize(1);
        dep_[2].dependencies.set(0, models[RName]);
        // Faster access R
        models[RName]->updateTable(RName);
        R_ = models[RName]->operator[](0);
    }
    else if (modelName == CpMCvName)
    {
        sMod_.set(R, models[RName]);
        dep_[3].model = models[CpMCvName];
        dep_[3].dependencies.setSize(1);
        dep_[3].dependencies.set(0, models[RName]);
        // Faster access R
        models[RName]->updateTable(RName);
        R_ = models[RName]->operator[](0);
    }
}


Foam::baseModels<Foam::scalar>* Foam::perfectFluidLaw::castScalarModel
(
    const word& modelName
)
{
    if (modelName == rhoModel::typeName)
    {
        return dynamic_cast<rhoModel*>(this);
    }
    else if (modelName == SDepartureModel::typeName)
    {
        return dynamic_cast<SDepartureModel*>(this);
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
    return nullptr;
}


bool Foam::perfectFluidLaw::read()
{
    p_ = constructOrReturnRefFieldPtr<scalar>("p");
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    if (sMod_(R) != nullptr)
    {
        sMod_[R].read();
        sMod_[R].updateTable(RModel::typeName);
        R_ = sMod_[R][0];
    }

    rho0_ = dict_->lookup<scalar>("rho0");

    return true;
}


// ************************************************************************* //
