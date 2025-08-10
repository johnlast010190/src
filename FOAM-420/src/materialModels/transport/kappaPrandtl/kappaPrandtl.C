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

#include "kappaPrandtl.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kappaPrandtl, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        kappaPrandtl,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kappaPrandtl::kappaPrandtl
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
    dep_.setSize(1);
    kappaPrandtl::read();
}


Foam::autoPtr<Foam::kappaPrandtl>
Foam::kappaPrandtl::clone() const
{
    return autoPtr<kappaPrandtl>
    (
        new kappaPrandtl(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kappaPrandtl::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    // Model names
    const word& kappaName = kappaModel::typeName;
    const word& CpName = CpModel::typeName;
    const word& muName = muModel::typeName;

    if (modelName == kappaName)
    {
        sMod_.set(Cp, models[CpName]);
        sMod_.set(muMod, models[muName]);
        dep_[0].model = models[kappaName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, models[CpName]);
        dep_[0].dependencies.set(1, models[muName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::kappaPrandtl::castScalarModel
(
    const word& modelName
)
{
    if (modelName == kappaModel::typeName)
    {
        return dynamic_cast<kappaModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::kappaPrandtl::kappaPatch
(
    const label patchi
) const
{
    return
        sMod_[muMod].boundaryField()[patchi]
       *sMod_[Cp].boundaryField()[patchi]
       *rPr_;
}


Foam::tmp<Foam::scalarField> Foam::kappaPrandtl::kappaInternal() const
{
    return
        sMod_[muMod].primitiveField()*sMod_[Cp].primitiveField()*rPr_;
}


Foam::tmp<Foam::volScalarField>
Foam::kappaPrandtl::kappaGeometric() const
{
    return sMod_[muMod]()*sMod_[Cp]()*rPr_;
}


Foam::scalar Foam::kappaPrandtl::kappaCell(const label celli) const
{
    return sMod_[muMod][celli]*sMod_[Cp][celli]*rPr_;
}


bool Foam::kappaPrandtl::read()
{
    rPr_ = 1.0/dict_->lookup<scalar>("Pr");
    return true;
}


// ************************************************************************* //
