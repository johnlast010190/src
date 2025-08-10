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

#include "vAlphahKappaOverCp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vAlphahKappaOverCp, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        vAlphahKappaOverCp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vAlphahKappaOverCp::vAlphahKappaOverCp
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
    vMod_.setSize(modelsEnumSize2_);
    dep_.setSize(2);
    vAlphahKappaOverCp::read();
}


Foam::autoPtr<Foam::vAlphahKappaOverCp>
Foam::vAlphahKappaOverCp::clone() const
{
    return autoPtr<vAlphahKappaOverCp>
    (
        new vAlphahKappaOverCp(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vAlphahKappaOverCp::updateTable(const word& modelName)
{
    const matScalarTable& sModels =
        materialTables_.sTable(phaseName_, specieName_);
    const matVectorTable& vModels =
        materialTables_.vTable(phaseName_, specieName_);

    if (modelName == vAlphahModel::typeName)
    {
        sMod_.set(Cp, sModels[CpModel::typeName]);
        vMod_.set(vKappaInd, vModels[vKappaModel::typeName]);

        dep_[0].model = vModels[vAlphahModel::typeName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, sModels[CpModel::typeName]);
        dep_[0].dependencies.set(1, vModels[vKappaModel::typeName]);
    }
}


Foam::baseModels<Foam::vector>* Foam::vAlphahKappaOverCp::castVectorModel
(
    const word& modelName
)
{
    if (modelName == vAlphahModel::typeName)
    {
        return dynamic_cast<vAlphahModel*>(this);
    }
    return nullptr;
}


Foam::vector Foam::vAlphahKappaOverCp::vAlphahCell
(
    const label celli
) const
{
    return vMod_[vKappaInd][celli]/sMod_[Cp][celli];
}


Foam::tmp<Foam::vectorField> Foam::vAlphahKappaOverCp::vAlphahPatch
(
    const label patchi
) const
{
    return
        vMod_[vKappaInd].boundaryField()[patchi]
       /sMod_[Cp].boundaryField()[patchi];
}


Foam::tmp<Foam::vectorField> Foam::vAlphahKappaOverCp::vAlphahInternal() const
{
    return
        vMod_[vKappaInd].primitiveField()
       /sMod_[Cp].primitiveField();
}


Foam::tmp<Foam::volVectorField>
Foam::vAlphahKappaOverCp::vAlphahGeometric() const
{
    return vMod_[vKappaInd]()/sMod_[Cp]();
}


bool Foam::vAlphahKappaOverCp::read()
{
    return true;
}


// ************************************************************************* //
