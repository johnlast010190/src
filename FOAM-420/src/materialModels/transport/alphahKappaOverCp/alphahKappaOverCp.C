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

#include "alphahKappaOverCp.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(alphahKappaOverCp, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        alphahKappaOverCp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphahKappaOverCp::alphahKappaOverCp
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
    alphahKappaOverCp::read();
}


Foam::autoPtr<Foam::alphahKappaOverCp>
Foam::alphahKappaOverCp::clone() const
{
    return autoPtr<alphahKappaOverCp>
    (
        new alphahKappaOverCp(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphahKappaOverCp::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    // Model names
    const word& kappaName = kappaModel::typeName;
    const word& alphahName = alphahModel::typeName;
    const word& CpName = CpModel::typeName;

    if (modelName == alphahName)
    {
        sMod_.set(Cp, models[CpName]);
        sMod_.set(kappaInd, models[kappaName]);
        dep_[0].model = models[alphahName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, models[CpName]);
        dep_[0].dependencies.set(1, models[kappaName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::alphahKappaOverCp::castScalarModel
(
    const word& modelName
)
{
    if (modelName == alphahModel::typeName)
    {
        return dynamic_cast<alphahModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::alphahKappaOverCp::alphahCell
(
    const label celli
) const
{
    return sMod_[kappaInd][celli]/sMod_[Cp][celli];
}


Foam::tmp<Foam::scalarField> Foam::alphahKappaOverCp::alphahPatch
(
    const label patchi
) const
{
    return
        sMod_[kappaInd].boundaryField()[patchi]
       /sMod_[Cp].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::alphahKappaOverCp::alphahInternal() const
{
    return
        sMod_[kappaInd].primitiveField()
       /sMod_[Cp].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::alphahKappaOverCp::alphahGeometric() const
{
    return sMod_[kappaInd]()/sMod_[Cp]();
}


bool Foam::alphahKappaOverCp::read()
{
    return true;
}


// ************************************************************************* //
