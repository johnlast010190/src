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

#include "kappaSutherland.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kappaSutherland, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        kappaSutherland,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kappaSutherland::kappaSutherland
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
    kappaSutherland::read();
}


Foam::autoPtr<Foam::kappaSutherland>
Foam::kappaSutherland::clone() const
{
    return autoPtr<kappaSutherland>
    (
        new kappaSutherland(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kappaSutherland::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    // Model names
    const word& kappaName = kappaModel::typeName;
    const word& CvName = CvModel::typeName;
    const word& muName = muModel::typeName;
    const word& RName = RModel::typeName;

    if (modelName == kappaName)
    {
        sMod_.set(Cv, models[CvName]);
        sMod_.set(mu, models[muName]);
        sMod_.set(R, models[RName]);
        dep_[0].model = models[kappaName];
        dep_[0].dependencies.setSize(3);
        dep_[0].dependencies.set(0, models[CvName]);
        dep_[0].dependencies.set(1, models[muName]);
        dep_[0].dependencies.set(2, models[RName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::kappaSutherland::castScalarModel
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


Foam::tmp<Foam::scalarField> Foam::kappaSutherland::kappaPatch
(
    const label patchi
) const
{
    tmp<scalarField> CvThermo = sMod_[Cv].boundaryField()[patchi];

    return
        sMod_[mu].boundaryField()[patchi]
       *CvThermo.ref()
       *(1.32 + 1.77*sMod_[R][0]/CvThermo.ref());
}


Foam::tmp<Foam::scalarField> Foam::kappaSutherland::kappaInternal() const
{
    tmp<scalarField> CvThermo = sMod_[Cv].primitiveField();

    return
        sMod_[mu].primitiveField()
       *CvThermo.ref()
       *(1.32 + 1.77*sMod_[R][0]/CvThermo.ref());
}


Foam::tmp<Foam::volScalarField>
Foam::kappaSutherland::kappaGeometric() const
{
    tmp<volScalarField> CvThermo = sMod_[Cv]();
    dimensionedScalar RR
    (
        "R",
        dimEnergy/(dimMass*dimTemperature),
        sMod_[R][0]
    );
    return
        sMod_[mu]()
       *CvThermo.ref()
       *(1.32 + 1.77*RR/CvThermo.ref());
}


Foam::scalar Foam::kappaSutherland::kappaCell(const label celli) const
{
    const scalar cellCv = sMod_[Cv][celli];
    return sMod_[mu][celli]*cellCv*(1.32 + 1.77*sMod_[R][0]/cellCv);
}


bool Foam::kappaSutherland::read()
{
    return true;
}


// ************************************************************************* //
