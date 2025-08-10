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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "materialsSpecie.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materialsSpecie, 0);
    addToRunTimeSelectionTable(materialModel, materialsSpecie, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialsSpecie::materialsSpecie
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& materialsSpecieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, materialsSpecieName, name)
{
    materialsSpecie::read();
    if (RR == 0.0 || Pstd == 0.0 || Tstd == 0.0)
    {
        const_cast<scalar&>(RR) = constant::physicoChemical::R.value()*1000;
        const_cast<scalar&>(Pstd) = constant::standard::Pstd.value();
        const_cast<scalar&>(Tstd) = constant::standard::Tstd.value();
        WarningInFunction
            << "Static constants not initialised!"
            << " Attmpting in-situ assignment."
            << nl << "Specie static variable values:"
            << nl << tab << "RR = " << RR
            << nl << tab << "Pstd = " << Pstd
            << nl << tab << "Tstd = " << Tstd << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::materialsSpecie::castScalarModel
(
    const word& modelName
)
{
    if (modelName == WModel::typeName)
    {
        return dynamic_cast<WModel*>(this);
    }
    else if (modelName == YModel::typeName)
    {
        return dynamic_cast<YModel*>(this);
    }
    else if (modelName == RModel::typeName)
    {
        return dynamic_cast<RModel*>(this);
    }
    return nullptr;
}


bool Foam::materialsSpecie::read()
{
    Y_ = dict_->lookupOrDefault("massFraction", 1.0);
    molWeight_ = dict_->lookup<scalar>("molWeight");
    R_ = dict_->lookupOrDefault("R", RR/molWeight_);

    return true;
}


// ************************************************************************* //
