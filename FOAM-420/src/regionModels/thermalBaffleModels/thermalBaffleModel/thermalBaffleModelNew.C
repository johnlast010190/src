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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermalBaffleModel/thermalBaffleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermalBaffleModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<thermalBaffleModel> thermalBaffleModel::New(const fvMesh& mesh)
{
    word modelType;
    {
        IOdictionary thermalBafflePropertiesDict
        (
            IOobject
            (
                "thermalBaffleProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        word modelType =
            thermalBafflePropertiesDict.lookupOrDefault<word>
            (
                "thermalBaffleModel",
                "thermalBaffle"
            );
    }

    const auto ctor = ctorTableLookup("thermalBaffleModel type", meshConstructorTable_(), modelType);
    return autoPtr<thermalBaffleModel>(ctor(modelType, mesh));
}


autoPtr<thermalBaffleModel> thermalBaffleModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word modelType =
        dict.lookupOrDefault<word>("thermalBaffleModel", "thermalBaffle");

    const auto ctor = ctorTableLookup("thermalBaffleModel type", dictionaryConstructorTable_(), modelType);
    return autoPtr<thermalBaffleModel>(ctor(modelType, mesh, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalBaffleModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
