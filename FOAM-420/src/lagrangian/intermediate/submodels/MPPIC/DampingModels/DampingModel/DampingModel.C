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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/DampingModels/DampingModel/DampingModel.H"
#include "submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DampingModel<CloudType>::DampingModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    timeScaleModel_(nullptr)
{}


template<class CloudType>
Foam::DampingModel<CloudType>::DampingModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    timeScaleModel_
    (
        TimeScaleModel::New
        (
            this->coeffDict().subDict(TimeScaleModel::typeName)
        )
    )
{}


template<class CloudType>
Foam::DampingModel<CloudType>::DampingModel(const DampingModel<CloudType>& cm)
:
    CloudSubModelBase<CloudType>(cm),
    timeScaleModel_(cm.timeScaleModel_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DampingModel<CloudType>::~DampingModel()
{}


// * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::DampingModel<CloudType>>
Foam::DampingModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting damping model " << modelType << endl;

    const auto ctor = ctorTableLookup("damping model type", dictionaryConstructorTable_(), modelType);
    return autoPtr<DampingModel<CloudType>>(ctor(dict, owner));
}


// ************************************************************************* //
