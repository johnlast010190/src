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

#include "submodels/MPPIC/PackingModels/PackingModel/PackingModel.H"
#include "submodels/MPPIC/AveragingMethods/AveragingMethod/AveragingMethod.H"
#include "submodels/MPPIC/ParticleStressModels/ParticleStressModel/ParticleStressModel.H"
#include "submodels/MPPIC/CorrectionLimitingMethods/CorrectionLimitingMethod/CorrectionLimitingMethod.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModel<CloudType>::PackingModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    particleStressModel_(nullptr)
{}


template<class CloudType>
Foam::PackingModel<CloudType>::PackingModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    particleStressModel_
    (
        ParticleStressModel::New
        (
            this->coeffDict().subDict(ParticleStressModel::typeName)
        )
    )
{}


template<class CloudType>
Foam::PackingModel<CloudType>::PackingModel(const PackingModel<CloudType>& cm)
:
    CloudSubModelBase<CloudType>(cm),
    particleStressModel_(cm.particleStressModel_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModel<CloudType>::~PackingModel()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::PackingModel<CloudType>>
Foam::PackingModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting packing model " << modelType << endl;

    const auto ctor = ctorTableLookup("packing model type", dictionaryConstructorTable_(), modelType);
    return autoPtr<PackingModel<CloudType>>(ctor(dict, owner));
}


// ************************************************************************* //
