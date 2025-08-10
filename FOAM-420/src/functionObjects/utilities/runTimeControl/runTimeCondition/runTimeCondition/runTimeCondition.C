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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(runTimeCondition, 0);
    defineRunTimeSelectionTable(runTimeCondition, dictionary);
}
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::dictionary&
Foam::functionObjects::runTimeControls::runTimeCondition::setConditionDict()
{
    dictionary& propertyDict = state_.propertyDict();

    if (!propertyDict.found(name_))
    {
        propertyDict.add(name_, dictionary());
    }

    return propertyDict.subDict(name_);
}


const Foam::dictionary&
Foam::functionObjects::runTimeControls::runTimeCondition::conditionDict() const
{
    return conditionDict_;
}


Foam::dictionary&
Foam::functionObjects::runTimeControls::runTimeCondition::conditionDict()
{
    return conditionDict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeCondition::runTimeCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    name_(name),
    obr_(obr),
    state_(state),
    active_(dict.lookupOrDefault<bool>("active", true)),
    conditionDict_(setConditionDict()),
    groupID_(dict.lookupOrDefault("groupID", -1)),
    log(state.log)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeCondition::~runTimeCondition()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

const Foam::word&
Foam::functionObjects::runTimeControls::runTimeCondition::name() const
{
    return name_;
}


bool Foam::functionObjects::runTimeControls::runTimeCondition::active() const
{
    return active_;
}


Foam::label
Foam::functionObjects::runTimeControls::runTimeCondition::groupID() const
{
    return groupID_;
}


// ************************************************************************* //
