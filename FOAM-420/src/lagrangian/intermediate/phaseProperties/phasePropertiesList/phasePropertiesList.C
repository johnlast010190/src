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

#include "phaseProperties/phasePropertiesList/phasePropertiesList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePropertiesList::phasePropertiesList()
:
    props_(),
    phaseTypeNames_(),
    stateLabels_()
{}


Foam::phasePropertiesList::phasePropertiesList
(
    Istream& is,
    const wordList& gasNames,
    const wordList& liquidNames,
    const wordList& solidNames
)
:
    props_(is),
    phaseTypeNames_(),
    stateLabels_()
{
    forAll(props_, i)
    {
        props_[i].reorder(gasNames, liquidNames, solidNames);
    }

    phaseTypeNames_.setSize(props_.size());
    stateLabels_.setSize(props_.size());
    forAll(props_, i)
    {
        phaseTypeNames_[i] = props_[i].phaseTypeName();
        stateLabels_[i] = props_[i].stateLabel();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePropertiesList::~phasePropertiesList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::phaseProperties>&
Foam::phasePropertiesList::props() const
{
    return props_;
}


const Foam::wordList& Foam::phasePropertiesList::phaseTypes() const
{
    return phaseTypeNames_;
}


const Foam::wordList& Foam::phasePropertiesList::stateLabels() const
{
    return stateLabels_;
}


Foam::label Foam::phasePropertiesList::size() const
{
    return props_.size();
}


const Foam::phaseProperties&
Foam::phasePropertiesList::operator[](const label phaseI) const
{
    return props_[phaseI];
}


// ************************************************************************* //
