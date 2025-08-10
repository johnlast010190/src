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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionData.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "submodels/Kinematic/PatchInteractionModel/PatchInteractionModel/PatchInteractionModel.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::patchInteractionData::patchInteractionData()
:
    interactionTypeName_("unknownInteractionTypeName"),
    patchName_("unknownPatch"),
    e_(0.0),
    mu_(0.0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::patchInteractionData::interactionTypeName() const
{
    return interactionTypeName_;
}


const Foam::word& Foam::patchInteractionData::patchName() const
{
    return patchName_;
}


Foam::scalar Foam::patchInteractionData::e() const
{
    return e_;
}


Foam::scalar Foam::patchInteractionData::mu() const
{
    return mu_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    patchInteractionData& pid
)
{
    is.check(FUNCTION_NAME);

    const dictionaryEntry entry(dictionary::null, is);

    pid.patchName_ = entry.keyword();
    entry.lookup("type") >> pid.interactionTypeName_;
    pid.e_ = entry.lookupOrDefault<scalar>("e", 1.0);
    pid.mu_ = entry.lookupOrDefault<scalar>("mu", 0.0);

    return is;
}


// ************************************************************************* //
