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
    (c) 2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermophysicalProperties/thermophysicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalProperties, 0);
    defineRunTimeSelectionTable(thermophysicalProperties,);
    defineRunTimeSelectionTable(thermophysicalProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalProperties::thermophysicalProperties(scalar W)
:
    W_(W)
{}


Foam::thermophysicalProperties::thermophysicalProperties(const dictionary& dict)
:
    W_(readScalar(dict.lookup("W")))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing thermophysicalProperties" << endl;
    }

    const auto ctor = ctorTableLookup("thermophysicalProperties type", ConstructorTable_(), name);
    return autoPtr<thermophysicalProperties>(ctor());
}


Foam::autoPtr<Foam::thermophysicalProperties>
Foam::thermophysicalProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing thermophysicalProperties" << endl;
    }

    const word& thermophysicalPropertiesTypeName = dict.dictName();

    const auto ctor = ctorTableLookup("thermophysicalProperties type", dictionaryConstructorTable_(), thermophysicalPropertiesTypeName);
    return autoPtr<thermophysicalProperties>(ctor(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermophysicalProperties::readIfPresent(const dictionary &dict)
{
    dict.readIfPresent("W", W_);
}


void Foam::thermophysicalProperties::writeData(Ostream& os) const
{
    os  << W_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const thermophysicalProperties& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
