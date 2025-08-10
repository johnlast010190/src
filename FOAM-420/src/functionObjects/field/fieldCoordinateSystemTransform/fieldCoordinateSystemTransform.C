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

\*---------------------------------------------------------------------------*/

#include "fieldCoordinateSystemTransform/fieldCoordinateSystemTransform.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldCoordinateSystemTransform, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fieldCoordinateSystemTransform,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
fieldCoordinateSystemTransform
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(),
    coordSysPtr_
    (
        coordinateSystem::New(mesh_, dict, coordinateSystem::typeName_())
    )
{
    read(dict);

    Info<< type() << " " << name << ":" << nl
        << "   Applying transformation from global Cartesian to local "
        << *coordSysPtr_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
~fieldCoordinateSystemTransform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word
Foam::functionObjects::fieldCoordinateSystemTransform::transformFieldName
(
    const word& fieldName
) const
{
    return fieldName + ":Transformed";
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("fields") >> fieldSet_;

    return true;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::execute()
{
    forAll(fieldSet_, fieldi)
    {
        transform<scalar>(fieldSet_[fieldi]);
        transform<vector>(fieldSet_[fieldi]);
        transform<sphericalTensor>(fieldSet_[fieldi]);
        transform<symmTensor>(fieldSet_[fieldi]);
        transform<tensor>(fieldSet_[fieldi]);
    }

    return true;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::write()
{
    forAll(fieldSet_, fieldi)
    {
        writeObject(transformFieldName(fieldSet_[fieldi]));
    }

    return true;
}


// ************************************************************************* //
