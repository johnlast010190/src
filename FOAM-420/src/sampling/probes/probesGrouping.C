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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "probes/probes.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/IOobjectList/IOobjectList.H"
#include "primitives/strings/lists/stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::probes::clearFieldGroups()
{
    scalarFields_.clear();
    vectorFields_.clear();
    sphericalTensorFields_.clear();
    symmTensorFields_.clear();
    tensorFields_.clear();

    surfaceScalarFields_.clear();
    surfaceVectorFields_.clear();
    surfaceSphericalTensorFields_.clear();
    surfaceSymmTensorFields_.clear();
    surfaceTensorFields_.clear();
}


Foam::label Foam::probes::classifyFields()
{
    label nFields = 0;
    clearFieldGroups();

    HashTable<wordHashSet> available =
    (
        loadFromFiles_
      ? IOobjectList(obr(), mesh_.time().timeName()).classes(fieldSelection_)
      : obr().classes(fieldSelection_)
    );

    forAllConstIters(available, iter)
    {
        const word& fieldType = iter.key();
        const wordList fieldNames = iter.object().sortedToc();

        const label count = fieldNames.size(); // pre-filtered, so non-empty

        if (fieldType == volScalarField::typeName)
        {
            scalarFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volVectorField::typeName)
        {
            vectorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            sphericalTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            symmTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volTensorField::typeName)
        {
            tensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceScalarField::typeName)
        {
            surfaceScalarFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceVectorField::typeName)
        {
            surfaceVectorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceSphericalTensorField::typeName)
        {
            surfaceSphericalTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceSymmTensorField::typeName)
        {
            surfaceSymmTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceTensorField::typeName)
        {
            surfaceTensorFields_.append(fieldNames);
            nFields += count;
        }
    }

    return nFields;
}


// ************************************************************************* //
