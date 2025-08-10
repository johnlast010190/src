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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "functionObjects/fieldSelections/solverFieldSelection/solverFieldSelection.H"
#include "fvMesh/fvMesh.H"
#include "volMesh/volMesh.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solverFieldSelection::solverFieldSelection
(
    const objectRegistry& obr,
    const bool includeComponents
)
:
    fieldSelection(obr, includeComponents)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "Registry must be of type " << fvMesh::typeName
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::solverFieldSelection::updateSelection()
{
    List<fieldInfo> oldSet(std::move(selection_));

    DynamicList<fieldInfo> newSelection(oldSet.size());

    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const dictionary& solverDict = mesh.solverPerformanceDict();

    const wordList solvedFieldNames(solverDict.sortedToc());

    for (const fieldInfo& fi : *this)
    {
        for (const word& solvedField : solvedFieldNames)
        {
            if (fi.name().match(solvedField))
            {
                newSelection.append
                (
                    fieldInfo(wordRe(solvedField), fi.component())
                );
                fi.found() = true;
            }
        }
    }

    selection_.transfer(newSelection);

    if (!fieldSelection::checkSelection())
    {
        WarningInFunction
            << "Valid solver fields are: " << solvedFieldNames;
    }

    return selection_ != oldSet;
}


// ************************************************************************* //
