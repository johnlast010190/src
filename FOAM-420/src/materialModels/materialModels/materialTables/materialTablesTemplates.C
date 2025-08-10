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
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "materialTables.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::materialTables::updateTable
(
    const HashTable<Type>& table
)
{
    const wordList tableNames = table.toc();
    forAll(tableNames, tableI)
    {
        const Type& modelTable = table[tableNames[tableI]];
        const wordList modelNames = modelTable.toc();
        forAll(modelNames, modelI)
        {
            const word& modelName = modelNames[modelI];

            // Using function cast type instead of name in the table
            modelTable[modelName]->updateTable
            (
                modelTable[modelName]->funcType()
            );
        }
    }
}


template<class Type>
void Foam::materialTables::reportModelsInfo
(
    const HashTable<Type>& table
) const
{
    Info<< endl;
    wordList topLevel = table.toc();
    forAll(topLevel, tableI)
    {
        const word tName(topLevel[tableI]);
        Info<< tName << endl;
        for (unsigned long i=0; i<tName.size(); ++i)
        {
            Info<< "-";
        }
        Info<< endl;
        const wordList modelsToUpdate = table[tName].sortedToc();
        forAll(modelsToUpdate, modelI)
        {
            const word modelName(modelsToUpdate[modelI]);
            Info<< "  " << modelName << " ";
            for (unsigned long i=0; i<(12 - modelName.size()); ++i)
            {
                Info<< " ";
            }
            Info<< "" << table[tName][modelName]->type()
                << ";" << endl;
        }
        Info<< endl;
    }
    Info<< endl;
}


template<class Type>
bool Foam::materialTables::foundModel
(
    const word& tName,
    const word& modelName
)
{
    if (typeid(Type) == typeid(scalar))
    {
        return
            scalarModels_.found(tName)
         && scalarModels_[tName].found(modelName);
    }
    else if (typeid(Type) == typeid(vector))
    {
        return
            vectorModels_.found(tName)
         && vectorModels_[tName].found(modelName);
    }
    else if (typeid(Type) == typeid(tensor))
    {
        return
            tensorModels_.found(tName)
         && tensorModels_[tName].found(modelName);
    }

    return false;
}


template<class objectType>
void Foam::materialTables::addSpeciesAndSpeciesMixtures
(
    const word& phaseName,
    const word& noDictModelType,
    const word& noDictDepGroup
)
{
    addSpeciesAndSpeciesMixtures
    (
        objectType::typeName,
        objectType::listModels(),
        phaseName,
        noDictModelType,
        noDictDepGroup
    );
}


// ************************************************************************* //
