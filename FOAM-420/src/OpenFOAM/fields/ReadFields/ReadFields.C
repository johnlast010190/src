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
    (c) 2015 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/ReadFields/ReadFields.H"
#include "db/IOobjectList/IOobjectList.H"
#include "db/objectRegistry/objectRegistry.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Read all fields of type. Returns names of fields read. Guarantees all
// processors to read fields in same order.
Foam::wordList Foam::fieldNames
(
    const IOobjectList& fieldObjects,
    const bool syncPar
)
{
    // Get sorted field names. Sorting needed in parallel since different
    // processors (using different file servers) might pick up the files
    // in different order.
    wordList masterNames(fieldObjects.sortedNames());

    if (syncPar && Pstream::parRun())
    {
        // Check that I have the same fields as the master
        const wordList localNames(masterNames);
        Pstream::scatter(masterNames);

        HashSet<word> localNamesSet(localNames);

        forAll(masterNames, i)
        {
            const word& masterFld = masterNames[i];

            HashSet<word>::iterator iter = localNamesSet.find(masterFld);

            if (iter == localNamesSet.end())
            {
                FatalErrorInFunction
                    << "Fields not synchronised across processors." << endl
                    << "Master has fields " << masterNames
                    << "  processor " << Pstream::myProcNo()
                    << " has fields " << localNames << exit(FatalError);
            }
            else
            {
                localNamesSet.erase(iter);
            }
        }

        forAllConstIter(HashSet<word>, localNamesSet, iter)
        {
            FatalErrorInFunction
                << "Fields not synchronised across processors." << endl
                << "Master has fields " << masterNames
                << "  processor " << Pstream::myProcNo()
                << " has fields " << localNames << exit(FatalError);
        }
    }

    return masterNames;
}


// ************************************************************************* //
