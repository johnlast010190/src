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

\*---------------------------------------------------------------------------*/

#include "fields/cloud/cloud.H"
#include "ensight/type/ensightPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newData
(
    const word& name
) const
{
    autoPtr<ensightFile> output;

    if (Pstream::master())
    {
        const ensight::VarName varName(name);
        output = createDataFile(varName);

        // description
        output().write
        (
            string
            (
                padded(timeIndex_) / varName
              + " <" + pTraits<Type>::typeName + ">"
            )
        );
        output().newline();

        // note variable for later use
        noteVariable(varName, ensightPTraits<Type>::typeName);
    }

    return output;
}


template<class Type>
Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newCloudData
(
    const word& cloudName,
    const word& name
) const
{
    autoPtr<Foam::ensightFile> output;

    if (Pstream::master())
    {
        const ensight::VarName varName(name);
        output = createCloudFile(cloudName, varName);

        // description
        output().write
        (
            string
            (
                padded(timeIndex_) / cloudName / varName
              + " <" + pTraits<Type>::typeName + ">"
            )
        );
        output().newline();

        // note cloud variable for later use
        noteCloud(cloudName, varName, ensightPTraits<Type>::typeName);
    }

    return output;
}


// ************************************************************************* //
