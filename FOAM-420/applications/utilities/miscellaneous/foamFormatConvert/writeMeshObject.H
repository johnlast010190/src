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
    (c) 2011-17 OpenFOAM Foundation

Global
    writeMeshObject

Description
    Write a mesh object in the format specified in controlDict.

\*---------------------------------------------------------------------------*/

#ifndef writeMeshObject_H
#define writeMeshObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CheckType = Type>
inline bool writeMeshObject
(
    const word& name,
    const fileName& meshDir,
    Time& runTime,
    const bool strictTypeChecking = true,
    const bool disableHeaderChecking = false
)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    bool writeOk = false;

    if (io.typeHeaderOk<CheckType>(strictTypeChecking, true, false))
    {
        Info<< "        Reading " << io.headerClassName()
            << " : " << name << endl;

        // Switch off type checking (for reading e.g. faceZones as
        // generic list of dictionaries).
        word oldTypeName;
        if (disableHeaderChecking)
        {
            oldTypeName = Type::typeName;
            const_cast<word&>(Type::typeName) = word::null;
        }

        Type meshObject(io);

        if (disableHeaderChecking)
        {
            const_cast<word&>(Type::typeName) = oldTypeName;
            // Fake type back to what was in field
            const_cast<word&>(meshObject.type()) = io.headerClassName();
        }

        Info<< "        Writing " << name << endl;
        writeOk = meshObject.regIOobject::write();
    }

    return writeOk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
