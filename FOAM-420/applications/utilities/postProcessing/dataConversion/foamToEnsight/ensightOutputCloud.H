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
    (c) 2016 OpenCFD Ltd.

Namespace
    Foam::ensightOutput

Description
    A collection of global functions for writing ensight file content.

SourceFiles
    ensightOutputCloud.C
    ensightOutputCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutputCloud_H
#define ensightOutputCloud_H

#include "ensight/file/ensightFile.H"
#include "ensight/mesh/ensightMesh.H"

#include "memory/autoPtr/autoPtr.H"
#include "db/IOobjects/IOField/IOField.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace ensightCloud
{

//- Write cloud positions
void writePositions
(
    const fvMesh& mesh,
    const word& cloudName,
    const bool exists,
    autoPtr<ensightFile>& output
);


//- Write cloud field, returning true if the field is non-empty.
template<class Type>
bool writeCloudField
(
    const IOField<Type>& field,
    ensightFile& os
);


//- Write cloud field from IOobject, always returning true.
template<class Type>
bool writeCloudField
(
    const IOobject& fieldObject,
    const bool exists,
    autoPtr<ensightFile>& output
);


} // namespace ensightCloud
} // namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
