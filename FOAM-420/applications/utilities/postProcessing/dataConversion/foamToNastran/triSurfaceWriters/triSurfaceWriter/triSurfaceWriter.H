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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.

Class
    triSurfaceWriter

Description

SourceFiles
    triSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceWriter_H
#define triSurfaceWriter_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/autoPtr/autoPtr.H"
#include "triSurface/fields/triSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fileName;

/*---------------------------------------------------------------------------*\
                           Class triSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceWriter
{

public:

    //- Runtime type information
    TypeName("triSurfaceWriter");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            triSurfaceWriter,
            word,
            (),
            ()
        );


    // Selectors

        //- Return a reference to the selected triSurfaceWriter
        static autoPtr<triSurfaceWriter> New(const word& writeFormat);


    // Constructors

        //- Construct null
        triSurfaceWriter();

    // Destructor

        virtual ~triSurfaceWriter();


    // Member Functions

        //- Writes surface + data to file.
        virtual void write
        (
            const fileName& path,           // <root>/<case>/sampleSurfaces
            const fileName& timeDir,        // time name
            const fileName& surfaceName,    // name of surface
            const triSurface& surface,
            PtrList<triSurfaceScalarField>& scalarFields,
            PtrList<triSurfaceVectorField>& vectorFields,
            PtrList<triSurfaceTensorField>& tensorFields
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
