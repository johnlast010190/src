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
    (c) 2011 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Class
    Foam::fileFormats::STARCDMeshWriter

Description
    Writes polyMesh in pro-STAR (v4) bnd/cel/vrt format

    The cellTableId and cellTable information are used (if available).
    Otherwise the cellZones are used (if available).

SourceFiles
    STARCDMeshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDMeshWriter_H
#define STARCDMeshWriter_H

#include "common/writer/meshWriter.H"
#include "starcd/STARCDCore.H"
#include "db/IOstreams/IOstreams/IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::STARCDMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class STARCDMeshWriter
:
    public meshWriter,
    protected STARCDCore
{
    // Private data

        //- Write boundary (bnd) file - default true
        bool writeBoundary_;


    // Private Member Functions

        //- Write cells
        void writeCells(const fileName& baseName) const;

        //- Write boundary
        void writeBoundary(const fileName& baseName) const;

        void getCellTable();

        label findDefaultBoundary() const;


        //- Disallow default bitwise copy construct
        STARCDMeshWriter(const STARCDMeshWriter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const STARCDMeshWriter&) = delete;


public:

    // Constructors

        //- Write mesh files in PROSTAR format
        STARCDMeshWriter
        (
            const polyMesh&,
            const scalar scaleFactor = 1.0,
            const bool writeBndFile = true
        );


    //- Destructor
    virtual ~STARCDMeshWriter();


    // Member Functions

    // Write

        //- Write volume mesh
        virtual bool write
        (
            const fileName& meshName = fileName::null
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
