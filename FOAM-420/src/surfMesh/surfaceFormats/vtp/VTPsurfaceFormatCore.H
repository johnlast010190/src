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

Class
    Foam::fileFormats::VTPsurfaceFormatCore

Description
    Internal class used by the VTPsurfaceFormat

SourceFiles
    VTPsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef VTPsurfaceFormatCore_H
#define VTPsurfaceFormatCore_H

#include "MeshedSurface/MeshedSurface.H"
#include "vtk/format/foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class VTPsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class VTPsurfaceFormatCore
{
protected:

    // Protected Member Functions

        //- Write file header information with points
        static void writeHeader
        (
            vtk::formatter& format,
            const pointField& pts,
            const label nFaces
        );

        //- Write file footer
        static void writeFooter(vtk::formatter& format);


        //- Write regions (zones) information as CellData
        static void writeCellData
        (
            vtk::formatter& format,
            const UList<surfZone>& zones
        );

        //- Write regions (zones) information as CellData
        static void writeCellData
        (
            vtk::formatter& format,
            const labelUList& zoneIds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
