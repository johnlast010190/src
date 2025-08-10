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

Class
    Foam::fileFormats::STARCDsurfaceFormatCore

Description
    Internal class used by the STARCDsurfaceFormat

SourceFiles
    STARCDsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDsurfaceFormatCore_H
#define STARCDsurfaceFormatCore_H

#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "MeshedSurface/MeshedSurface.H"
#include "starcd/STARCDCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class STARCDsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class STARCDsurfaceFormatCore
:
    public STARCDCore
{
protected:

    // Protected Member Functions

    static Map<word> readInpCellTable(IFstream&);

    static void writeCase
    (
        Ostream&,
        const pointField&,
        const label nFaces,
        const UList<surfZone>&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
