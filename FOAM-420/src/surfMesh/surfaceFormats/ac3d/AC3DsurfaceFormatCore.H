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

Class
    Foam::fileFormats::AC3DsurfaceFormatCore

Description
    Internal class used by the AC3DsurfaceFormat

SourceFiles
    AC3DsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef AC3DsurfaceFormatCore_H
#define AC3DsurfaceFormatCore_H

#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "MeshedSurface/MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class AC3DsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class AC3DsurfaceFormatCore
{
protected:

    // Protected Member Functions

        //- Read a type via IStringStream
        template<class Type>
        static Type parse(const string&);

        //- Read cmd, args from IFstream
        static bool readCmd(IFstream&, string& cmd, string& args);

        //- Cue up to cmd, reading args
        static bool cueTo(IFstream&, const string& cmd, string& args);

        //- Cue up to cmd, reading args or exit with a FatalError
        //  returns the command args
        static string cueToOrDie
        (
            IFstream&,
            const string& cmd,
            const string& errorMsg=string::null
        );

        //- Write header with materials
        static void writeHeader(Ostream&, const UList<surfZone>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFormats/ac3d/AC3DsurfaceFormatCoreTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
