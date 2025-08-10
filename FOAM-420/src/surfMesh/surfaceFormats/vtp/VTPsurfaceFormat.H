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
    Foam::fileFormats::VTPsurfaceFormat

Description
    Provide a means of writing VTP (xml) format.
    The output is never sorted by zone.

SourceFiles
    VTPsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTPsurfaceFormat_H
#define VTPsurfaceFormat_H

#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"
#include "surfaceFormats/vtp/VTPsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class VTPsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class VTPsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTPsurfaceFormatCore
{
    // Private Member Functions

        //- Write polygons
        static void writePolys
        (
            vtk::formatter& format,
            const UList<Face>& faces
        );


        //- Disallow default bitwise copy construct
        VTPsurfaceFormat(const VTPsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const VTPsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct null
        VTPsurfaceFormat();


    //- Destructor
    virtual ~VTPsurfaceFormat()
    {}


    // Member Functions

      // Write

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf
        );

        //- Write UnsortedMeshedSurface, the output remains unsorted
        static void write
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf
        );

        //- Write object file
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFormats/vtp/VTPsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
