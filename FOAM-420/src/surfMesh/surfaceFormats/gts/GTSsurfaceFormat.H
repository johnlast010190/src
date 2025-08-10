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

Class
    Foam::fileFormats::GTSsurfaceFormat

Description
    Provide a means of reading/writing GTS format.
    The output is never sorted by zone and is only written if it consists
    entirely of triangles.

SourceFiles
    GTSsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef GTSsurfaceFormat_H
#define GTSsurfaceFormat_H

#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class GTSsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class GTSsurfaceFormat
:
    public UnsortedMeshedSurface<Face>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        GTSsurfaceFormat(const GTSsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const GTSsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct from file name
        GTSsurfaceFormat(const fileName&);


    // Selectors

        //- Read file and return surface
        static autoPtr<UnsortedMeshedSurface<Face>> New(const fileName& name)
        {
            return autoPtr<UnsortedMeshedSurface<Face>>
            (
                new GTSsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~GTSsurfaceFormat()
    {}


    // Member Functions

        //- Write MeshedSurface
        static void write
        (
            const fileName&,
            const MeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface, the output remains unsorted
        static void write
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFormats/gts/GTSsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
