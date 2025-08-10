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
    Foam::fileFormats::OFFsurfaceFormat

Description
    Provide a means of reading/writing Geomview OFF polyList format.


See also
    The <a href="http://www.geoview.org">Geoview</a>
    file format information:
    http://www.geomview.org/docs/html/OFF.html#OFF

Note
    When reading, the optional \a colorspec is ignored.
    When writing, it is set to the zone number (integer).

SourceFiles
    OFFsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef OFFsurfaceFormat_H
#define OFFsurfaceFormat_H

#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class OFFsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class OFFsurfaceFormat
:
    public MeshedSurface<Face>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        OFFsurfaceFormat(const OFFsurfaceFormat&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const OFFsurfaceFormat&) = delete;


public:

    // Constructors

        //- Construct from file name
        OFFsurfaceFormat(const fileName&);


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face>> New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face>>
            (
                new OFFsurfaceFormat(name)
            );
        }


    //- Destructor
    virtual ~OFFsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
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
    #include "surfaceFormats/off/OFFsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
