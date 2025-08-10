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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::fileFormats::FLMAsurfaceFormat

Description
    Provide a means of writing AVL/FIRE FLMA format.

SourceFiles
    FLMAsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef FLMAsurfaceFormat_H
#define FLMAsurfaceFormat_H

#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "fire/FIRECore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class FLMAsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FLMAsurfaceFormat
:
    public MeshedSurface<Face>,
    public FIRECore
{
    // Private Member Functions

        static inline void writeShell(OSstream&, const Face&);
        static inline void writeType(OSstream&, const Face&);

        //- Disallow default bitwise copy construct
        FLMAsurfaceFormat(const FLMAsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const FLMAsurfaceFormat<Face>&) = delete;


protected:

    // Protected Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            OSstream&,
            const MeshedSurfaceProxy<Face>&
        );


        //- Write surface mesh components by proxy with/without compression
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&,
            bool compress
        );


public:

    // Constructors

        //- Construct null
        FLMAsurfaceFormat()
        {}


    //- Destructor
    virtual ~FLMAsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write flma file
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }

};


/*---------------------------------------------------------------------------*\
                     Class FLMAZsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FLMAZsurfaceFormat
:
    public FLMAsurfaceFormat<Face>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        FLMAZsurfaceFormat(const FLMAZsurfaceFormat<Face>&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const FLMAZsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct null
        FLMAZsurfaceFormat()
        {}


    //- Destructor
    virtual ~FLMAZsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write flmaz file
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
#include "surfaceFormats/fire/FLMAsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
