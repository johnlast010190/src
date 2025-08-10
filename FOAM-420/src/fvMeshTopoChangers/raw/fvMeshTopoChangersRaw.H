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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::fvMeshTopoChangers::raw

Description
    Prototype fvMeshTopoChanger without any added functionality.

    Can be used as the starting point for any polyTopoChanger based
    topo-changer.

    Note: run without FOAM_SETNAN. Temporary has unitialised patch
    data when faces are created out of nothing.

SourceFiles
    fvMeshTopoChangersRaw.C
    fvMeshTopoChangersRawTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersRaw_H
#define fvMeshTopoChangersRaw_H

#include "fvMesh/fvMeshTopoChangers/fvMeshTopoChanger/fvMeshTopoChanger.H"
#include "polyTopoChange/polyTopoChanger/polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                   Class fvMeshTopoChangers::raw Declaration
\*---------------------------------------------------------------------------*/

class raw
:
    public fvMeshTopoChanger
{
    // Private Data

        polyTopoChanger topoChanger_;


    // Private Member Functions

        //- Set unmapped values
        template<class Type, template<class> class PatchField, class GeoMesh>
        static void setUnmappedValues
        (
            GeometricField<Type, PatchField, GeoMesh>& fld,
            const PackedBoolList& mappedFace,
            const GeometricField<Type, PatchField, GeoMesh>& baseFld
        );

        template<class Type, template<class> class PatchField, class GeoMesh>
        void zeroUnmappedValues(const PackedBoolList&) const;


public:

    //- Runtime type information
    TypeName("raw");


    // Constructors

        //- Construct from fvMesh
        explicit raw(fvMesh& io);

        //- Disallow default bitwise copy construction
        raw(const raw&) = delete;


    //- Destructor
    virtual ~raw();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const raw&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "raw/fvMeshTopoChangersRawTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
