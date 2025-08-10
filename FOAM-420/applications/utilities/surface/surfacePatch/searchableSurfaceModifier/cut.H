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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::searchableSurfaceModifier::cut

Description
    Patchify triangles based on orientation w.r.t other (triangulated or
    triangulatable) surfaces

SourceFiles
    cut.C

\*---------------------------------------------------------------------------*/

#ifndef cut_H
#define cut_H

#include "searchableSurfaceModifier.H"
#include "primitives/strings/lists/wordReList.H"
#include "meshes/meshShapes/face/faceList.H"
#include "meshes/primitiveShapes/point/pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class edgeIntersections;

namespace searchableSurfaceModifiers
{

/*---------------------------------------------------------------------------*\
                             Class cut Declaration
\*---------------------------------------------------------------------------*/

class cut
:
    public searchableSurfaceModifier
{
    // Private data

        //- Name of surfaces to cut with
        const wordReList cutterNames_;


   // Private Member Functions

        //- Triangulate faces around 0th vertex
        void triangulate(const faceList&, pointField&, triSurface&) const;

        //- Triangulate searchableSurface (currently only supported for
        //  searchableBox and triSurfaceMesh)
        triSurface& triangulate(const searchableSurface&, triSurface&) const;

        //- Intersect surfaces. Perturb to avoid degenerates.
        bool intersectSurfaces
        (
            triSurface& surf1,
            edgeIntersections& edgeCuts1,
            triSurface& surf2,
            edgeIntersections& edgeCuts2
        ) const;


public:

    //- Runtime type information
    TypeName("cut");


    // Constructors

        //- Construct from dictionary
        cut(const searchableSurfaces&, const dictionary&);

        //- Clone
        autoPtr<searchableSurfaceModifier> clone() const
        {
            NotImplemented;
        }


    //- Destructor
    virtual ~cut();


    // Member Functions

        //- Apply this selector
        virtual bool modify(const labelList& regions, searchableSurface&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace searchableSurfaceModifiers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
