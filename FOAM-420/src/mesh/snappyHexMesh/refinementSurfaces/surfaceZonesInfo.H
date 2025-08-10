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
    (c) 2015 OpenCFD Ltd.
    (c) 2014 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::surfaceZonesInfo

Description

SourceFiles
    surfaceZonesInfo.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceZonesInfo_H
#define surfaceZonesInfo_H

#include "primitives/enums/NamedEnum.H"
#include "meshes/primitiveShapes/point/point.H"
#include "primitives/strings/word/word.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "primitives/ints/lists/labelList.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;
class searchableSurfaces;
class polyMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                      Class surfaceZonesInfo Declaration
\*---------------------------------------------------------------------------*/

class surfaceZonesInfo
{
public:

    //- Types of selection of area
    enum areaSelectionAlgo
    {
        INSIDE,
        OUTSIDE,
        INSIDEPOINT,
        NONE
    };

    static const NamedEnum<areaSelectionAlgo, 4> areaSelectionAlgoNames;

    //- What to do with faceZone faces
    enum faceZoneType
    {
        INTERNAL,
        BAFFLE,
        BOUNDARY
    };

    static const NamedEnum<faceZoneType, 3> faceZoneTypeNames;


private:

    // Private data

        //- Per 'interface' surface : name of faceZone to put faces into
        word faceZoneName_;

        //- Per 'interface' surface : name of cellZone to put cells into
        word cellZoneName_;

        //- Per 'interface' surface : (only used if surface is closed)
        //  How to select zone cells : surface inside or outside or given
        //  inside location.
        areaSelectionAlgo zoneInside_;

        //- If zoneInside=location gives the corresponding inside point
        point zoneInsidePoint_;

        //- Per 'interface' surface :
        //  What to do with outside
        faceZoneType faceType_;

        //- Per 'interface' surface :
        //  How to handle free standing zone faces
        bool freeStanding_;

        //- Per 'interface' surface :
        //  Use walking zoneing method for surface
        bool zoneByWalk_;

        //- Per 'interface' surface :
        //  Whether to treat faces as internal for error checks
        bool baffleCheck_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const surfaceZonesInfo&);


public:

    // Constructors

        //- Construct from surfaces and dictionary
        surfaceZonesInfo
        (
            const searchableSurface& surface,
            const dictionary& surfacesDict,
            const bool globalFreeStanding,
            const bool globalZoneWalk
        );

        //- Construct from components
        surfaceZonesInfo
        (
            const word& faceZoneNames,
            const word& cellZoneNames,
            const areaSelectionAlgo& zoneInside,
            const point& zoneInsidePoints,
            const faceZoneType& faceType
        );

        //- Copy constructor
        surfaceZonesInfo(const surfaceZonesInfo&);

        //- Return clone
        autoPtr<surfaceZonesInfo> clone() const
        {
            return autoPtr<surfaceZonesInfo>(new surfaceZonesInfo(*this));
        }


    // Member Functions

        // Access

            //- Per 'interface' surface : empty or name of faceZone to put
            //  faces into
            const word& faceZoneName() const
            {
                return faceZoneName_;
            }

            //- Per 'interface' surface : empty or name of cellZone to put
            //  cells into
            const word& cellZoneName() const
            {
                return cellZoneName_;
            }

            const areaSelectionAlgo& zoneInside() const
            {
                return zoneInside_;
            }

            //- Get specified inside locations for surfaces with a cellZone
            const point& zoneInsidePoint() const
            {
                return zoneInsidePoint_;
            }

            //- How to handle face of surfaces with a faceZone
            const faceZoneType& faceType() const
            {
                return faceType_;
            }

            //- Whether to allow free-standing zone face
            const bool& freeStanding() const
            {
                return freeStanding_;
            }

            //- Whether to allow zoneing walk methood
            const bool& zoneByWalk() const
            {
                return zoneByWalk_;
            }

            //- Whether to treat faces as internal for error checks
            const bool& baffleCheck() const
            {
                return baffleCheck_;
            }


        // Query

            //- Get indices of baffle and internal surface
            static labelList getBoundaryNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of baffle and internal surface
            static labelList getNonBoundaryNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of unnamed surfaces (surfaces without faceZoneName)
            static labelList getUnnamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of unnamed surfaces (surfaces without faceZoneName) and
            // boundary named surfaces
            static labelList getUnnamedAndBoundaryNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of named surfaces (surfaces with faceZoneName)
            static labelList getNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of named non-boundary free-standing surfaces
            static labelList getNonBoundaryFreeStandingSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of named surfaces without a cellZone
            static labelList getStandaloneNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            //- Get indices of surfaces with a cellZone that are closed and
            //  have 'inside' or 'outside' selection.
            static labelList getClosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that are unclosed
            static labelList getUnclosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that are closed.
            static labelList getAllClosedNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const searchableSurfaces& allGeometry,
                const labelList& surfaces
            );

            //- Get indices of surfaces with a cellZone that have 'insidePoint'
            //  section.
            static labelList getInsidePointNamedSurfaces
            (
                const PtrList<surfaceZonesInfo>& surfList
            );

            static label addCellZone
            (
                const word& name,
                const labelList& addressing,
                polyMesh& mesh
            );

            static labelList addCellZonesToMesh
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const labelList& namedSurfaces,
                polyMesh& mesh
            );

            static label addFaceZone
            (
                const word& name,
                const labelList& addressing,
                const boolList& flipMap,
                polyMesh& mesh
            );

            static labelList addFaceZonesToMesh
            (
                const PtrList<surfaceZonesInfo>& surfList,
                const labelList& namedSurfaces,
                polyMesh& mesh
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
