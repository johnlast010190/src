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
    (c) 2010-2012 Esi Ltd.
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::shellSurfaces

Description
    Encapsulates queries for volume refinement ('refine all cells within
    shell').

SourceFiles
    shellSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef shellSurfaces_H
#define shellSurfaces_H

#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "primitives/Tuple2/Tuple2.H"
#include "fvMesh/fvMesh.H"
#include "primitives/Vector/labelVector/labelVector.H"
#include "coordinate/systems/coordinateSystem.H"
#include "meshControl/meshControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;

/*---------------------------------------------------------------------------*\
                           Class shellSurfaces Declaration
\*---------------------------------------------------------------------------*/

class shellSurfaces
{
public:

    // Public data types

        //- Volume refinement controls
        enum refineMode
        {
            INSIDE,         // Refine all inside shell
            OUTSIDE,        // ,,         outside
            DISTANCE,       // Refine based on distance to shell
            DIRECTION       // Refine based on directional distances
        };

        //- Volume refinement type
        enum refineType
        {
            ISO,         // isotropic
            ANISO        // an-isotropic refinement
        };


private:

    // Private data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are shells
        labelList shells_;

        //- Per shell whether to refine inside or outside
        List<refineMode> modes_;

        //- Per shell whether isotropic or anisotropic refinement
        List<refineType> types_;

        //- Per shell the list of ranges (stored as Tuple2 to
        // cover directional refinement)
        List<List<Tuple2<point,point>>> distances_;

        //- Per shell per distance the refinement level
        List<List<labelVector>> levels_;

        //- Per shell whether to refine to the surface (anisotrpic refinement)
        List<bool> refineSurface_;

        //- Per shell whether to perform more expensive inside/outside
        // check on interface cells
        List<bool> additionalInsideCheck_;

        // Optional co-ordinate system used for directional refinement
        List<autoPtr<coordinateSystem>> coord_;

        // Gap level refinement

            //- Per shell, per region the small-gap level specification
            List<List<FixedList<label, 3>>> extendedGapLevel_;

            //- Per shell, per region the small-gap level specification
            List<List<volumeType>> extendedGapMode_;


    // Private data

        //- refineMode names
        static const NamedEnum<refineMode, 4> refineModeNames_;

        //- refineType names
        static const NamedEnum<refineType, 2> refineTypeNames_;


    // Private Member Functions

        //Read coordinate rotation dictonary or set to default
        dictionary getCoordinateDict
        (
            const dictionary& dict
        );

        //- Helper function for initialisation of levels (direction)
        void setAndCheckLevels
        (
            const label shellI,
            const List<Tuple2<Tuple2<point,point>,label>> &
        );

        //- Helper function for initialisation of levels
        template<class T> void setAndCheckLevels
        (
            const label shellI,
            const List<Tuple2<scalar, T>>&
        );

        //- Helper function for checking of gap information
        void checkGapLevels
        (
            const dictionary&,
            const label shellI,
            const List<FixedList<label, 3>>& levels
        );

        //- Find first shell with a level higher than maxLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& maxLevel,
            const bool threaded = false
        ) const;

        //- Find first shell with a level higher than maxLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label shellI,
            const labelList& regions,
            labelList& maxLevel,
            labelList& minCmptLevel,
            labelList& shellCheck,
            const refineType,
            const direction cmpt
        ) const;

        void orient();

        //- Update highest min gap level
        void findHigherGapLevel
        (
            const pointField& pt,
            const labelList& ptLevel,
            const label shellI,
            labelList& gapShell,
            List<FixedList<label, 3>>& gapInfo,
            List<volumeType>& gapMode
        ) const;

        //- Find first shell with a level lower or equal to minLevel. Update
        //  minLevel and shell.
        void findLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& minLevel,
            labelList& shell,
            const refineType=ISO,
            const direction cmpt=0
        ) const;

public:

    // Constructors

        //- Construct from geometry and dictionary
        shellSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary& refineDict,
            const meshControl& controller,
            bool globalAddedCheck = false
        );


    // Member Functions

        // Access

            //- Indices of surfaces that are shells
            const labelList& shells() const
            {
                return shells_;
            }

            //- Per shell per distance the refinement level
            const List<List<labelVector>>& levels() const
            {
                return levels_;
            }

            //- Per shell the refinement mode
            const List<refineMode>& modes() const
            {
                return modes_;
            }

        // Query

            //- Convert distance Tuple field as isotropic distance scalar field
            const scalarField shellIsoDistances(label shelli) const
            {
                const List<Tuple2<point,point>>& dist = distances_[shelli];
                scalarField isodist(dist.size());
                forAll(dist, disti)
                {
                    isodist[disti] = dist[disti].first().x();
                }
                return isodist;
            }

            //- Whether to refine down to the surface
            const boolList& refineSurface() const
            {
                return refineSurface_;
            }

            //- Whether to perform additional inside check
            const boolList& additionalInsideCheck() const
            {
                return additionalInsideCheck_;
            }

            //- Maximum level excluding distance based levels
            label maxInsideOutsideLevel() const;

            //- Highest isotropic refinement shell level
            label maxIsoLevel() const;

            //- Highest anisotropic refinement shell level
            label maxAnisoLevel(direction cmpt) const;

            bool hasAnisotropicRefinement() const;

            //- Highest shell gap level
            labelList maxGapLevel() const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel,
                const bool threaded = false
            ) const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const fvMesh& mesh,
                const pointField& pt,
                const labelList& ptLevel,
                const labelList& ptMap,
                labelList& maxLevel,
                labelList& minCmptLevel,
                labelList& shellCheck,
                const refineType=ISO,
                const direction cmpt=0
            ) const;

            //- Find a shell whose minimum gap level is >= ptLevel
            void findHigherGapLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& gapShell,
                List<FixedList<label, 3>>& gapInfo,
                List<volumeType>& gapMode
            ) const;

            //- Find a shell whose minimum gap level is >= ptLevel. gapInfo
            //  is (0 0 0) if no shell found
            void findHigherGapLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                List<FixedList<label, 3>>& gapInfo,
                List<volumeType>& gapMode
            ) const;

            //- Find first shell (or -1) with level equal or lower than ptLevel.
            void findLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& shell
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
