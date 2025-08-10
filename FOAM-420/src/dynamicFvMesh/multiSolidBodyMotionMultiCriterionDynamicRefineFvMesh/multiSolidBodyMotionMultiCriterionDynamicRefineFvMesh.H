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
    (c) 2011-2014 OpenFOAM Foundation
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.
    (c) 2016 Esi Ltd.

Class
    multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh


Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function


SourceFiles
    multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh_H
#define multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh_H

#include "dynamicMultiCriterionRefineFvMesh/dynamicMultiCriterionRefineFvMesh.H"

#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/pointIOField.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiSolidBodyMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh
:
    public dynamicMultiCriterionRefineFvMesh
{
    // Private data

        //- The reference points which are transformed
        pointIOField undisplacedPoints_;

        //- The motion control function
        PtrList<motionSolver> motionPtr_;

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh
        (
            const multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh&
        );

        //- update the pointID list
        void updatePointIDs();


public:

    //- Runtime type information
    TypeName("multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh");


    // Constructors

        //- Construct from IOobject
        multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh(const IOobject& io);


    //- Destructor
    ~multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
