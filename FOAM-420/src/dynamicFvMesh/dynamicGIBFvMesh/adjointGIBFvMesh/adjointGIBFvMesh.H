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
    (c) 2018 Esi Ltd.

Class
    Foam::adjointGIBFvMesh

Description
    Class for snapping internal point to an STL and making the faceZone needed
    for the GIB

SourceFiles
    adjointGIBFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef adjointGIBFvMesh_H
#define adjointGIBFvMesh_H

#include "dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.H"
#include "dynamicGIBFvMesh/deformingBodyGIBFvMesh/deformingBodyGIBFvMesh.H"
#include "fields/Fields/scalarField/scalarIOField.H"
#include "meshes/primitiveShapes/objectHit/pointHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adjointGIBFvMesh Declaration
\*---------------------------------------------------------------------------*/

class adjointGIBFvMesh
:
    public deformingBodyGIBFvMesh
{

    // Private data

        //- Reference to adjoint properties
        const IOdictionary* adjPropertiesPtr_;

        //- old velocities of the interface
    //    mutable scalarField* pGOld_;

    // Private functions

        //- courant based scaling based on a max and a mean courant
        void courantScaling(scalarField&);

        //- point curvature calculation
        virtual scalarField calcCurvature();

        //- if new points pf the polypatch are moving out of the  basemesh
        //  bounds, find the intersections at the baseMesh boundary and
        //  then snap there (!!! round errors careful)
        virtual void constraintPolyPointsOutMotion
        (
            vectorField&
        );

        virtual void nearBoundaryIntersectionsChecking
        (
            vectorField&
        );

        virtual void checkConcaveBoundaryPatchPoints
        (
            vectorField&
        );

        //-positive point sens at  boundary points
        // which are on the GIB with positive point sensitivity
        // fixed to zero
        virtual void filterBoundaryPoints
        (
            scalarField&,
            const indirectPolyPatch&
        );

        //- mark the boundary faces of the GIB that have not moved and
        //  the sensitivity want to move outwards.
        //  these faces are constrained inside the smoother
        //  COUPLED needs correction
        virtual boolList findConstraintFaces();

        //- fix the points of the constraint patches (if moved)
        //  inlet outlet currently only supported
        virtual void fixConstraintPatches
        (
            pointField&
        );

        //- if processor patch point changed location, then sync the other side
        virtual void syncProcBoundaryPoints
        (
            pointField&,
            const pointField&
        );


        //- Disallow default bitwise copy construct
        adjointGIBFvMesh(const adjointGIBFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const adjointGIBFvMesh&);


protected:

        virtual bool includeWalls() const
        {
            return true;
        }

        //- do tangential Boundary Motion
        virtual bool tangentialBoundaryMotion() const
        {
            return false;
        }

        //- move old GIB based on the curvatture and the sensitivities
        virtual tmp<vectorField> computeNewPoints
        (
            primitivePatch&,
            const scalarField&
        );


public:

    //- Runtime type information
    TypeName("adjointGIBFvMesh");


    // Constructors

        //- Construct from IOobject
        adjointGIBFvMesh(const IOobject& io);


    //- Destructor
    ~adjointGIBFvMesh();


    // Member Functions

        //- faceZone initialization using the stl information
        virtual void updateInit();

        //- update the gib
        virtual bool update();

        //- reference to the adjoint properties
        const IOdictionary& adjProperties();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
