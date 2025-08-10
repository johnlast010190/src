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
    (c) 2015 Esi Ltd.

Class
    Foam::deformingBodyGIBFvMesh

Description
    Class for snapping internal point to an STL and making the faceZone needed
    for the GIB

SourceFiles
    deformingBodyGIBFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef deformingBodyGIBFvMesh_H
#define deformingBodyGIBFvMesh_H

#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "fields/volFields/volFields.H"
#include "dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.H"
#include "dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIB.H"
#include "dynamicGIBFvMesh/movingGIBTools/constraintBoundaryPoints/constraintBoundaryPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class deformingBodyGIBFvMesh Declaration
\*---------------------------------------------------------------------------*/

class deformingBodyGIBFvMesh
:
    public dynamicGIBFvMesh
{


        virtual void initialization();

        virtual void makeRecAllPoints0() const;

protected:

    // Private data

        //- stl interface stl location
        fileName postProFolder_;

        //- interface output switch
        Switch writeInterface_;

        //- interface properties
        autoPtr<deformingBodyMotionFunction> DBMFPtr_;

    // Protected Member Functions

        virtual bool includeWalls() const
        {
            return false;
        }

        //- do tangential Boundary Motion
        virtual bool tangentialBoundaryMotion() const
        {
            return true;
        }


        //- move old GIB
        virtual tmp<vectorField> movePolyPatch
        (
            primitivePatch&,
            const vectorField&
        );

        //- postPro folder
        void writeInterface(const indirectPolyPatch&);

        //- postPro folder
        void postPro();

        //- Disallow default bitwise copy construct
        deformingBodyGIBFvMesh(const deformingBodyGIBFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const deformingBodyGIBFvMesh&);


public:

    //- Runtime type information
    TypeName("deformingBodyGIBFvMesh");


    // Constructors

        //- Construct from IOobject
        deformingBodyGIBFvMesh(const IOobject&, const word&);

        //- Construct from IOobject
        deformingBodyGIBFvMesh(const IOobject&);

        //- Construct from IOobject
        deformingBodyGIBFvMesh(const IOobject&, const dictionary);

    //- Destructor
    ~deformingBodyGIBFvMesh();


    // Member Functions

        virtual void updateInit(const word&);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- do tangential Boundary Motion
        virtual bool tangentialBoundaryMotion()
        {
            return true;
        }

        //- Correct boundaryConditions
        virtual tmp<vectorField> velocityCorrect(const vectorField&) const;

        //- Calculate the old boundary location
        virtual tmp<vectorField> oldBoundaryLocation() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
