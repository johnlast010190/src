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
    (c) 2017 Esi Ltd.

Class
    surfaceAreaSmoothing

Description
    surface area implicit smoothing

SourceFiles
    surfaceAreaSmoothing.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceAreaSmoothing_H
#define surfaceAreaSmoothing_H

#include "db/typeInfo/typeInfo.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "faMesh/faMesh.H"
#include "include/faCFD.H"
#include "faMesh/geodeticWallDist/geodeticWallDist/geodeticWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class surfaceAreaSmoothing Declaration
\*---------------------------------------------------------------------------*/

class surfaceAreaSmoothing
{

    // Private Member

        const fvMesh& mesh_;

        faMesh areaMesh_;

        const label& patchID_;

        boolList constraintFaces_;

        dimensionedScalar r_;

        scalar standardDeviation_;

        scalar cutOffAngle_;

        //create zero band for smooth transition
        bool smoothTransition_;

        const volScalarField& G_;

        mutable areaScalarField* Gsm_;

        mutable areaScalarField* Gsource_;

        mutable edgeScalarField* DsF_;

        mutable volScalarField* Gsmoothed_;

        mutable geodeticWallDist* y_;


    // Private Member Functions

        labelList labelListPatch(const label&);

        void smooth();

        void createFAMesh();

        void createFAFields();

        labelList findIsolatedAreas();

        void setSmoothingOptions(const dictionary&);

        //- Clear all demand-driven data
        virtual void clearOutData();


public:


    // Constructors

        //- Construct from components
        surfaceAreaSmoothing
        (
            const fvMesh&,
            const dictionary&,
            const label&,
            boolList
        );


    //- Destructor
    virtual ~surfaceAreaSmoothing()
    {
        clearOutData();
    };


    // Member Functions


        //- update
        void update();

        tmp<scalarField> smoothSens() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
