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
    (c) 2011-2012 OpenFOAM Foundation

InClass
    cfMeshOptimize

\*---------------------------------------------------------------------------*/

#include "autoOptimize/cfMeshOptimize/cfMeshOptimize.H"
#include "utilities/smoothers/geometry/meshOptimizer/meshOptimizer.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "global/unitConversion/unitConversion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cfMeshOptimize, 0);
    addToRunTimeSelectionTable
    (
        autoOptimize,
        cfMeshOptimize,
        dictionary
    );
}


void Foam::cfMeshOptimize::movePoints(Foam::pointField& newPoints)
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const Time& runTime = mesh_.time();

    polyMeshGen pmg
    (
        runTime,
        newPoints,
        mesh_.faces(),
        mesh_.cells(),
        bMesh
    );

    //- construct the smoother
    meshOptimizer mOpt(pmg);

    //- clear geometry information before volume smoothing
    pmg.clearAddressingData();

    //- perform optimisation using the laplace smoother and
    mOpt.optimizeMeshFV
    (
        numLaplaceIterations_,
        maxNumGlobalIterations_,
        maxNumIterations_,
        maxNumSurfaceIterations_,
        errorBufferLayers_,
        relaxedCheck_,
        relaxedBoundaryCheck_,
        checkWarped_,
        minFaceArea_
    );

    //- perform optimisation of worst quality faces
    mOpt.optimizeMeshFVBestQuality(maxNumGlobalIterations_, qualityThreshold_);

    //- perform orthogonality optimisation
    if (maxOrth_ < 180 - SMALL)
    {
        scalar orthoThreshold  = Foam::cos(degToRad(maxOrth_));

        mOpt.optimizeMeshFVOrthogonality
        (
            maxNumGlobalIterations_,
            orthoThreshold
        );
    }

    //- check the mesh again and untangle bad regions if any of them exist
    if (qualityThreshold_ > SMALL)
    {
        mOpt.untangleMeshFV
        (
            maxNumGlobalIterations_,
            maxNumIterations_,
            maxNumSurfaceIterations_,
            errorBufferLayers_,
            relaxedCheck_,
            relaxedBoundaryCheck_,
            checkWarped_,
            minFaceArea_
        );
    }

    if (syncPts_)
    {
        vectorField dispVec( pmg.points() - newPoints );
        syncTools::syncPointList
        (
            mesh_,
            dispVec,
            maxMagSqrEqOp<point>(),
            vector::zero
        );
        newPoints += dispVec;
    }
    else
    {
        newPoints = pmg.points();
    }
}


void Foam::cfMeshOptimize::optimize()
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const Time& runTime = mesh_.time();

    polyMeshGen pmg
    (
        runTime,
        mesh_.points(),
        mesh_.faces(),
        mesh_.cells(),
        bMesh
    );

    //- construct the smoother
    meshOptimizer mOpt(pmg);

    //- clear geometry information before volume smoothing
    pmg.clearAddressingData();

    //- perform optimisation using the laplace smoother and
    mOpt.optimizeMeshFV
    (
        numLaplaceIterations_,
        maxNumGlobalIterations_,
        maxNumIterations_,
        maxNumSurfaceIterations_,
        errorBufferLayers_,
        relaxedCheck_,
        relaxedBoundaryCheck_,
        checkWarped_,
        minFaceArea_
    );

    //- perform optimisation of worst quality faces
    mOpt.optimizeMeshFVBestQuality(maxNumGlobalIterations_, qualityThreshold_);

    //- perform orthogonality optimisation
    if (maxOrth_ < 180 - SMALL)
    {
        scalar orthoThreshold  = Foam::cos(degToRad(maxOrth_));

        mOpt.optimizeMeshFVOrthogonality
        (
            maxNumGlobalIterations_,
            orthoThreshold
        );
    }

    //- check the mesh again and untangl bad regions if any of them exist
    if (qualityThreshold_ > SMALL)
    {
        mOpt.untangleMeshFV
        (
            maxNumGlobalIterations_,
            maxNumIterations_,
            maxNumSurfaceIterations_,
            errorBufferLayers_,
            relaxedCheck_,
            relaxedBoundaryCheck_,
            checkWarped_,
            minFaceArea_
        );
    }

    if (syncPts_)
    {
        //sync the displacement
        vectorField updatedPts( pmg.points() - mesh_.points() );
        syncTools::syncPointList
        (
            mesh_,
            updatedPts,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        updatedPts += mesh_.points();
        mesh_.movePoints(updatedPts);
    }
    else
    {
        mesh_.movePoints(pmg.points());
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cfMeshOptimize::cfMeshOptimize
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    autoOptimize(mesh, dict),
    mesh_(mesh)
{
    if (dict.found("cfMeshOptimizeCoeffs"))
    {
        dictionary coeffsDict = dict.subDict("cfMeshOptimizeCoeffs");

        qualityThreshold_ =
            coeffsDict.lookupOrDefault<scalar>("qualityThreshold", 0.0);
        numLaplaceIterations_ =
            coeffsDict.lookupOrDefault<label>("numLaplaceIterations", 0);
        maxNumGlobalIterations_ =
            coeffsDict.lookupOrDefault<label>("maxNumGlobalIterations", 10);
        maxNumIterations_ =
            coeffsDict.lookupOrDefault<label>("maxNumIterations", 25);
        maxNumSurfaceIterations_ =
            coeffsDict.lookupOrDefault<label>("maxNumSurfaceIterations", 2);
        errorBufferLayers_ =
            coeffsDict.lookupOrDefault<label>("errorBufferLayers", 5);
        relaxedCheck_ =
            coeffsDict.lookupOrDefault<Switch>("relaxedCheck", false);
        relaxedBoundaryCheck_ =
            coeffsDict.lookupOrDefault<Switch>("relaxedBoundaryCheck", false);
        checkWarped_ =
            coeffsDict.lookupOrDefault<Switch>("checkWarped", false);
        syncPts_ =
            coeffsDict.lookupOrDefault<Switch>("syncPts", true);
        maxOrth_ =
            coeffsDict.lookupOrDefault<scalar>("maxNonOrtho", 180.);
        minFaceArea_ =
           coeffsDict.lookupOrDefault<scalar>("minFaceArea", VSMALL);
    }
    else
    {
        qualityThreshold_ = 0.0;
        numLaplaceIterations_ = 0;
        maxNumGlobalIterations_ = 10;
        maxNumIterations_ = 25;
        maxNumSurfaceIterations_ = 2;
        errorBufferLayers_ = 5;
        relaxedCheck_ = false;
        relaxedBoundaryCheck_ = false;
        checkWarped_ = false;
        syncPts_ = true;
        maxOrth_ = 180.;
        minFaceArea_ = VSMALL;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cfMeshOptimize::~cfMeshOptimize()
{}

// ************************************************************************* //
