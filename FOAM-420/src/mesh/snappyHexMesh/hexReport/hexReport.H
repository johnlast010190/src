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
    (c) 2014-1014 Esi Ltd.

Class
    Foam::hexReport

Description
    All to do with reporting hex mesh statistics

SourceFiles
    hexReport.C

\*---------------------------------------------------------------------------*/

#ifndef hexReport_H
#define hexReport_H

#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class hexReport Declaration
\*---------------------------------------------------------------------------*/

class hexReport
{

private:

    // Private data

        //- Mesh+surface
        const meshRefinement& meshRefiner_;

        //- initialisation time
        scalar initialisationTime_;

        //- Refinement time
        scalar refinementTime_;

        //- Dual pre-processing time
        scalar dualTime_;

        //- Snapping time
        scalar snappingTime_;

        //- Layer time
        scalar layerTime_;

        //- Number of layer cells
        label numLayerCells_;

        //- Percentile layer coverage
        scalar layerCoverage_;

        //-Layer information per patch
        List<Tuple2<word,scalarField>> patchLayerInfo_;

public:

    //- Runtime type information
    TypeName("hexReport");

    // Constructors

        //- Construct from components
        hexReport
        (
            const meshRefinement& meshRefiner
        );

    // Destructor

    virtual ~hexReport(){};

    // Member Functions

        void setInitialisationTime(const scalar time)
        {
            initialisationTime_ = time;
        }

        void setRefineTime(const scalar time)
        {
            refinementTime_ = time;
        }

        void setDualSetupTime(const scalar time)
        {
            dualTime_ = time;
        }

        void setSnapTime(const scalar time)
        {
            snappingTime_ = time;
        }

        void setLayerTime(const scalar time)
        {
            layerTime_ = time;
        }

        void setNumLayerCells(const label nCells)
        {
            numLayerCells_ = nCells;
        }

        void setLayerCoverage(const scalar coverage)
        {
            layerCoverage_ = coverage;
        }

        void setLayerPatchInfo(const List<Tuple2<word,scalarField>> patchLayerInfo)
        {
            patchLayerInfo_ = patchLayerInfo;
        }

        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
