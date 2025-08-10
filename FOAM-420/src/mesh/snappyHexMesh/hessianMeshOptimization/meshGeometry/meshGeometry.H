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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef meshGeometry_H
#define meshGeometry_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "hessianMeshOptimization/layerCellHandling/cellMetric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

namespace Foam
{

class meshGeometry
{
    // Private data

        //- Description of data_
        const fvMesh& mesh_;

        const pointField& state_;

        vectorField fCtrs_;

        vectorField fAreas_;

        List<List<scalar>> magfAreas_;

        vectorField cellCtrs_;

        scalarField cellVols_;

        scalarField cellSurface_;

        List<tensor> faceCentreDerivatives_;

        List<vector> cellSurfaceDerivatives_;

        List<vector> cellVolumeDerivatives_;

        List<tensor> hessian_;

        labelList pointCellSize_;

        boolList boundaryPoint_;

        labelList processorPoints_;

        boolList boundaryProcessorPoint_;

        labelListList procPointsLists_;

        scalarList cellTypeCoef_;

        void updateAreas();

        void updateVolumes();

        void updateCellSurface();

        void initializeCellSurface();

        void calculateFaceCenterDerivatives();

        void calculateCellDerivatives();

        int getTopologyIndex(const label& cellLabel);

        void numberOfCellsPerPoint();

        void findBoundaryPoints();

        void findBoundaryProcessorPoints();

        void splitProcessorPoints();

        void calculateCellTypes();

protected:

        cellMetric metric_;

        label pLabel_;

        void updateGeometry(const label& pI);

        void calculateDerivatives(const label& pI);

public:

    meshGeometry
    (
        const fvMesh& mesh,
        const pointField& state,
        const dictionary& dict
    );

    ~meshGeometry();

    void limitDisplacement(const label& pI, vector& ds);

    scalar findMaxDisplacement(const label& pI);

    void processorPointsMaxDisplacement
    (
        scalarField& maxDisp,
        const labelList& procPoints
    );

    const scalarField& cellVols() const;

    const scalarField& cellSurface() const;

    const vectorField& cellCtrs() const;

    const vectorField& fCtrs() const;

    const vectorField& fAreas() const;

    const List<tensor>& faceCentreDerivatives() const;

    const List<vector>& cellVolumeDerivatives() const;

    const List<vector>& cellSurfaceDerivatives() const;

    const List<tensor>& surfaceHessian() const;

    const labelList& pointCellsSize() const
    {
        return pointCellSize_;
    }

    const fvMesh& getMesh() const
    {
        return mesh_;
    }

    const pointField& getPoints() const
    {
        return state_;
    }

    const boolList& boundaryPoint() const
    {
        return boundaryPoint_;
    }

    const boolList& boundaryProcessorPoint() const
    {
        return boundaryProcessorPoint_;
    }

    const labelList& processorPoints() const
    {
        return processorPoints_;
    }

    const scalarList& cellTypeCoef() const
    {
        return cellTypeCoef_;
    }

    const cellMetric& metric() const
    {
        return metric_;
    }

    const labelListList& procPointsLists() const
    {
        return procPointsLists_;
    }
};

} /* namespace Foam */

#endif
