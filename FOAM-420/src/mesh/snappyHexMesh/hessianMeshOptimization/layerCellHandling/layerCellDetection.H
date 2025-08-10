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
    Foam::layerCellDetection

Description
    Detects the layers cells of a mesh based on the aspect ratio.
    Creates vertical stacks of layer cells.

SourceFiles
    layerCellDetection.C

\*---------------------------------------------------------------------------*/

#ifndef layerCellDetection_H
#define layerCellDetection_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class layerCellDetection Declaration
\*---------------------------------------------------------------------------*/

class layerCellDetection
{
        const fvMesh& mesh_;

        const scalar& aspectRatioTrigger_;

        const polyBoundaryMesh& patches_;

        scalarField cellMaxDistance_;

        scalarField cellMinDistance_;

        volScalarField aspectRatio_;

        volScalarField layerCells_;

        vectorField layerDirection_;

        List<bool> firstRowLayerCell_;

        List<bool> isALayerCell_;

        List<bool> isAnActualLayerCell_;

        List<vector> layerGrowningDirection_;

        labelList numLayers_;

        scalarField rowOfLayerCells_;

        label layerCounter_;

        List<bool> isALayerProcessorBoundaryCell_;

        List<bool> isALayerProcessorBoundaryFace_;

        vectorField layerDirectionProcessorFaces_;

        labelField numberOflayersProcFace_;
        labelField layerStackCounter_;

        void calculateNumberOfLayers();
        void findCellExtremesDistances();
        void aspectRatioCalc();
        void findLayerCells();
        void calculateGrowingDirection();
        void approximateGrowingDirection();
        pointField calcSmoothFaceNormals();
        autoPtr<indirectPrimitivePatch> getLayerPatch();
        void findFirstRowLayerCells();
        void buildStack(const label& cI, const label& pI);
        void buildStackFromProc(const label& cI, const label& layerNu);

        void handlePentahedraCells();
        bool isAPyramidCell(const label& cell);
        void excludeCell(const label& cell);
        void excludePyramidCells();
        bool isABoundaryFace(const label& face);
        bool isAPrismCell(const label& cellLabel);
        void findProcessorLayerCells();

        void clearData();

        //- Disallow default bitwise copy construct
        layerCellDetection(const layerCellDetection&);

        //- Disallow default bitwise assignment
        void operator=(const layerCellDetection&);

public:

        TypeName("layerCellDetection");

        layerCellDetection
        (
            const fvMesh& mesh,
            const scalar& aspectRatioTrigger
        );

        static autoPtr<layerCellDetection> New();

        virtual ~layerCellDetection();

        const volScalarField& aspectRatio() const;

        const bool& isALayerCell(const label& cell) const;

        const vectorField& layerDirection() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
