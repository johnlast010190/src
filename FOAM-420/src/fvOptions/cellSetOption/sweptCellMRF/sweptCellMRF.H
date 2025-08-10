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
    Foam::sweptCellMRF

Description
    Automatically creates MRF zone based on a set of given MRF patches.

SourceFiles
    sweptCellMRF.C

\*---------------------------------------------------------------------------*/


#ifndef sweptCellMRF_H
#define sweptCellMRF_H

#include "cfdTools/general/include/fvCFD.H"
#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sweptCellMRF Declaration
\*---------------------------------------------------------------------------*/

class sweptCellMRF
{

    // Private Data

        //- MRF dictionary
        const dictionary& MRFdict_;

        //- fvMesh
        const fvMesh& mesh_;

        //- List of MRF patches
        wordList MRFPatchNames_;

        //- labelList of MRF patches
        labelList MRFPatchLabels_;

        //- Convection Field
        volScalarField tracerField_;

        //- Coordinate frame pointer
        coordinateFrame* coorFramePtr_;

        //- Axes of rotation
        vector axisOfRotation_;

        //- Center of rotation
        vector CofR_;

        //- Tracer threshold value
        scalar tracerThreshold_;

        bool writeTracerField_;

        label numberOfIterations_;
    // Private Functions
        void calculateConvectionField();

        void setFvSolutionDict();

        void setFvSchemesDict();
public:

    // Constructors

        //- Construct from componets
        sweptCellMRF
        (
            const dictionary& dict_,
            const fvMesh& mesh_
        );

        //- Destructor
        ~sweptCellMRF();

        // Member Functions
        labelList getMRFCells();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
