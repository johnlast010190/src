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
    (c) 2019 Esi Ltd.


Class
    Foam::npGrid

Group
    grpFoamMap

Description
   A class designed to convert Numpy array into OpenFOAM field for visualization
   in Paraview

SourceFile
    npGrid.H

\*---------------------------------------------------------------------------*/

#pragma once

#include "gridField.H"

namespace Foam
{

class npGrid: public gridField
{
    public:
        std::map<word,nodeField> scalarFields_;
        std::map<word,nodeVectors> vectorFields_;
        std::map<word,nodeField> vectorComps_;

        npGrid();
        static void read_text_file
        (
            word &input,
            const word &fname
        );

        void getScalarField
        (
            const word& fname,
            DynamicList<scalar> &scfield
        );

        void combineVectorField();
        void writeMeshDict();
        void writeScalarField(label k=0);
        void writeVectorField(label k=0);
        void writeVisualization();

        word timeName()
        {
            return "1";
        }

        npGrid
        (
            const word &blkdef,
            const word &npfiles
        );

        virtual ~npGrid(){;}
};

}//Foam
