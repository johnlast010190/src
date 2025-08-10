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
    (c) 2020 Esi Ltd.


Class
    Foam::smBoundBox

Group
    grpfoamMap

Description
    Smart bounding box with irregular boundary. The purpose of this class is
    to efficiently determine whether a given point (x,y,z) is within the
    boundBox. The determination is approximate, to the error of maximum cell
    dimension*ampcoef, where ampcoef  is amplication coefficient.

SourceFile
    smBoundBox.C

\*---------------------------------------------------------------------------*/


#pragma once
#include "cfdTools/general/include/fvCFD.H"
#include "array3d.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class smBoundBox : public boundBox
{
    protected:
        // type of bounding box, simple=0, smart=1
        label boxType_;

        // number of cells in a grid
        Array3d<bool> numCells_;

        // simple bounding box
       // autoPtr<boundBox> bbox_;

        // biggest cell dimension amplication coef
        scalar ampcoef_;

        //grid dimension in x
        label nx_;

        //grid dimension in y
        label ny_;

        //grid dimension in z
        label nz_;

        scalar deltaX_;
        scalar deltaY_;
        scalar deltaZ_;

    public:
        void createBox
        (
            const fvMesh *mesh,
            label ncells=40000
        );

        void defaults();

        inline bool pointInBox (const point &pt) const;

        smBoundBox
        (
            const fvMesh *mesh,
            const word type="boundBox",
            scalar ampcoef=1.5,
            label cellNum=40000
        );

        smBoundBox
        (
            const point &pmin,
            const point &pmax,
            const word type="boundBox",
            scalar ampcoef=1.5
        );


        smBoundBox();

        // get number of cells in each grid
        void getCellList
        (
            const fvMesh *mesh,
            label
        );

        virtual ~smBoundBox(){;}

};



bool smBoundBox::pointInBox (const point &pt) const
{

    bool inBoundBox=contains(pt);
    if (boxType_==0)
    {
        return inBoundBox;
    }

    if (!inBoundBox)
    {
        return false;
    }

    const point &pmin=boundBox::min();

    label ix=label((pt.x()-pmin.x())/deltaX_);
    label iy=label((pt.y()-pmin.y())/deltaY_);
    label iz=label((pt.z()-pmin.z())/deltaZ_);
    if (ix>nx_-1 || iy>ny_-1 || iz>nz_-1)
    {
        return false;
    }
    return numCells_(ix,iy,iz);

}



}//Foam
