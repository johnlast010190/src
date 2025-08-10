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


\*---------------------------------------------------------------------------*/

#include "smBoundBox.H"

using namespace Foam;

smBoundBox::smBoundBox()
:
boundBox()
{
   defaults();
}

void smBoundBox::defaults()
{

    ampcoef_ = 1.5;
    boxType_ = 0;
    nx_      = 100;
    ny_      = 100;
    nz_      = 1;
}


void smBoundBox::getCellList
(
    const fvMesh *mesh,
    label ncells
)
{
    const point& pmin = boundBox::min();
    const point& pmax = boundBox::max();

    Info<< "smBoundBox min" << pmin << "; smBoundBox max" << pmax << endl;
    Info<< "ncells: " << ncells << endl;

    scalar deltaX = pmax.x() - pmin.x();
    scalar deltaY = pmax.y() - pmin.y();
    scalar deltaZ = pmax.z() - pmin.z();

    Info<< "deltaXYZ: " << deltaX << " " << deltaY << " " << deltaZ << endl;

    scalar dmax = Foam::max(deltaX, deltaY);
    dmax = Foam::max(dmax, deltaZ);

    scalar ax = deltaX/(dmax + SMALL);
    scalar ay = deltaY/(dmax + SMALL);
    scalar az = deltaZ/(dmax + SMALL);

    nx_ = 1;
    ny_ = 1;
    nz_ = 1;

    scalar ann = sqrt(scalar(ncells));

    if (ax*ann < 1.0 && ay*az > 0.0)
    {
        // 2d in x direction
        nz_ = label(sqrt(az*ncells/ay) + 0.5);
        ny_ = label((ay/az)*nz_ + 0.5);
        nx_ = 1;
    }
    else if (ay*ann < 1.0 && ax*az > 0.0)
    {
        // 2d in y direction
        nz_ = label(sqrt(az*ncells/ax) + 0.5);
        nx_ = label((ax/az)*nz_ + 0.5);
        ny_ = 1;
    }
    else if (az*ann < 1.0 && ax*ay > 0.0)
    {
        // 2d in z direction
        ny_ = label(sqrt(ay*ncells/ax) + 0.5);
        nx_ = label((ax/ay)*ny_ + 0.5);
        nz_ = 1;
    }
    else if (ax*ay*az > 0.0)
    {
        scalar coef = az*az/(ax*ay);
        nz_ = label(Foam::pow(ncells*coef, 1.0/3.0) + 0.5);
        nx_ = label((ax/az)*nz_ + 0.5);
        ny_ = label((ay/az)*nz_ + 0.5);
    }

    Info<< "nx,ny,nz: " << nx_<< " " << ny_<< " " << nz_<<endl;

    numCells_.setSize(nx_ + 1, ny_ + 1, nz_ + 1);
    numCells_ = 0;

    deltaX_ = deltaX/scalar(nx_) + SMALL;
    deltaY_ = deltaY/scalar(ny_) + SMALL;
    deltaZ_ = deltaZ/scalar(nz_) + SMALL;

    Info<< "deltas: " << deltaX_ << " " << deltaY_ << " " << deltaZ_ << endl;

    ampcoef_ = Foam::max(ampcoef_, 1);
    ampcoef_ = Foam::min(ampcoef_, 1.5);

    //determine the number of cells in each grid
    forAll(mesh->C(), i)
    {
        const point& pt = mesh->C()[i];

        if (!contains(pt)) continue;

        label ix = label(Foam::min(labelMax, (pt.x() - pmin.x())/deltaX_));
        label iy = label(Foam::min(labelMax, (pt.y() - pmin.y())/deltaY_));
        label iz = label(Foam::min(labelMax, (pt.z() - pmin.z())/deltaZ_));

        ix = Foam::min(nx_ - 1, ix);
        iy = Foam::min(ny_ - 1, iy);
        iz = Foam::min(nz_ - 1, iz);

        //Info<<"ix iy iz:"<<ix<<" "<<iy<<" "<<iz<<endl;
        numCells_.setval(ix, iy, iz, true);
    }

    //second level
    const vectorField verts = mesh->points();

    forAll(mesh->cellPoints(), ic)
    {
        const labelUList& cverts = mesh->cellPoints()[ic];

        scalar xminc=  GREAT;
        scalar xmaxc= -GREAT;
        scalar yminc=  GREAT;
        scalar ymaxc= -GREAT;
        scalar zminc=  GREAT;
        scalar zmaxc= -GREAT;

        forAll(cverts, j)
        {
            label jv = cverts[j];

            const point& pv = verts[jv];

            xminc = Foam::min(xminc, pv.x());
            yminc = Foam::min(yminc, pv.y());
            zminc = Foam::min(zminc, pv.z());

            xmaxc = Foam::max(xmaxc, pv.x());
            ymaxc = Foam::max(ymaxc, pv.y());
            zmaxc = Foam::max(zmaxc, pv.z());
        }

        xmaxc = xminc + ampcoef_*(xmaxc - xminc);
        ymaxc = yminc + ampcoef_*(ymaxc - yminc);
        zmaxc = zminc + ampcoef_*(zmaxc - zminc);

        label i0 = label(Foam::max(0,   Foam::min(labelMax, (xminc - pmin.x())/deltaX_)));
        label i1 = label(Foam::min(nx_, Foam::max(labelMin, (xmaxc - pmin.x())/deltaX_)));

        label j0 = label(Foam::max(0,   Foam::min(labelMax, (yminc - pmin.y())/deltaY_)));
        label j1 = label(Foam::min(ny_, Foam::max(labelMin, (ymaxc - pmin.y())/deltaY_)));

        label k0 = label(Foam::max(0,   Foam::min(labelMax, (zminc - pmin.z())/deltaZ_)));
        label k1 = label(Foam::min(nz_, Foam::max(labelMin, (zmaxc - pmin.z())/deltaZ_)));

        for (label i = i0; i < i1; i++)
        {
            for (label j = j0; j < j1; j++)
            {
                for (label k = k0; k < k1; k++)
                {
                    if (!numCells_.getval(i, j, k))
                    {
                        numCells_.setval(i, j, k, true);
                    }
                }
            }
        }
    }
} //getCellList


void smBoundBox::createBox
(
    const fvMesh* mesh,
    label cellNum
)
{
    if (boxType_ == 1)
    {
        getCellList(mesh, cellNum);
    }

}


smBoundBox::smBoundBox
(
    const point &pmin,
    const point &pmax,
    const word type,/* = "boundBox" */
    scalar ampcoef   /* = 1.5*/
)
:
boundBox(pmin, pmax)
{
    defaults();

    if (type == "boundBox")
    {
        boxType_ = 0;
        ampcoef_ = 1;
    }
    else
    {
        boxType_ = 1;
        ampcoef_ = ampcoef;
    }
}


smBoundBox::smBoundBox
(
    const fvMesh* mesh,
    const word type,
    scalar ampcoef,
    label cellNum
)
:
boundBox
(
    mesh->C(),
    false
)
{
    defaults();

    ampcoef_ = ampcoef;

    if (type == "simple")
    {
        boxType_ = 0;
    }
    else
    {
        boxType_ = 1;
    }

    createBox(mesh, cellNum);
}
