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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    controlBox

Description
    Defines a control lattice

\*---------------------------------------------------------------------------*/

#include "vNurbsDeformation/controlBox/controlBox.H"
#include "db/IOobject/IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void controlBox::computePoints(const pointField& blockPoints)
{
    if (blockPoints.size() != 8)
        FatalErrorInFunction
            << "Block constructor must receive 8 points.\n"
            << exit(FatalError);

    const point& p0 = blockPoints[0];
    const point& p1 = blockPoints[1];
    const point& p2 = blockPoints[2];
    const point& p3 = blockPoints[3];
    const point& p4 = blockPoints[4];
    const point& p5 = blockPoints[5];
    const point& p6 = blockPoints[6];
    const point& p7 = blockPoints[7];

    // calculate two pointFields for top and bottom of lattice
    pointField top(nXCPs_*nYCPs_, vector::zero);
    pointField bot(nXCPs_*nYCPs_, vector::zero);

    for (label i = 0; i < nXCPs_; i++)
    {
        bot[i] = p0 + (p1-p0)*scalar(i)/scalar(nXCPs_-1);
        top[i] = p4 + (p5-p4)*scalar(i)/scalar(nXCPs_-1);

        bot[i + (nYCPs_-1)*nXCPs_] = p3 + (p2-p3)*scalar(i)/scalar(nXCPs_-1);
        top[i + (nYCPs_-1)*nXCPs_] = p7 + (p6-p7)*scalar(i)/scalar(nXCPs_-1);
    }

    for (label i = 0; i < nXCPs_; i++)
    {
        point bF = bot[i];
        point bL = bot[i + (nYCPs_-1)*nXCPs_];
        point tF = top[i];
        point tL = top[i + (nYCPs_-1)*nXCPs_];

        for (label j = 1; j < nYCPs_-1; j++)
        {
            bot[i + j*nXCPs_] = bF + (bL-bF)*scalar(j)/scalar(nYCPs_-1);
            top[i + j*nXCPs_] = tF + (tL-tF)*scalar(j)/scalar(nYCPs_-1);
        }
    }

    // fill points array using top & bottom

    for (label j = 0; j < nYCPs_; j++)
    {
        for (label i = 0; i < nXCPs_; i++)
        {
            point pF = bot[i+j*nXCPs_];
            point pL = top[i+j*nXCPs_];

            for (label k = 0; k < nZCPs_; k++)
            {
                points_[this->calculateIndex(i, j, k)] =
                    pF + (pL-pF)*scalar(k)/scalar(nZCPs_-1);
            }
        }
    }
}

void controlBox::init(const pointField& blockPoints)
{
        cBoxMesh = nullptr;
    // Initialize fields

    points_ = pointField(nXCPs_*nYCPs_*nZCPs_, vector::zero);
    displacements_ = pointField(nXCPs_*nYCPs_*nZCPs_, vector::zero);

    // Calculate point positions

    this->computePoints(blockPoints);
}

label controlBox::calculateIndex
(
    const label  iCP,
    const label  jCP,
    const label  kCP
) const
{
    label index = iCP + jCP*nXCPs_ + kCP*nXCPs_*nYCPs_;

    return index;
}

face controlBox::createFace
(
    const label& index0,
    const label& index1,
    const label& index2,
    const label& index3
) const
{
    labelList list(4, label(-1));
    list[0] = index0;
    list[1] = index1;
    list[2] = index2;
    list[3] = index3;

    return face(list);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

controlBox::controlBox
(
    const dictionary& dict
)
{

    nXCPs_ = readLabel(dict.lookup("nUCPs"));
    nYCPs_ = readLabel(dict.lookup("nVCPs"));
    nZCPs_ = readLabel(dict.lookup("nWCPs"));

    pointField blockPoints = dict.lookup("block");

    this->init(blockPoints);
}

controlBox::controlBox
(
    const controlBox& box
)
{
    nXCPs_ = box.nXCPs_;
    nYCPs_ = box.nYCPs_;
    nZCPs_ = box.nZCPs_;
    points_ = box.points();
    displacements_ = box.displacements();
    cBoxMesh = box.cBoxMesh;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const point& controlBox::getPoint
(
    const label  iCP,
    const label  jCP,
    const label  kCP
) const
{
    const label index = this->calculateIndex(iCP, jCP, kCP);

    const point& p = points_[index];

    return p;
}

const vector& controlBox::getDisplacement
(
    const label  iCP,
    const label  jCP,
    const label  kCP
) const
{
    const label index = this->calculateIndex(iCP, jCP, kCP);

    const vector& dv = displacements_[index];

    return dv;
}

void controlBox::updateBoxPosition()
{
    forAll(points_, pI)
    {
        vector dv = displacements_[pI];

    point p = points_[pI];

    p += dv;

    points_[pI] = p;

    displacements_[pI] = vector::zero;
    }

    cBoxMesh->movePoints(points_);
}

void controlBox::setDisplacement
(
    const label  iCP,
    const label  jCP,
    const label  kCP,
    const vector dV
)
{
    const label index = this->calculateIndex(iCP, jCP, kCP);

    displacements_[index] = dV;
}

void controlBox::setDisplacement
(
    const label  index,
    const vector dV
)
{
    displacements_[index] = dV;
}

void controlBox::addDisplacement
(
    const label  index,
    const vector dV
)
{
    displacements_[index] += dV;
}

void controlBox::boxToPolyMesh(const IOobject& io)
{
    bool inverted = false;

    label up = std::floor((nXCPs_-1)/2);
    label vp = std::floor((nYCPs_-1)/2);
    label wp = std::floor((nZCPs_-1)/2);

    Foam::vector du = this->getPoint(up+1, vp, wp)-this->getPoint(up, vp, wp);
    Foam::vector dv = this->getPoint(up, vp+1, wp)-this->getPoint(up, vp, wp);
    Foam::vector dw = this->getPoint(up, vp, wp+1)-this->getPoint(up, vp, wp);

    inverted = ( ( (du^dv)&dw ) < 0 );

    // Will use the points_ list as a base

    // Generate faces: Boundary faces will have a normal looking out of the box
    // Internal faces will have a normal looking to the cell with the larger label

    faceList fList;

    cellList cList;

    labelListList dummy;

    label index(0);

    for (label k = 0; k < nZCPs_-1; k++)
    {
        for (label j = 0; j < nYCPs_-1; j++)
        {
            for (label i = 0; i < nXCPs_-1; i++)
            {
                // grab 8 vertices

                label index0 = this->calculateIndex(i, j, k);
                label index1 = this->calculateIndex(i+1, j, k);
                label index2 = this->calculateIndex(i+1, j+1, k);
                label index3 = this->calculateIndex(i, j+1, k);
                label index4 = this->calculateIndex(i, j, k+1);
                label index5 = this->calculateIndex(i+1, j, k+1);
                label index7 = this->calculateIndex(i, j+1, k+1);

                // Make 6 faces looking out of the cell
                face f0 = this->createFace(index0, index3, index2, index1);
                face f1 = this->createFace(index0, index1, index5, index4);
                face f2 = this->createFace(index0, index4, index7, index3);

                if (inverted)
                {
                    f0.flip();
                    f1.flip();
                    f2.flip();
                }

                if (i != 0) f2.flip();
                if (j != 0) f1.flip();
                if (k != 0) f0.flip();

                labelList c;

                if (k != 0)
                {
                    fList.append(f0);
                    c.append(index);
                    index++;
                }
                if (j != 0)
                {
                    fList.append(f1);
                    c.append(index);
                    index++;
                }
                if (i != 0)
                {
                    fList.append(f2);
                    c.append(index);
                    index++;
                }

                dummy.append(c);
            }
        }
    }

    for (label k = 0; k < nZCPs_-1; k++)
    {
        for (label j = 0; j < nYCPs_-1; j++)
        {
            for (label i = 0; i < nXCPs_-1; i++)
            {
                // grab 8 vertices

                label index0 = this->calculateIndex(i, j, k);
                label index1 = this->calculateIndex(i+1, j, k);
                label index2 = this->calculateIndex(i+1, j+1, k);
                label index3 = this->calculateIndex(i, j+1, k);
                label index4 = this->calculateIndex(i, j, k+1);
                label index5 = this->calculateIndex(i+1, j, k+1);
                label index6 = this->calculateIndex(i+1, j+1, k+1);
                label index7 = this->calculateIndex(i, j+1, k+1);

                // Make 6 faces looking out of the cell
                face f0 = this->createFace(index0, index3, index2, index1);
                face f1 = this->createFace(index0, index1, index5, index4);
                face f2 = this->createFace(index0, index4, index7, index3);
                face f3 = this->createFace(index1, index2, index6, index5);
                face f4 = this->createFace(index4, index5, index6, index7);
                face f5 = this->createFace(index2, index3, index7, index6);

                if (inverted)
                {
                    f0.flip();
                    f1.flip();
                    f2.flip();
                    f3.flip();
                    f4.flip();
                    f5.flip();
                }

                labelList c = dummy[i+j*(nXCPs_-1)+k*(nXCPs_-1)*(nYCPs_-1)];

                if (i != 0) f2.flip();
                if (j != 0) f1.flip();
                if (k != 0) f0.flip();

                if (k == 0)
                {
                    fList.append(f0);
                    c.append(index);
                    index++;
                }
                if (j == 0)
                {
                    fList.append(f1);
                    c.append(index);
                    index++;
                }
                if (i == 0)
                {
                    fList.append(f2);
                    c.append(index);
                    index++;
                }

                if (k == nZCPs_-2)
                {
                    fList.append(f4);
                    c.append(index);
                    index++;
                }
                if (j == nYCPs_-2)
                {
                    fList.append(f5);
                    c.append(index);
                    index++;
                }
                if (i == nXCPs_-2)
                {
                    fList.append(f3);
                    c.append(index);
                    index++;
                }

                dummy[i+j*(nXCPs_-1)+k*(nXCPs_-1)*(nYCPs_-1)] = c;
            }
        }
    }

    for (label k = 0; k < nZCPs_-1; k++)
    {
        for (label j = 0; j < nYCPs_-1; j++)
        {
            for (label i = 0; i < nXCPs_-1; i++)
            {
                label cIndex = i+j*(nXCPs_-1)+k*(nXCPs_-1)*(nYCPs_-1);

                labelList c = dummy[cIndex];

                if (k != nZCPs_-2)
                {
                    label nIndex = i+j*(nXCPs_-1)+(k+1)*(nXCPs_-1)*(nYCPs_-1);


                    labelList nextCell = dummy[nIndex];

                    c.append(nextCell[0]);
                }
                if (j != nYCPs_-2)
                {
                    label nIndex = i+(j+1)*(nXCPs_-1)+k*(nXCPs_-1)*(nYCPs_-1);

                    labelList nextCell = dummy[nIndex];

                    label fj(0);

                    if (k != 0)
                        fj++;

                    c.append(nextCell[fj]);
                }
                if (i != nXCPs_-2)
                {
                    label nIndex = i+1+j*(nXCPs_-1)+k*(nXCPs_-1)*(nYCPs_-1);

                    labelList nextCell = dummy[nIndex];

                    label fi(0);

                    if (k != 0)
                        fi++;
                    if (j != 0)
                        fi++;

                    c.append(nextCell[fi]);
                }

                dummy[cIndex] = c;
            }
        }
    }


    // Generate cells

    forAll(dummy, dI)
    {
        cell c(dummy[dI]);
        cList.append(c);
    }

    Xfer <pointField> p(points_);
    Xfer <faceList> faces(fList);
    Xfer <cellList> cells(cList);

    cBoxMesh = new polyMesh
    (
        io,
        p,
        faces,
        cells
    );

    const pointField& points = cBoxMesh->points();
    maxpt_ = max(points);
    minpt_ = min(points);
    return;
}

const polyMesh* controlBox::getPolyMesh() const
{
    return cBoxMesh;
}

void controlBox::smoothDisplacements()
{
    for (label k = 2; k < nZCPs_-2; k++)
    {
        for (label j = 2; j < nYCPs_-2; j++)
        {
            for (label i = 2; i < nXCPs_-2; i++)
            {
                const point& p = points_[this->calculateIndex(i, j, k)];

                scalar l1jk = mag(p-points_[this->calculateIndex(i+1, j, k)]);
                scalar l0jk = mag(p-points_[this->calculateIndex(i-1, j, k)]);
                scalar li1k = mag(p-points_[this->calculateIndex(i, j+1, k)]);
                scalar li0k = mag(p-points_[this->calculateIndex(i, j-1, k)]);
                scalar lij1 = mag(p-points_[this->calculateIndex(i, j, k+1)]);
                scalar lij0 = mag(p-points_[this->calculateIndex(i, j, k-1)]);

                vector dv = displacements_[this->calculateIndex(i+1, j, k)]/l1jk
                          + displacements_[this->calculateIndex(i-1, j, k)]/l0jk
                          + displacements_[this->calculateIndex(i, j+1, k)]/li1k
                          + displacements_[this->calculateIndex(i, j-1, k)]/li0k
                          + displacements_[this->calculateIndex(i, j, k+1)]/lij1
                          + displacements_[this->calculateIndex(i, j, k-1)]/lij0;

                dv /= (1.0/l1jk + 1.0/l0jk + 1.0/li1k + 1.0/li0k + 1.0/lij1 + 1.0/lij0);

                displacements_[this->calculateIndex(i, j, k)] = dv;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
