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
    (c) 2012 OpenFOAM Foundation

Description
    Test the construction, insertion and removal of points from the dynamic
    indexed octree.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "cpuTime/cpuTime.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "algorithms/dynamicIndexedOctree/dynamicIndexedOctree.H"
#include "algorithms/dynamicIndexedOctree/dynamicTreeDataPoint.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    cpuTime timer;

    treeBoundBox overallBb
    (
        point(0,0,0),
        point(1,1,1)
    );

    label size = 5e6;

    Info<< "Creating indexed octrees with " << size
        << " points in bounding box: " << overallBb << endl;


// Create the data

    DynamicList<point> pointList(size);
    pointField pointFieldList(size);

    for (label pI = size - 1; pI >= 0; --pI)
    {
        scalar factor = pI/scalar(size);

        pointList.append(0.99*point(factor, factor, factor));
        pointFieldList[pI] = 0.99*point(factor, factor, factor);
    }

    for (label i=0; i<5; ++i)
    {
        pointList.append(point(0.95, 0.95,0.95));
        pointFieldList.append(point(0.95, 0.95,0.95));
    }

    Info<< "Time to construct lists of points: "
        << timer.cpuTimeIncrement() << endl;


// Construct the trees

    dynamicIndexedOctree<dynamicTreeDataPoint> tree
    (
        dynamicTreeDataPoint(pointList),
        overallBb,  // overall search domain
        20,         // max levels
        100,        // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity; n/a since no bounding boxes.
    );

    Info<< "Time to construct dynamic tree:    "
        << timer.cpuTimeIncrement() << endl;

    indexedOctree<treeDataPoint> tree2
    (
        treeDataPoint(pointFieldList),
        overallBb,  // overall search domain
        20,         // max levels
        100,       // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity; n/a since no bounding boxes.
    );

    Info<< "Time to construct normal tree:     "
        << timer.cpuTimeIncrement() << endl;

    Info<< "Statistics of the dynamic tree:" << endl;
    tree.writeTreeInfo();


// Test calls to findNearest()

    point p(0.94,0.94,0.94);

    Info<< "Nearest point to " << p << " = "
        << tree.findNearest(p, 0.4) << endl;

    for (label i = 0; i < size; ++i)
    {
        tree.findNearest
        (
            point(scalar(i)/size, scalar(i)/size, scalar(i)/size),
            0.4
        );
    }

    Info<< "Time to perform " << size
        << " findNearest() calls on the dynamic tree: "
        << timer.cpuTimeIncrement() << endl;

    for (label i = 0; i < size; ++i)
    {
        tree2.findNearest
        (
            point(scalar(i)/size, scalar(i)/size, scalar(i)/size),
            0.4
        );
    }

    Info<< "Time to perform " << size
        << " findNearest() calls on the normal tree: "
        << timer.cpuTimeIncrement() << endl;


// Test point insertion

    label index = pointList.size();
    pointList.append(p);

    Info<< nl << "Inserting point " << p << " with index " << index << endl;

    // Now insert a point into the tree.
    tree.insert(index, pointList.size());

    Info<< "Nearest point to " << p << " = "
        << tree.findNearest(p, 0.4) << endl;

    index = pointList.size();
    pointList.append(p);

    Info<< "Inserting same point " << p << " with index " << index << endl;

    tree.insert(index, pointList.size());

    Info<< "Nearest point to " << p << " = "
        << tree.findNearest(p, 0.4) << endl;

    Info<< nl << "Tree after insertion:" << endl;
    tree.writeTreeInfo();


// Test point removal

    label nRemoved = 0;
    for (label pI = size - 5; pI >= 5; pI--)
    {
        tree.remove(pI);
        ++nRemoved;
    }

    Info<< "Time to remove " << nRemoved << " points from the tree: "
        << timer.cpuTimeIncrement() << endl;

    Info<< nl << "Tree after removal:" << endl;
    tree.writeTreeInfo();

    Info<< "Nearest point to " << p << " = "
        << tree.findNearest(p, 0.4) << endl;

    return 0;
}


// ************************************************************************* //
