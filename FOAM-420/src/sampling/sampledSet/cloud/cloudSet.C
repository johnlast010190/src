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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSet/cloud/cloudSet.H"
#include "sampledSet/sampledSet/sampledSet.H"
#include "meshSearch/meshSearch.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/strings/word/word.H"
#include "fields/Fields/DynamicField/DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, cloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    const meshSearch& queryMesh = searchEngine();

    labelList foundProc(sampleCoords_.size(), -1);

    forAll(sampleCoords_, sampleI)
    {
        label celli = queryMesh.findCell(sampleCoords_[sampleI]);

        if (celli != -1)
        {
            samplingPts.append(sampleCoords_[sampleI]);
            samplingCells.append(celli);
            samplingFaces.append(-1);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0 * sampleI);

            foundProc[sampleI] = Pstream::myProcNo();
        }
    }

    // Check that all have been found
    labelList maxFoundProc(foundProc);
    Pstream::listCombineReduce(maxFoundProc, maxOp<label>());

    labelList minFoundProc(foundProc.size(), labelMax);
    forAll(foundProc, i)
    {
        if (foundProc[i] != -1)
        {
            minFoundProc[i] = foundProc[i];
        }
    }
    Pstream::listCombineReduce(minFoundProc, minOp<label>());

    DynamicField<point> missingPoints(sampleCoords_.size());

    forAll(sampleCoords_, sampleI)
    {
        if (maxFoundProc[sampleI] == -1)
        {
            // No processor has found the location.
            missingPoints.append(sampleCoords_[sampleI]);
        }
        else if (minFoundProc[sampleI] != maxFoundProc[sampleI])
        {
            WarningInFunction
                << "For sample set " << name()
                << " location " << sampleCoords_[sampleI]
                << " seems to be on multiple domains: "
                << minFoundProc[sampleI] << " and " << maxFoundProc[sampleI]
                << nl
                << "This might happen if the location is on"
                << " a processor patch. Change the location slightly"
                << " to prevent this." << endl;
        }
    }


    if (missingPoints.size() > 0)
    {
        if (missingPoints.size() < 100 || debug)
        {
            WarningInFunction
                << "For sample set " << name()
                << " did not found " << missingPoints.size()
                << " points out of " << sampleCoords_.size()
                << nl
                << "Missing points:" << missingPoints << endl;
        }
        else
        {
            WarningInFunction
                << "For sample set " << name()
                << " did not found " << missingPoints.size()
                << " points out of " << sampleCoords_.size()
                << nl
                << "Print missing points by setting the debug flag"
                << " for " << cloudSet::typeName << endl;
        }
    }
}


void Foam::cloudSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloudSet::cloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::cloudSet::cloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.lookup("points"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloudSet::~cloudSet()
{}


// ************************************************************************* //
