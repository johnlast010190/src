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
    (c) 2011-2013 OpenFOAM Foundation

Application
    pdfPlot

Group
    grpPostProcessingUtilitie

Description
    Generates a graph of a probability distribution function.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "distributionModel/distributionModel.H"
#include "graphField/makeGraph.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "createFields.H"

    label iCheck = 100;
    for (label i=1; i<=nSamples; i++)
    {
        scalar ps = p->sample();
        label n = label((ps - xMin)*nIntervals/(xMax - xMin));
        samples[n]++;

        if (writeData)
        {
            filePtr() << ps << nl;
        }

        if (i % iCheck == 0)
        {
            Info<< "    processed " << i << " samples" << endl;

            if (i == 10*iCheck)
            {
                iCheck *= 10;
            }
        }
    }

    scalarField x(nIntervals);

    forAll(x, i)
    {
        x[i] = xMin + i*(xMax - xMin)/(nIntervals - 1);
    }

    makeGraph(x, samples, p->type(), pdfPath, runTime.graphFormat());

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
