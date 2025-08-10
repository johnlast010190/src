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
    (c) 2018 OpenCFD Ltd.

Description
    Code chunk for converting face and point sets - included by foamToVTK.

\*---------------------------------------------------------------------------*/

// flag to top-level code to signal early stop.
bool wroteTopoSet = false;

// Write faceSet (as PolyData)
if (faceSetName.size())
{
    // Load
    faceSet set(mesh, faceSetName);

    fileName outputName
    (
        outputDir/regionPrefix/"face-set"
      / set.name()/set.name() + timeDesc
    );

    Info<< "    faceSet   : "
        << outputName << nl;

    vtk::writeFaceSet
    (
        mesh,
        set,
        writeOpts,
        outputName,
        Pstream::parRun()
    );

    wroteTopoSet = true;
}


// Write pointSet (as PolyData)
if (pointSetName.size())
{
    // Load
    pointSet set(mesh, pointSetName);

    fileName outputName
    (
        outputDir/regionPrefix/"point-set"
      / set.name()/set.name() + timeDesc
    );

    Info<< "    pointSet  : "
        << outputName << nl;

    vtk::writePointSet
    (
        mesh,
        set,
        writeOpts,
        outputName,
        Pstream::parRun()
    );

    wroteTopoSet = true;
}


// ************************************************************************* //
