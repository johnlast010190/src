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
    (c) 2017-2018 OpenCFD Ltd.

InNamespace
    Foam::vtk

Description
    Write topoSet in VTK format

SourceFiles
    foamVtkWriteTopoSet.C
    foamVtkWriteCellSetFaces.C
    foamVtkWriteFaceSet.C
    foamVtkWritePointSet.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkWriteTopoSet_H
#define foamVtkWriteTopoSet_H

#include "primitives/strings/fileName/fileName.H"
#include "db/IOstreams/Pstreams/UPstream.H"
#include "vtk/output/foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMesh;
class topoSet;
class faceSet;
class cellSet;
class pointSet;

namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Dispatch to vtk::writeCellSetFaces, vtk::writeFaceSet, vtk::writePointSet
//  The file name is with/without an extension.
//  \return True on successful dispatch
bool writeTopoSet
(
    const polyMesh& mesh,
    const topoSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = Pstream::parRun()
);


//- Write faceSet as VTK polydata file.
//  Only one CELL_DATA, which is the original faceID
//  The file name is with/without an extension.
bool writeFaceSet
(
    const polyMesh& mesh,
    const faceSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = Pstream::parRun()
);


//- Write perimeter faces of cellset to vtk polydata file.
//  The data are the original cell ids
//  The file name is with/without an extension.
bool writeCellSetFaces
(
    const polyMesh& mesh,
    const cellSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = Pstream::parRun()
);


//- Write pointSet to VTK polydata file.
//  Only one CELL_DATA, which is the original pointID
//  The file name is with/without an extension.
bool writePointSet
(
    const polyMesh& mesh,
    const pointSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = Pstream::parRun()
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
