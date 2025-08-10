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
    (c) 2017 OpenCFD Ltd.

InNamespace
    Foam::vtk

Description
    Write faces of cellSet to vtk polydata file.

    The data are the original cell ids

SourceFiles
    foamVtkWriteCellSetFaces.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkWriteCellSetFaces_H
#define foamVtkWriteCellSetFaces_H

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "meshes/primitiveMesh/primitivePatch/uindirectPrimitivePatch.H"
#include "vtk/output/foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class primitiveMesh;
class cellSet;
class fileName;

namespace vtk
{

//- Write perimeter faces of cellset to vtk polydata file.
//  The data are the original cell ids
void writeCellSetFaces
(
    const primitiveMesh& mesh,
    const cellSet& set,
    const fileName& baseName,
    const vtk::outputOptions outOpts,
    bool parallel = Pstream::parRun()
);

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
