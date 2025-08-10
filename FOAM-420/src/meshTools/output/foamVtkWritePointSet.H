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
    (c) 2011 OpenFOAM Foundation

InNamespace
    Foam::vtk

Description
    Write pointSet to vtk polydata file.
    The data are the original point ids.

SourceFiles
    foamVtkWritePointSet.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkWritePointSet_H
#define foamVtkWritePointSet_H

#include "vtk/output/foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class primitiveMesh;
class pointSet;
class fileName;

namespace vtk
{

//- Write pointSet to vtk polydata file.
//  The data are the original point ids.
void writePointSet
(
    const primitiveMesh& mesh,
    const pointSet& set,
    const fileName& baseName,
    const vtk::outputOptions outOpts
);


} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
