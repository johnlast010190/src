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

InNamespace
    Foam

Description
    The bandCompression function renumbers the addressing such that the
    band of the matrix is reduced. The algorithm uses a simple search
    through the neighbour list in order of connectivity.
    (CutHill-McKee algorithm)

SourceFiles
    bandCompression.C

\*---------------------------------------------------------------------------*/

#ifndef bandCompression_H
#define bandCompression_H

#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Renumbers the addressing to reduce the band of the matrix.
//  The algorithm uses a simple search through the neighbour list
//  Returns the order in which the cells need to be visited (i.e. ordered to
//  original)
labelList bandCompression(const labelListList& addressing);

//- Renumber with addressing in losort form (neighbour + start in neighbour)
labelList bandCompression(const labelList& cellCells, const labelList& offsets);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
