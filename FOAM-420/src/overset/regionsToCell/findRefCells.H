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

InNamespace
    Foam

Description
    Find the reference cells per region. Regions given as cellSet names.

SourceFiles
    findRefCell.C

SeeAlso
    Foam::regionSplit

\*---------------------------------------------------------------------------*/

#ifndef findRefCell_H
#define findRefCell_H

#include "primitives/ints/label/label.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Set reference cells for multi-region domains. Returns per region:
//  - whether a reference is needed
//  - if so a reference cell (or -1)
//  - and a reference value
void setRefCells
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,

    const label refCelli,
    const scalar refValue,

    boolList& needReference,
    labelList& refCells,
    scalarField& refValues,
    const bool forceReference = false
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
