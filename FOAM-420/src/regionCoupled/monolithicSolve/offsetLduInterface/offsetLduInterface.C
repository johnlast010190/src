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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "offsetLduInterface/offsetLduInterface.H"
#include "containers/Lists/SubList/SubList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(offsetLduInterface, 0);
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::offsetLduInterface::offsetLduInterface
(
    const lduInterface& slaveInterface,
    const label internalOffset,
    const label internalSize,
    const label foreignOffset,
    const label foreignSize
)
:
    slaveInterface_(slaveInterface),
    internalOffset_(internalOffset),
    internalSize_(internalSize),
    foreignOffset_(foreignOffset),
    foreignSize_(foreignSize)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::offsetLduInterface::~offsetLduInterface()
{}


// * * * * * * * * * * * * * * Public members  * * * * * * * * * * * * * * * //

const Foam::labelUList& Foam::offsetLduInterface::faceCells() const
{
    faceCells_ = slaveInterface_.faceCells();
    forAll(faceCells_, i)
    {
        faceCells_[i] += internalOffset_;
    }
    return faceCells_;
}

//TODO: for these three funcitons we should really offset the field before calling and then offset
// it back agian afterwards just in case the slave interface functions actually try to interpret the data
// as indices.
Foam::tmp<Foam::labelField> Foam::offsetLduInterface::interfaceInternalField
(
    const labelUList& internalData
) const
{
    tmp<labelField> tlf =
        slaveInterface_.interfaceInternalField
        (
            labelList::subList(internalData, internalSize_, internalOffset_)
        );
    return tlf;
}

void Foam::offsetLduInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    slaveInterface_.initInternalFieldTransfer
    (
        commsType,
        labelList::subList(iF, foreignSize_, foreignOffset_)
    );
}

Foam::tmp<Foam::labelField> Foam::offsetLduInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    tmp<labelField> tlf =
        slaveInterface_.internalFieldTransfer
        (
            commsType,
            labelList::subList(iF, foreignSize_, foreignOffset_)
        );
    return tlf;
}

// ************************************************************************* //
