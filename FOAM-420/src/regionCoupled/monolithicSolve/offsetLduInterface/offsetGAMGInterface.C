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
    (c) 2010-2015 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "offsetLduInterface/offsetGAMGInterface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(offsetGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        offsetGAMGInterface,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        offsetGAMGInterface,
        Istream
    );
}


// * * * * * * * * * * * Private Member functions  * * * * * * * * * * * * * //

void Foam::offsetGAMGInterface::grabSlaveCoarseInterfaces() const
{
    if (!slaveCoarseInterfaces_.size())
    {
        slaveCoarseInterfaces_.resize(coarseInterfaces_.size());
        forAll(coarseInterfaces_, i)
        {
            if (coarseInterfaces_.set(i))
            {
                slaveCoarseInterfaces_.set
                (
                    i,
                    &refCast<const offsetGAMGInterface>
                    (
                        coarseInterfaces_[i]
                    ).slaveInterface()
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::offsetGAMGInterface::offsetGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces
    ),
    slaveGAMGInterface_
    (
        GAMGInterface::New
        (
            index,
            slaveCoarseInterfaces_,
            refCast<const offsetLduInterfaceBase>(fineInterface).slaveInterface(),
            localRestrictAddressing,
            neighbourRestrictAddressing,
            fineLevelIndex,
            coarseComm
        )
    )
{}


Foam::offsetGAMGInterface::offsetGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    GAMGInterface(index, coarseInterfaces, is),
    slaveGAMGInterface_
    (
        GAMGInterface::New
        (
            word(is),
            index,
            slaveCoarseInterfaces_,
            is
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::offsetGAMGInterface::~offsetGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::offsetGAMGInterface::interfaceInternalField
(
    const labelUList& internalData
) const
{
    tmp<labelField> tlf =
        slaveGAMGInterface().interfaceInternalField(internalData);
    return tlf;
}


void Foam::offsetGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    slaveGAMGInterface().initInternalFieldTransfer(commsType, iF);
}


Foam::tmp<Foam::labelField> Foam::offsetGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    tmp<labelField> tlf =
        slaveGAMGInterface().internalFieldTransfer(commsType, iF);
    return tlf;
}


Foam::tmp<Foam::scalarField> Foam::offsetGAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    return slaveGAMGInterface().agglomerateCoeffs(fineCoeffs);
}


void Foam::offsetGAMGInterface::write(Ostream& os) const
{
    GAMGInterface::write(os);
    os << slaveGAMGInterface().type() << endl;
    slaveGAMGInterface().write(os);
}


// ************************************************************************* //
