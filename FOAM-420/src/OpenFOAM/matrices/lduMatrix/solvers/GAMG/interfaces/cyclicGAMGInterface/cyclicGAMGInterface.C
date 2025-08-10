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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/solvers/GAMG/interfaces/cyclicGAMGInterface/cyclicGAMGInterface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/Pair/labelPairHashes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        Istream
    );

    // Add under name cyclicSlip
    addNamedToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        lduInterface,
        cyclicSlip
    );
    addNamedToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        Istream,
        cyclicSlip
    );

    // Add under name nonConformalCyclic
    addNamedToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        lduInterface,
        nonConformalCyclic
    );
    addNamedToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        Istream,
        nonConformalCyclic
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterface::cyclicGAMGInterface
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
    GAMGInterface(index, coarseInterfaces),
    nbrPatchID_
    (
        refCast<const cyclicLduInterface>(fineInterface).nbrPatchID()
        // Add the offset from the patch numbers in the original matrix to those
        // in the super-matrix by looking at the difference in our own patch no.
      + (
            index
          - refCast<const cyclicLduInterface>
            (
                fineInterface
            ).nbrPatch().nbrPatchID()
        )
    ),
    owner_(refCast<const cyclicLduInterface>(fineInterface).owner()),
    transform_(refCast<const cyclicLduInterface>(fineInterface).transform())
{
    // From coarse face to coarse cell
    DynamicList<label> dynFaceCells(localRestrictAddressing.size());
    // From fine face to coarse face
    DynamicList<label> dynFaceRestrictAddressing
    (
        localRestrictAddressing.size()
    );

    // From coarse cell pair to coarse face
    labelPairLookup cellsToCoarseFace(2*localRestrictAddressing.size());

    forAll(localRestrictAddressing, ffi)
    {
        labelPair cellPair;

        // Do switching on master/slave indexes based on the owner/neighbour of
        // the processor index such that both sides get the same answer.
        if (owner())
        {
            // Master side
            cellPair = labelPair
            (
                localRestrictAddressing[ffi],
                neighbourRestrictAddressing[ffi]
            );
        }
        else
        {
            // Slave side
            cellPair = labelPair
            (
                neighbourRestrictAddressing[ffi],
                localRestrictAddressing[ffi]
            );
        }

        labelPairLookup::const_iterator fnd = cellsToCoarseFace.find(cellPair);

        if (fnd == cellsToCoarseFace.end())
        {
            // New coarse face
            label coarseI = dynFaceCells.size();
            dynFaceRestrictAddressing.append(coarseI);
            dynFaceCells.append(localRestrictAddressing[ffi]);
            cellsToCoarseFace.insert(cellPair, coarseI);
        }
        else
        {
            // Already have coarse face
            dynFaceRestrictAddressing.append(fnd());
        }
    }

    faceCells_.transfer(dynFaceCells);
    faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
}


Foam::cyclicGAMGInterface::cyclicGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    GAMGInterface(index, coarseInterfaces, is),
    nbrPatchID_(readLabel(is)),
    owner_(readBool(is)),
    transform_(is)
{}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterface::~cyclicGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::cyclicGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    const cyclicGAMGInterface& nbr = nbrPatch();
    const labelUList& nbrFaceCells = nbr.faceCells();

    tmp<labelField> tpnf(new labelField(size()));
    labelField& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = iF[nbrFaceCells[facei]];
    }

    return tpnf;
}


void Foam::cyclicGAMGInterface::write(Ostream& os) const
{
    GAMGInterface::write(os);
    os  << token::SPACE << nbrPatchID_
        << token::SPACE << owner_
        << token::SPACE << transform_;
}


// ************************************************************************* //
