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
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"
#include "mappedPatches/GAMG/mappedGAMGInterface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/HashTables/Map/Map.H"
#include "meshes/polyMesh/polyMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        mappedGAMGInterface,
        lduInterface
    );
}


void Foam::mappedGAMGInterface::makeIndicesUnique
(
    labelField& fullRestrictAddressing,
    const labelListList& constructMap
) const
{
    label maxSoFar = 0;
    label maxPrevProc = 0;
    forAll(constructMap, procI)
    {
        forAll(constructMap[procI], i)
        {
            fullRestrictAddressing[constructMap[procI][i]] +=
                maxPrevProc+1;
            maxSoFar = max
                (
                    maxSoFar,
                    fullRestrictAddressing[constructMap[procI][i]]
                );
        }
        maxPrevProc = maxSoFar;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedGAMGInterface::mappedGAMGInterface
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
    fineMappedLduInterface_
    (
        refCast<const mappedLduInterface>(fineInterface)
    ),
    nbrPatchID_(-1),
    nbrPatchSearched_(false),
    nbrAddrSize_(neighbourRestrictAddressing.size())
{
    bool coupled =
        (
            !neighbourRestrictAddressing.size() ||
            neighbourRestrictAddressing[0] != -1
        );
    reduce(coupled, andOp<bool>(), coarseComm);

    if (usingAMI() || !coupled)
    {
        // Construct face agglomeration from cell agglomeration
        {
            // From coarse face to cell
            DynamicList<label> dynFaceCells(localRestrictAddressing.size());

            // From face to coarse face
            DynamicList<label> dynFaceRestrictAddressing
            (
                localRestrictAddressing.size()
            );

            Map<label> masterToCoarseFace(localRestrictAddressing.size());

            forAll(localRestrictAddressing, ffi)
            {
                label curMaster = localRestrictAddressing[ffi];

                Map<label>::const_iterator fnd =
                    masterToCoarseFace.find
                    (
                        curMaster
                    );

                if (fnd == masterToCoarseFace.end())
                {
                    // New coarse face
                    label coarseI = dynFaceCells.size();
                    dynFaceRestrictAddressing.append(coarseI);
                    dynFaceCells.append(curMaster);
                    masterToCoarseFace.insert(curMaster, coarseI);
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
    }

    // On the owner side construct the AMI
    if (coupled && usingAMI() && owner())
    {
        // Construct the neighbour side agglomeration (as the neighbour would
        // do it so it the exact loop above using neighbourRestrictAddressing
        // instead of localRestrictAddressing)

        labelList nbrFaceRestrictAddressing;
        {
            // From face to coarse face
            DynamicList<label> dynNbrFaceRestrictAddressing
            (
                neighbourRestrictAddressing.size()
            );

            Map<label> masterToCoarseFace
            (
                neighbourRestrictAddressing.size()
            );

            forAll(neighbourRestrictAddressing, ffi)
            {
                label curMaster = neighbourRestrictAddressing[ffi];

                Map<label>::const_iterator fnd =
                    masterToCoarseFace.find
                    (
                        curMaster
                    );

                if (fnd == masterToCoarseFace.end())
                {
                    // New coarse face
                    label coarseI = masterToCoarseFace.size();
                    dynNbrFaceRestrictAddressing.append(coarseI);
                    masterToCoarseFace.insert(curMaster, coarseI);
                }
                else
                {
                    // Already have coarse face
                    dynNbrFaceRestrictAddressing.append(fnd());
                }
            }

            nbrFaceRestrictAddressing.transfer
            (
                dynNbrFaceRestrictAddressing
            );
        }

        amiPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                fineMappedLduInterface_.AMI(),
                faceRestrictAddressing_,
                nbrFaceRestrictAddressing
            )
        );
    }

    if (coupled && !usingAMI())
    {
        {
            // We need to get all slave-side agglomeration info onto our processor
            // so that we can generate the same number of coarse faces on each
            // side.
            labelField fullNeiRestrictAddressing(neighbourRestrictAddressing);
            // Map from both sides instead of using reverseDistribute, as gaps
            // are possible
            fineMappedLduInterface_.map().distribute
            (
                fullNeiRestrictAddressing
            );

            makeIndicesUnique
            (
                fullNeiRestrictAddressing,
                fineMappedLduInterface_.map().constructMap()
            );

            // Now perform agglomeration

            // From coarse face to coarse cell
            DynamicList<label> dynFaceCells(localRestrictAddressing.size());
            // From fine face to coarse face
            DynamicList<label> dynFaceRestrictAddressing
                (
                    localRestrictAddressing.size()
                );

            // From coarse cell pair to coarse face
            HashTable<label, labelPair, labelPair::Hash<>> cellsToCoarseFace
                (
                    2 * localRestrictAddressing.size()
                );

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
                        fullNeiRestrictAddressing[ffi]
                    );
                }
                else
                {
                    // Slave side
                    cellPair = labelPair
                    (
                        fullNeiRestrictAddressing[ffi],
                        localRestrictAddressing[ffi]
                    );
                }

                HashTable<label, labelPair, labelPair::Hash<>>::const_iterator
                    fnd = cellsToCoarseFace.find(cellPair);

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

        // Generate the processor map for the agglomeration
        {
            // Perform same agglomeration as if from the neighbour side to get
            // the neighbour's addressing restriction map.

            labelField nbrFaceRestrictAddressing;
            {
                labelField fullLocalRestrictAddressing(localRestrictAddressing);
                fineMappedLduInterface_.nbrPatch().map().distribute
                (
                    fullLocalRestrictAddressing
                );

                makeIndicesUnique
                (
                    fullLocalRestrictAddressing,
                    fineMappedLduInterface_.nbrPatch().map().constructMap()
                );

                // From coarse face to coarse cell
                DynamicList<label> dynFaceCells(neighbourRestrictAddressing.size());
                // From fine face to coarse face
                DynamicList<label> dynFaceRestrictAddressing
                (
                    neighbourRestrictAddressing.size()
                );

                // From coarse cell pair to coarse face
                HashTable<label, labelPair, labelPair::Hash<>> cellsToCoarseFace
                (
                    2 * neighbourRestrictAddressing.size()
                );

                forAll(neighbourRestrictAddressing, ffi)
                {
                    labelPair cellPair;

                    // Slave side
                    cellPair = labelPair
                        (
                            fullLocalRestrictAddressing[ffi],
                            neighbourRestrictAddressing[ffi]
                        );

                    HashTable<label, labelPair, labelPair::Hash<>>::const_iterator
                        fnd = cellsToCoarseFace.find(cellPair);

                    if (fnd == cellsToCoarseFace.end())
                    {
                        // New coarse face
                        label coarseI = dynFaceCells.size();
                        dynFaceRestrictAddressing.append(coarseI);
                        dynFaceCells.append(neighbourRestrictAddressing[ffi]);
                        cellsToCoarseFace.insert(cellPair, coarseI);
                    }
                    else
                    {
                        // Already have coarse face
                        dynFaceRestrictAddressing.append(fnd());
                    }
                }

                nbrFaceRestrictAddressing.transfer(dynFaceRestrictAddressing);
            }

            // subMap - send from other patch (this processor) to this
            // patch (all processors). Index into other patch.
            // constructMap - gather from other patch (all processors) to
            // this patch (this processor). Index into this patch.
            const labelListList& fineSubMap =
                fineMappedLduInterface_.map().subMap();
            const labelListList& fineConstructMap =
                fineMappedLduInterface_.map().constructMap();
            labelListList subMap(fineSubMap.size());
            labelListList constructMap(fineConstructMap.size());

            // This needs to be re-worked because it only works for 1-1 mapping.
            // Everything needs to be referenced to the agglomeration on this
            // side to differentiate between a face that is repeated in
            // the mapping and two faces that agglomerate to the same one
            forAll(subMap, procI)
            {
                DynamicList<label> dymCoarseSubMap;
                boolList coarseDone
                    (max(0, max(nbrFaceRestrictAddressing) + 1), false);
                boolList fineDone(nbrFaceRestrictAddressing.size(), false);
                forAll(fineSubMap[procI], i)
                {
                    label fineI = fineSubMap[procI][i];
                    label coarseI = nbrFaceRestrictAddressing[fineI];
                    if (!coarseDone[coarseI])
                    {
                        coarseDone[coarseI] = true;
                        dymCoarseSubMap.append(coarseI);
                    }
                    else if (fineDone[fineI])
                    {
                        FatalErrorInFunction
                            << "GAMG agglomeration works with mapped "
                            << "boundaries having 1-1 mapping only"
                            << exit(FatalError);
                    }
                    fineDone[fineI] = true;
                }
                subMap[procI].transfer(dymCoarseSubMap);
            }

            label constructSize = 0;
            forAll(constructMap, procI)
            {
                DynamicList<label> dymCoarseConstructMap;
                boolList coarseDone(size(), false);
                forAll(fineConstructMap[procI], i)
                {
                    label coarseI =
                        faceRestrictAddressing_[fineConstructMap[procI][i]];
                    if (!coarseDone[coarseI])
                    {
                        coarseDone[coarseI] = true;
                        dymCoarseConstructMap.append(coarseI);
                        constructSize = max(constructSize, coarseI + 1);
                    }
                }
                constructMap[procI].transfer(dymCoarseConstructMap);
            }

            mapPtr_.reset
            (
                new mapDistribute
                (
                    constructSize,
                    subMap.xfer(),
                    constructMap.xfer()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedGAMGInterface::~mappedGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::mappedGAMGInterface::nbrPatchID() const
{
    if (!nbrPatchSearched_)
    {
        nbrPatchID_ = -1;
        if (isInterface() && nbrMeshFound())
        {
            forAll(coarseInterfaces_, coarsePatchI)
            {
                if
                (
                    coarseInterfaces_.set(coarsePatchI) &&
                    isA<mappedGAMGInterface>
                    (
                        coarseInterfaces_[coarsePatchI]
                    )
                )
                {
                    const mappedGAMGInterface& coarseInterface
                    (
                        refCast<const mappedGAMGInterface>
                        (
                            coarseInterfaces_[coarsePatchI]
                        )
                    );
                    if
                    (
                        coarseInterface.patchName() == nbrPatchName() &&
                        coarseInterface.regionName() == nbrRegionName()
                    )
                    {
                        nbrPatchID_ = coarsePatchI;
                        break;
                    }
                }
            }
        }
        nbrPatchSearched_ = true;
    }
    return nbrPatchID_;
}


Foam::tmp<Foam::labelField> Foam::mappedGAMGInterface::
internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    // Must allow for no neighbour found as this is the case for fields solved
    // non-monolithically or is not an interface
    if (nbrPatchID() == -1)
    {
        // Just return a dummy field of size 1 as the AMG agglomerates
        // to 1 when the target isn't found
        tmp<labelField> tpnf(new labelField(1, -1));
        return tpnf;
    }
    else
    {
        const labelUList& nbrFaceCells = nbrPatch().faceCells();

        tmp<labelField> tpnf(new labelField(nbrFaceCells.size()));
        labelField &pnf = tpnf.ref();

        forAll(pnf, facei)
        {
            pnf[facei] = iF[nbrFaceCells[facei]];
        }

        return tpnf;
    }
}


// ************************************************************************* //
