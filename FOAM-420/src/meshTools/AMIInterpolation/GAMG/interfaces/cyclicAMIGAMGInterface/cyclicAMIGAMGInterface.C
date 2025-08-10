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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"
#include "AMIInterpolation/GAMG/interfaces/cyclicAMIGAMGInterface/cyclicAMIGAMGInterface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/HashTables/Map/Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicAMIGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicAMIGAMGInterface::cyclicAMIGAMGInterface
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
    fineCyclicAMIInterface_
    (
        refCast<const cyclicAMILduInterface>(fineInterface)
    ),
    superMatrixPatchOffset_
    (
        index - fineCyclicAMIInterface_.nbrPatch().nbrPatchID()
    )
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

            Map<label>::const_iterator fnd = masterToCoarseFace.find
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


    // On the owner side construct the AMI

    if (fineCyclicAMIInterface_.owner())
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

            Map<label> masterToCoarseFace(neighbourRestrictAddressing.size());

            forAll(neighbourRestrictAddressing, ffi)
            {
                label curMaster = neighbourRestrictAddressing[ffi];

                Map<label>::const_iterator fnd = masterToCoarseFace.find
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

            nbrFaceRestrictAddressing.transfer(dynNbrFaceRestrictAddressing);
        }

        amiPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                fineCyclicAMIInterface_.AMI(),
                faceRestrictAddressing_,
                nbrFaceRestrictAddressing
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::cyclicAMIGAMGInterface::~cyclicAMIGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::cyclicAMIGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    const cyclicAMIGAMGInterface& nbr =
        dynamic_cast<const cyclicAMIGAMGInterface&>(nbrPatch());
    const labelUList& nbrFaceCells = nbr.faceCells();

    tmp<labelField> tpnf(new labelField(nbrFaceCells.size()));
    labelField& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = iF[nbrFaceCells[facei]];
    }

    return tpnf;
}


// ************************************************************************* //
