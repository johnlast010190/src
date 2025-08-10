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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/findRefCell/findRefCell.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::setRefCell
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    if (fieldRef.needReference() || forceReference)
    {
        word refCellName = field.name() + "RefCell";
        word refProcessorName = field.name() + "RefProcessor";
        word refPointName = field.name() + "RefPoint";

        word refValueName = field.name() + "RefValue";

        if (dict.found(refCellName))
        {
            label refProc = dict.lookupOrDefault<label>(refProcessorName, 0);
            refCelli = readLabel(dict.lookup(refCellName));
            if (refCelli == -1)
            {
                // Special case: find first processor which has a cell
                // and use its cell 0 as ref
                refProc =
                (
                    field.mesh().nCells()
                    ? Pstream::myProcNo()
                    : Pstream::nProcs()
                );
                reduce(refProc, minOp<label>());
                if (Pstream::myProcNo() != refProc)
                {
                    refCelli = -1;
                }
                else
                {
                    refCelli = 0;
                }
            }
            if (Pstream::myProcNo() == refProc)
            {
                if (refCelli < 0 || refCelli >= field.mesh().nCells())
                {
                    FatalIOErrorInFunction
                    (
                        dict
                    )   << "Illegal cellID " << refCelli
                        << " on processor " << refProc
                        << ". Should be 0.." << field.mesh().nCells()
                        << exit(FatalIOError);
                }
            }
            else
            {
                refCelli = -1;
            }
        }
        else if (dict.found(refPointName))
        {
            point refPointi(dict.lookup(refPointName));

            // Try fast approximate search avoiding octree construction
            refCelli = field.mesh().findCell(refPointi, polyMesh::FACE_PLANES);

            label hasRef = (refCelli >= 0 ? 1 : 0);
            label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

            // If reference cell no found use octree search
            // with cell tet-decompositoin
            if (sumHasRef != 1)
            {
                refCelli = field.mesh().findCell(refPointi);

                hasRef = (refCelli >= 0 ? 1 : 0);
                sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
            }

            if (sumHasRef != 1)
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Unable to set reference cell for field " << field.name()
                    << nl << "    Reference point " << refPointName
                    << " " << refPointi
                    << " found on " << sumHasRef << " domains (should be one)"
                    << nl << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unable to set reference cell for field " << field.name()
                << nl
                << "    Please supply either " << refCellName
                << " or " << refPointName << nl << exit(FatalIOError);
        }

        refValue = readScalar(dict.lookup(refValueName));

        return true;
    }
    else
    {
        refCelli = -1;
        return false;
    }
}


bool Foam::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    return setRefCell(field, field, dict, refCelli, refValue, forceReference);
}


Foam::scalar Foam::getRefCellValue
(
    const volScalarField& field,
    const label refCelli
)
{
    scalar refCellValue = (refCelli >= 0 ? field[refCelli] : 0.0);
    return returnReduce(refCellValue, sumOp<scalar>());
}


// ************************************************************************* //
