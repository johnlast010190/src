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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, mappedFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::mappedFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return this->patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::mappedFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    if (isInterface())
    {
        if (iF.size() != nbrPatch().boundaryMesh().mesh().cells().size())
        {
            // We are being called from a segregated solver; coupling should not
            // be used so just return dummy array of the correct length
            return tmp<labelField>(new labelField(nbrPatch().size(), -1));
        }
        else
        {
            return nbrPatch().patchInternalField(iF);
        }
    }
    else
    {
        // Coupling BC should not be used for this type of sampling, or being
        // called by segregated solver where neighbour mesh does not exist.
        // Return dummy array
        return tmp<labelField>(new labelField(this->size(), -1));
    }

}

// ************************************************************************* //
