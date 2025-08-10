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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvPatches/constraint/cyclicACMI/cyclicACMIFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvMesh.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicACMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicACMIFvPatch::updateAreas() const
{
    if (cyclicACMIPolyPatch_.updated())
    {
        if (debug)
        {
            Pout<< "cyclicACMIFvPatch::updateAreas() : updating fv areas for "
                << name() << " and " << this->nonOverlapPatch().name()
                << endl;
        }

        // owner couple
        const_cast<vectorField&>(Sf()) = patch().faceAreas();
        const_cast<scalarField&>(magSf()) = patch().magFaceAreas();

        // owner non-overlapping
        const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
        const_cast<vectorField&>(nonOverlapPatch.Sf()) =
            nonOverlapPatch.patch().faceAreas();
        const_cast<scalarField&>(nonOverlapPatch.magSf()) =
            nonOverlapPatch.patch().magFaceAreas();

        // neighbour couple
        const cyclicACMIFvPatch& nbrACMI = nbrPatch();
        const_cast<vectorField&>(nbrACMI.Sf()) =
            nbrACMI.patch().faceAreas();
        const_cast<scalarField&>(nbrACMI.magSf()) =
            nbrACMI.patch().magFaceAreas();

        // neighbour non-overlapping
        const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();
        const_cast<vectorField&>(nbrNonOverlapPatch.Sf()) =
            nbrNonOverlapPatch.patch().faceAreas();
        const_cast<scalarField&>(nbrNonOverlapPatch.magSf()) =
            nbrNonOverlapPatch.patch().magFaceAreas();

        // set the updated flag
        cyclicACMIPolyPatch_.setUpdated(false);
    }
}


void Foam::cyclicACMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();
        const scalarField deltas(nf() & coupledFvPatch::delta());

        // These deltas are of the cyclic part alone - they are
        // not affected by the amount of overlap with the nonOverlapPatch
        scalarField nbrDeltas
        (
            interpolate
            (
                nbrPatch.nf() & nbrPatch.coupledFvPatch::delta()
            )
        );

        scalar tol = cyclicACMIPolyPatch::tolerance();


        forAll(deltas, facei)
        {
            scalar di = deltas[facei];
            scalar dni = nbrDeltas[facei];

            if (dni < tol)
            {
                // Avoid zero weights on disconnected faces. This value
                // will be weighted with the (zero) face area so will not
                // influence calculations.
                w[facei] = 1.0;
            }
            else
            {
                w[facei] = dni/(di + dni);
            }
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicACMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


Foam::tmp<Foam::vectorField> Foam::cyclicACMIFvPatch::delta() const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();

        const vectorField patchD(coupledFvPatch::delta());

        vectorField nbrPatchD(interpolate(nbrPatch.coupledFvPatch::delta()));


        tmp<vectorField> tpdv(new vectorField(patchD.size()));
        vectorField& pdv = tpdv.ref();

        // Do the transformation if necessary
        if (transform().transforms())
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - transform().transform(dni);
            }
        }
        else
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - dni;
            }
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// ************************************************************************* //
