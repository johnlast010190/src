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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvBoundaryMesh/fvBoundaryMesh.H"
#include "fvMesh/fvMesh.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvPatch, 0);
    defineRunTimeSelectionTable(fvPatch, polyPatch);
    addToRunTimeSelectionTable(fvPatch, fvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvPatch::fvPatch(const polyPatch& p, const fvBoundaryMesh& bm)
:
    polyPatch_(p),
    boundaryMesh_(bm),
    active_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvPatch::~fvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvPatch::constraintType(const word& pt)
{
    return fvPatchField<scalar>::patchConstructorTable_().count(pt) > 0;
}


Foam::wordList Foam::fvPatch::constraintTypes()
{
    wordList cTypes(polyPatchConstructorTable_().size());

    label i = 0;

    for
    (
        polyPatchConstructorTable::iterator cstrIter =
            polyPatchConstructorTable_().begin();
        cstrIter != polyPatchConstructorTable_().end();
        ++cstrIter
    )
    {
        if (constraintType(cstrIter->first))
        {
            cTypes[i++] = cstrIter->first;
        }
    }

    cTypes.setSize(i);

    return cTypes;
}


const Foam::vectorField& Foam::fvPatch::Cf() const
{
    return boundaryMesh().mesh().Cf().boundaryField()[index()];
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::Cn() const
{
    tmp<vectorField> tcc(new vectorField(size()));
    vectorField& cc = tcc.ref();

    const labelUList& faceCells = this->faceCells();

    // get reference to global cell centres
    const vectorField& gcc = boundaryMesh().mesh().cellCentres();

    forAll(faceCells, facei)
    {
        cc[facei] = gcc[faceCells[facei]];
    }

    return tcc;
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::nf() const
{
    return Sf()/magSf();
}


const Foam::vectorField& Foam::fvPatch::Sf() const
{
    return boundaryMesh().mesh().Sf().boundaryField()[index()];
}


const Foam::scalarField& Foam::fvPatch::magSf() const
{
    return boundaryMesh().mesh().magSf().boundaryField()[index()];
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::delta() const
{
    Switch version22
    (
        debug::debugSwitches().lookupOrDefault<Switch>
        ("version22", false)
    );

    tmp<vectorField> deltast(new vectorField(size(), vector::zero));
    vectorField& deltas = deltast.ref();

    const vectorField cn(Cn());
    const vectorField& cf = Cf();
    deltas = cf - cn;

    if (version22)
    {
        return deltast;
    }

    // Use patch-normal delta for all non-coupled BCs

    const vectorField nHat(nf());
    scalarField nHatCCf(nHat & deltas);

    // Calculate the normal deltas
    deltas = nHat*nHatCCf;

    // If the face has zero magnitute (all points at at the same line) - nHat=0
    // Then we use the standard expression
    // In the if-statement condition we use the scalar instead of mag(vector)
    // to avoid extra computational cost
    // this operation will be done rarely
    forAll(nHatCCf, fI)
    {
        if (nHatCCf[fI] == 0)
        {
            deltas[fI] = cf[fI] - cn[fI];
        }
    }

    return deltast;
}


void Foam::fvPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


void Foam::fvPatch::initMovePoints()
{}


void Foam::fvPatch::movePoints()
{}


const Foam::scalarField& Foam::fvPatch::deltaCoeffs() const
{
    return boundaryMesh().mesh().deltaCoeffs().boundaryField()[index()];
}


const Foam::scalarField& Foam::fvPatch::weights() const
{
    return boundaryMesh().mesh().weights().boundaryField()[index()];
}


// ************************************************************************* //
