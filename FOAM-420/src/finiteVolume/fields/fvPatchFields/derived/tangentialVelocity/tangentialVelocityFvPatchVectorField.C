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
    (c) 2010-2012 Esi Ltd.
    (c) 2006-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/tangentialVelocity/tangentialVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "containers/Lists/PackedList/PackedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tangentialVelocityFvPatchVectorField::
tangentialVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    magVel_(0),
    crossVectorPtr_(nullptr),
    directionPtr_(nullptr)
{}


Foam::tangentialVelocityFvPatchVectorField::
tangentialVelocityFvPatchVectorField
(
    const tangentialVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    magVel_(ptf.magVel_),
    crossVectorPtr_(nullptr),
    directionPtr_(nullptr)
{
    if (ptf.crossVectorPtr_.valid())
    {
        crossVectorPtr_.reset(new vector(ptf.crossVectorPtr_()));
    }
    else if (ptf.directionPtr_.valid())
    {
        directionPtr_.reset(new vector(ptf.directionPtr_()));
    }
}


Foam::tangentialVelocityFvPatchVectorField::
tangentialVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    magVel_(readScalar(dict.lookup("magnitude"))),
    crossVectorPtr_(nullptr),
    directionPtr_(nullptr)
{
    if (dict.found("crossVector"))
    {
        crossVectorPtr_.reset(new vector(dict.lookup("crossVector")));
    }
    else if (dict.found("direction"))
    {
        directionPtr_.reset(new vector(dict.lookup("direction")));
    }
}


Foam::tangentialVelocityFvPatchVectorField::
tangentialVelocityFvPatchVectorField
(
    const tangentialVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    magVel_(ptf.magVel_),
    crossVectorPtr_(nullptr),
    directionPtr_(nullptr)
{
    if (ptf.crossVectorPtr_.valid())
    {
        crossVectorPtr_.reset(new vector(ptf.crossVectorPtr_()));
    }
    else if (ptf.directionPtr_.valid())
    {
        directionPtr_.reset(new vector(ptf.directionPtr_()));
    }
}


Foam::tangentialVelocityFvPatchVectorField::
tangentialVelocityFvPatchVectorField
(
    const tangentialVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    magVel_(ptf.magVel_),
    crossVectorPtr_(nullptr),
    directionPtr_(nullptr)
{
    if (ptf.crossVectorPtr_.valid())
    {
        crossVectorPtr_.reset(new vector(ptf.crossVectorPtr_()));
    }
    else if (ptf.directionPtr_.valid())
    {
        directionPtr_.reset(new vector(ptf.directionPtr_()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tangentialVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();
    vectorField n(p.nf()());

    vectorField vdir(p.size(), Zero);

    if (crossVectorPtr_.valid())
    {
        vector globalCrossVector(crossVectorPtr_());
        makeVectorGlobal(globalCrossVector);
        vdir = globalCrossVector^n;
        vdir /= max(mag(vdir), SMALL);
        vdir -= n*(n & vdir);
    }
    else if (directionPtr_.valid())
    {
        vector globalDirection(directionPtr_());
        makeVectorGlobal(globalDirection);
        globalDirection = normalised(globalDirection);
        vdir = globalDirection - n*(n & globalDirection);
    }
    else
    {
        WarningInFunction
            << "No crossVector or direction specified."
            << exit(FatalError);
    }

    // tangential velocity
    vectorField ft(magVel_*vdir);
    frameFieldUpdate(ft, true, false);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::tangentialVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeEntry("magnitude", magVel_);

    if (crossVectorPtr_.valid())
    {
        os.writeEntry("crossVector", crossVectorPtr_());
    }
    else if (directionPtr_.valid())
    {
        os.writeEntry("direction", directionPtr_());
    }

    writeEntry("value", os);

    referenceFrameFvPatch<vector>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       tangentialVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
