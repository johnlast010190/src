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
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/flowRateInletVelocity/flowRateInletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "primitives/one/one.H"
#include "global/constants/mathematical/mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRateBase(patch(), db(), false),
    rotating_(false),
    axis_(),
    origin_(),
    angle_(),
    extrapolateProfile_(false),
    relax_(0.2),
    fluxCorrection_(false)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    flowRateBase(patch(), db(), dict, false),
    rotating_(dict.lookupOrDefault<Switch>("rotating", false)),
    extrapolateProfile_
    (
        dict.lookupOrDefault<Switch>("extrapolateProfile", false)
    ),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.2)),
    fluxCorrection_
    (
        dict.lookupOrDefault<Switch>("fluxCorrection", false)
    )
{

    // Value field require if mass based
    if (dict.found("value"))
    {
        forceAssign(vectorField("value", dict, p.size()));
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }

    // Load rotating settings
    if (rotating_)
    {
        axis_.reset(new vector(dict.lookup("axis")));
        axis_() /= mag(axis_());
        origin_.reset(new vector(dict.lookup("origin")));
        angle_ = Function1<scalar>::New("angle", dict);
    }
}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flowRateBase(patch(), db(), ptf, false),
    rotating_(ptf.rotating_),
    axis_(),
    origin_(),
    angle_(),
    extrapolateProfile_(ptf.extrapolateProfile_),
    relax_(ptf.relax_),
    fluxCorrection_(ptf.fluxCorrection_)
{
    if (ptf.rotating_)
    {
        axis_.reset(new vector(ptf.axis_()));
        origin_.reset(new vector(ptf.origin_()));
        angle_.reset(ptf.angle_().clone().ptr());
    }
}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    flowRateBase(patch(), db(), ptf, false),
    rotating_(ptf.rotating_),
    axis_(),
    origin_(),
    angle_(),
    extrapolateProfile_(ptf.extrapolateProfile_),
    relax_(ptf.relax_),
    fluxCorrection_(ptf.fluxCorrection_)
{
    if (ptf.rotating_)
    {
        axis_.reset(new vector(ptf.axis_()));
        origin_.reset(new vector(ptf.origin_()));
        angle_.reset(ptf.angle_().clone().ptr());
    }
}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    flowRateBase(patch(), db(), ptf, false),
    rotating_(ptf.rotating_),
    axis_(),
    origin_(),
    angle_(),
    extrapolateProfile_(ptf.extrapolateProfile_),
    relax_(ptf.relax_),
    fluxCorrection_(ptf.fluxCorrection_)
{
    if (ptf.rotating_)
    {
        axis_.reset(new vector(ptf.axis_()));
        origin_.reset(new vector(ptf.origin_()));
        angle_.reset(ptf.angle_().clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField::autoMap(m);
}


void Foam::flowRateInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField::rmap(ptf, addr);
}


void Foam::flowRateInletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField::autoMapGIB(mapper);
}


Foam::tmp<Foam::vectorField>
Foam::flowRateInletVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<vectorField> vict(new vectorField(this->size(), Zero));
    vectorField& vic = vict.ref();

    vic = neg(currentFlowRate());

    return vict;
}


Foam::tmp<Foam::vectorField>
Foam::flowRateInletVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<vectorField> vbct(new vectorField(*this));
    vectorField& vbc = vbct.ref();

    scalar curFlowRate = currentFlowRate();

    if (curFlowRate < 0)
    {
        if (fluxCorrection_)
        {
            vbc -= patchInternalField();
        }
        else
        {
            vbc = Zero;
        }
    }

    return vbct;
}


Foam::tmp<Foam::vectorField>
Foam::flowRateInletVelocityFvPatchVectorField::valueDivInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<vectorField>(new vectorField(this->size(), Zero));
}


Foam::tmp<Foam::vectorField>
Foam::flowRateInletVelocityFvPatchVectorField::valueDivBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<vectorField>(new vectorField(*this));
}


void Foam::flowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<vectorField> n = patch().nf();
    if (extrapolateProfile_)
    {
        vectorField newValues(this->patchInternalField());

        makeRelative(newValues);

        scalar flowRate = currentFlowRate();

        scalar estimatedFlowRate =
            -gSum(multiplyByRhoOrOne(this->patch().Sf() & newValues, false));

        if (estimatedFlowRate/flowRate > 0.5)
        {
            newValues *= (mag(flowRate)/mag(estimatedFlowRate));
        }
        else
        {
            newValues -=
                (
                    (
                        flowRate - estimatedFlowRate
                    )/gSum(multiplyByRhoOrOne(patch().magSf(), false))
                )*n;
        }

        forceAssign(newValues);
    }
    else
    {
        const scalar avgU =
           -currentFlowRate()
           /gSum(multiplyByRhoOrOne(patch().magSf(), true));
        forceAssign(n*avgU);
    }

    if (rotating_)
    {
        vectorField radial(this->patch().Cf() - origin_());
        //radial -= (radial & axis_()) * axis_();

        // Get angle
        scalar angle = angle_->value(db().time().timeOutputValue());
        angle *= constant::mathematical::pi/180.0;

        // Tangential magnitude
        scalarField tanMag(mag(*this)*tan(angle));

        // Rotational unit vector
        vectorField tanDir(axis_() ^ radial);
        tanDir /= mag(tanDir);

        // Remove normal component that could change massflow rate
        tanDir -= patch().nf()*(tanDir & patch().nf());

        forceAssign(*this + tanMag*tanDir);
    }
    else
    {
        vectorField Up(*this);
        frameFieldUpdate(Up);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletVelocityFvPatchVectorField::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>& bEq
) const
{
    bool isNegative = neg(currentFlowRate());
    if (relax_ == 1 || !isNegative)
    {
        return;
    }
    bEq.boundaryRelax(relax_, this->patch().index());
}


void Foam::flowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    flowRateBase::write(os);

    os.writeEntry("extrapolateProfile", extrapolateProfile_);
    os.writeEntry("rotating", rotating_);

    if (axis_.valid())
    {
        os.writeEntry("axis", axis_());
    }

    if (origin_.valid())
    {
        os.writeEntry("origin", origin_());
    }

    if (angle_.valid())
    {
        angle_->writeData(os);
    }

    writeEntryIfDifferent<scalar>(os, "relax", 0.2, relax_);
    writeEntryIfDifferent<Switch>
    (
        os,
        "fluxCorrection",
        false,
        fluxCorrection_
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        flowRateInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
