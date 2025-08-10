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
    (c) 2011 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/pressureCylindricalInletVelocity/pressureCylindricalInletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "global/constants/mathematical/mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureCylindricalInletVelocityFvPatchVectorField::
pressureCylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureDirectedInletVelocityFvPatchVectorField(p, iF),
    origin_(vector::zero),
    axis_(vector(0,0,1)),
    cylindricalVector_(p.size()),
    relax_(0.25)
{}


Foam::pressureCylindricalInletVelocityFvPatchVectorField::
pressureCylindricalInletVelocityFvPatchVectorField
(
    const pressureCylindricalInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    pressureDirectedInletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    cylindricalVector_(mapper(ptf.cylindricalVector_)),
    relax_(ptf.relax_)
{}


Foam::pressureCylindricalInletVelocityFvPatchVectorField::
pressureCylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    pressureDirectedInletVelocityFvPatchVectorField(p, iF),
    origin_(dict.lookupOrDefault<vector>("origin", vector::zero)),
    axis_(dict.lookupOrDefault<vector>("axis", vector(0,0,1))),
    cylindricalVector_("cylindricalInletDirection", dict, p.size()),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.25))
{

    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    //populate base class protected data members
    {
        convertCylindricalToCartesian();

        phiName_ = dict.lookupOrDefault<word>("phi", "phi");
        rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    }

    //fvPatchVectorField::operator=(inletDir_);
}


Foam::pressureCylindricalInletVelocityFvPatchVectorField::
pressureCylindricalInletVelocityFvPatchVectorField
(
    const pressureCylindricalInletVelocityFvPatchVectorField& ptf
)
:
    pressureDirectedInletVelocityFvPatchVectorField(ptf),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    cylindricalVector_(ptf.cylindricalVector_),
    relax_(ptf.relax_)
{}


Foam::pressureCylindricalInletVelocityFvPatchVectorField::
pressureCylindricalInletVelocityFvPatchVectorField
(
    const pressureCylindricalInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureDirectedInletVelocityFvPatchVectorField(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    cylindricalVector_(ptf.cylindricalVector_),
    relax_(ptf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureCylindricalInletVelocityFvPatchVectorField::
convertCylindricalToCartesian()
{
    //1. convert each face centre to cylindrical coordinates
    //2. obtain absolute cylindrical coordinates for inlet vector
    //3. convert absolute cylindrical to absolute cartesian
    //4. calculate relative cartesian vector norm

    //normalise axis
    vector e3 = axis_/mag(axis_);

    vectorField Cfr(patch().Cf());
    Cfr -= origin_;

    scalarField z (Cfr & e3);

    vectorField rxy (Cfr - e3*z);

    scalarField r (mag(rxy));

    //find the most aligned cartesian axis, flip as necessary
    vector e1(vector(1,0,0));
    {
        bool magXgtmagY = mag(e3.x()) > mag(e3.y());
        bool magXgtmagZ = mag(e3.x()) > mag(e3.z());
        bool magYgtmagZ = mag(e3.y()) > mag(e3.z());


        if (magXgtmagZ || magYgtmagZ)
        {
            if (magXgtmagY)
            {
                if (e3.x() >= 0)
                {
                    e1 = vector(0,1,0);
                }
                else
                {
                    e1 = vector(0,0,-1);
                }
            }
            else
            {
                if (e3.y() >= 0)
                {
                    e1 = vector(0,0,1);
                }
                else
                {
                    e1 = vector(-1,0,0);
                }
            }
        }
        else
        {
            if (e3.z() < 0)
            {
                e1 = vector(0,-1,0);
            }
        }
    }

    e1 -= e3 * (e1 & e3);
    e1 /= mag(e1);

    vector e2 = (e3^e1);
    e2 /= mag(e2);

    //theta proceed from e1 via right hand rule
    scalarField theta(this->size(), 0.0);
    {
        scalarField e1xy(rxy & e1);
        scalarField e2xy(rxy & e2);

        theta = acos(e1xy/mag(rxy));

        theta += neg(e2xy)*2*(constant::mathematical::pi - theta);
    }

    scalar degsToRads(constant::mathematical::pi/180);

    forAll(theta, i)
    {
        z[i] += cylindricalVector_[i].z();
        r[i] += cylindricalVector_[i].x();
        theta[i] += cylindricalVector_[i].y()*degsToRads;

        vector local(r[i] * cos(theta[i]), r[i] * sin(theta[i]), z[i]);

        inletDir_[i] = (tensor(e1, e2, e3) & local) - Cfr[i];

        inletDir_[i] /= mag(inletDir_[i]);
    }

}

void Foam::pressureCylindricalInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    pressureDirectedInletVelocityFvPatchVectorField::autoMap(m);
    m(cylindricalVector_, cylindricalVector_);
}


void Foam::pressureCylindricalInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    pressureDirectedInletVelocityFvPatchVectorField::rmap(ptf, addr);

    const pressureCylindricalInletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureCylindricalInletVelocityFvPatchVectorField>(ptf);

    cylindricalVector_.rmap(tiptf.cylindricalVector_, addr);
}


void Foam::pressureCylindricalInletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    pressureDirectedInletVelocityFvPatchVectorField::autoMapGIB(mapper);
    mapper.map(cylindricalVector_, vector::zero);
}


void Foam::pressureCylindricalInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    tmp<vectorField> n = patch().nf();

    tmp<vectorField> nid = n() * (n() & inletDir_);
    tmp<vectorField> tid = inletDir_ - nid();

    tid.ref() /= mag(nid());
    nid.ref() /= mag(nid());

    tmp<vectorField> Utan = (*this & tid())/magSqr(tid())*tid();

    //tmp<scalarField> ndmagS = (n & inletDir_)*patch().magSf();

    scalarField Uflux;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        Uflux = -phip/patch().magSf();
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                db(),
                rhoName_
            );

        Uflux = -phip/(rhop*patch().magSf());
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    //relax tangential component update
    //improves stability for high angle inlets

    forceAssign(Uflux*(nid() + relax_*tid()) + (1-relax_)*Utan());

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pressureCylindricalInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    //write non-default cylindrical coordinate properties
    writeEntryIfDifferent<vector>(os, "origin", vector::zero, origin_);
    writeEntryIfDifferent<vector>(os, "axis", vector(0, 0, 1), axis_);
    cylindricalVector_.writeEntry("cylindricalInletDirection", os);
    writeEntryIfDifferent<scalar>(os, "relax", 0.25, relax_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureCylindricalInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
