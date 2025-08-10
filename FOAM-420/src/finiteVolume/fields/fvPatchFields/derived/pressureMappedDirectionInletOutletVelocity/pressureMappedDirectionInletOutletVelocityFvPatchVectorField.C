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
    (c) 2017 OpenCFD Ltd.
    (c) 2020-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/pressureMappedDirectionInletOutletVelocity/pressureMappedDirectionInletOutletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/IOstreams/Fstreams/IFstream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
checkAndMap()
{
    // Initialise
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        fileName samplePointsFile
        (
            this->db().time().path()
           /this->db().time().caseConstant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        pointField samplePoints((IFstream(samplePointsFile)()));

        DebugInformation
            << "pressureMappedDirectionInletOutletVelocityFvPatchVectorField :"
            << " Read " << samplePoints.size() << " sample points from "
            << samplePointsFile << endl;


        if (cylindricalCoords_)
        {
            // x is r(m), y is theta(radians) and z is z(m)
            forAll(samplePoints, pI)
            {
                point& pntI = samplePoints[pI];

                scalar x = pntI.x() * Foam::cos(pntI.y());
                scalar y = pntI.x() * Foam::sin(pntI.y());

                pntI.x() = x;
                pntI.y() = y;
            }

            DebugInformation
                << "timeVaryingMappedFixedValueFvPatchField :"
                << " Converting cylindrical to cartesian coordinates "
                << endl;
        }


        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                this->patch().patch().faceCentres(),
                perturb_,
                delabella_,
                nearestOnly
            )
        );
    }

    fileName valsFile
    (
        this->db().time().path()
       /this->db().time().caseConstant()
       /"boundaryData"
       /this->patch().name()
       /fieldTableName_+"dir"
    );

    Field<vector> vals;

    IFstream(valsFile).operator()() >> vals;

    if (vals.size() != mapperPtr_().sourceSize())
    {
        FatalErrorInFunction
            << "Number of values (" << vals.size()
            << ") differs from the number of points ("
            <<  mapperPtr_().sourceSize()
            << ") in file " << valsFile << exit(FatalError);
    }

    inletDir_ = mapperPtr_().interpolate(vals);

    //- Normalisation
    forAll(inletDir_, fI)
    {
        scalar magDir = mag(inletDir_[fI]);
        if (magDir>SMALL)
        {
            inletDir_[fI] /= magDir;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
pressureMappedDirectionInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    fieldTableName_(iF.name()),
    phiName_("phi"),
    rhoName_("rho"),
    inletDir_(p.size()),
    perturb_(0),
    delabella_(true),
    mapMethod_(word::null),
    cylindricalCoords_(false),
    mapperPtr_(nullptr),
    preserveSign_(false)
{
    checkAndMap();
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
pressureMappedDirectionInletOutletVelocityFvPatchVectorField
(
    const pressureMappedDirectionInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldTableName_(ptf.fieldTableName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletDir_(mapper(ptf.inletDir_)),
    perturb_(ptf.perturb_),
    delabella_(ptf.delabella_),
    mapMethod_(ptf.mapMethod_),
    cylindricalCoords_(ptf.cylindricalCoords_),
    mapperPtr_(nullptr),
    preserveSign_(ptf.preserveSign_)
{
    checkAndMap();
}


Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
pressureMappedDirectionInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF, dict, false),
    fieldTableName_(iF.name()),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    inletDir_(this->size(), vector::zero),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    delabella_(dict.lookupOrDefault<Switch>("delabella", true)),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    cylindricalCoords_(dict.lookupOrDefault("cylindricalCoords", false)),
    mapperPtr_(nullptr),
    preserveSign_(dict.lookupOrDefault("preserveSign", false))
{
    dict.readIfPresent("fieldTable", fieldTableName_);
    patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
    checkAndMap();
}


Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
pressureMappedDirectionInletOutletVelocityFvPatchVectorField
(
    const pressureMappedDirectionInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    fieldTableName_(pivpvf.fieldTableName_),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_),
    perturb_(pivpvf.perturb_),
    delabella_(pivpvf.delabella_),
    mapMethod_(pivpvf.mapMethod_),
    cylindricalCoords_(pivpvf.cylindricalCoords_),
    mapperPtr_(nullptr),
    preserveSign_(pivpvf.preserveSign_)
{
    checkAndMap();
}


Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
pressureMappedDirectionInletOutletVelocityFvPatchVectorField
(
    const pressureMappedDirectionInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    fieldTableName_(pivpvf.fieldTableName_),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_),
    perturb_(pivpvf.perturb_),
    delabella_(pivpvf.delabella_),
    mapMethod_(pivpvf.mapMethod_),
    cylindricalCoords_(pivpvf.cylindricalCoords_),
    mapperPtr_(nullptr),
    preserveSign_(pivpvf.preserveSign_)
{
    checkAndMap();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
    mapperPtr_.clear();
    checkAndMap();
}


void Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    mapperPtr_.clear();
    checkAndMap();
}


void
Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchVectorField::autoMapGIB(mapper);
    mapperPtr_.clear();
    checkAndMap();
}


void Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
updateCoeffs()
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
    const scalarField ndmagS( (n & inletDir_) * patch().magSf() );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        forAll(inletDir_, fI)
        {
            if (ndmagS[fI]!=0)
            {
                refValue()[fI] = inletDir_[fI]*phip[fI]/ndmagS[fI];
            }
            else
            {
                refValue()[fI] = vector::zero;
            }
        }
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), rhoName_);

         const scalarField rhoNdS(rhop * ndmagS);

        forAll(inletDir_, fI)
        {
            if (rhoNdS[fI]!=0)
            {
                refValue()[fI] = inletDir_[fI]*phip[fI]/rhoNdS[fI];
            }
            else
            {
                refValue()[fI] = vector::zero;
            }
        }
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

    valueFraction() = 1.0 - pos0(phip);

    if (preserveSign_)
    {
        forAll(ndmagS, fI)
        {
            if ((ndmagS[fI]*phip[fI]) < 0)
            {
                refValue()[fI] = vector::zero;
                valueFraction()[fI] = 1.0;
            }
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    this->writeEntryIfDifferent
    (
        os,
        "fieldTable",
        this->internalField().name(),
        fieldTableName_
    );
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    this->writeEntryIfDifferent(os, "perturb", scalar(1e-5), perturb_);
    this->writeEntryIfDifferent(os, "delabella", Switch(true), delabella_);
    this->writeEntryIfDifferent
    (
        os, "cylindricalCoords", Switch(false), cylindricalCoords_
    );
    this->writeEntryIfDifferent
    (
        os, "preserveSign", Switch(false), preserveSign_
    );
    this->writeEntryIfDifferent
    (
        os,
        "mapMethod",
        word("planarInterpolation"),
        mapMethod_
    );
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureMappedDirectionInletOutletVelocityFvPatchVectorField::
operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=
    (
        valueFraction()*(inletDir_*(inletDir_ & pvf))
      + (1 - valueFraction())*pvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureMappedDirectionInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
