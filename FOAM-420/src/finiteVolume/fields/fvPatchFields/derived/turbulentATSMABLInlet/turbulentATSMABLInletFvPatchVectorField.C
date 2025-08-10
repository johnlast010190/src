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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/turbulentATSMABLInlet/turbulentATSMABLInletFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// namespace incompressible
// {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void turbulentATSMABLInletFvPatchVectorField::initialiseBaseClassData
(
    const fvPatch& p,
    const dictionary& dict
)
{
    vectorField::operator=
    (
        vectorField("value", dict, p.size())
    );

    density_ = dict.lookupOrDefault<scalar>("density", 1.0);
    perturb_ = dict.lookupOrDefault<scalar>("perturb", 1e-5);
    mapMethod_ = dict.lookupOrDefault<word>("mapMethod", "nearestCell");

    interpolateR_ = false;
    interpolateL_ = false;
    interpolateU_ = false;

    periodicInY_ = dict.lookupOrDefault<bool>("periodicInY", false);
    periodicInZ_ = dict.lookupOrDefault<bool>("periodicInZ", false);

    vortonType_  = dict.lookupOrDefault<word>("vortonType", "typeR");
    nVortonGlobal_ = dict.lookupOrDefault<label>("nVorton", 0);
    nVortonLocal_ = dict.lookupOrDefault<label>("nVortonLocal", 0);

    isCleanRestart_ = dict.lookupOrDefault<bool>("cleanRestart", false);

    nOutputFace_ = dict.lookupOrDefault<label>("nOutputFace", 0);
    outputFaceIndices_ = dict.lookupOrDefault<labelList>("outputFaceIndices", labelList(nOutputFace_, 0));

    //set U,R,L fields in ATSM base class
    profile_->setProfiles(U_, L_, R_);

    // Set UMean as patch area average value
    scalarField UMag(mag(U_));
    UMean_ = gSum(UMag*patch().magSf())/(gSum(patch().magSf()) + ROOTVSMALL);
    UMax_ = gMax(UMag);

    if (nVortonLocal_ && !isCleanRestart_)
    {
        isRestart_ = true;

        ITstream& is = dict.lookup("vortonLabel");
        is >> static_cast<List<label>&>(vortonLabel_);

        is = dict.lookup("vortonPosition");
        is >> static_cast<List<vector>&>(vortonPosition_);

        is = dict.lookup("vortonDistance");
        is >> static_cast<List<scalar>&>(vortonDistance_);

        is = dict.lookup("vortonScale");
        is >> static_cast<List<vector>&>(vortonScale_);

        is = dict.lookup("vortonIntensity");
        is >> static_cast<List<vector>&>(vortonIntensity_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentATSMABLInletFvPatchVectorField::
turbulentATSMABLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    turbulentATSMInletFvPatchVectorField(p, iF),
    profile_()
{}

turbulentATSMABLInletFvPatchVectorField::
turbulentATSMABLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentATSMInletFvPatchVectorField(p, iF),
    profile_(ABLProfile::New(p, dict))
{
    initialiseBaseClassData(p, dict);
}


turbulentATSMABLInletFvPatchVectorField::
turbulentATSMABLInletFvPatchVectorField
(
    const turbulentATSMABLInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentATSMInletFvPatchVectorField(ptf, p, iF, mapper),
    profile_(ptf.profile_->clone())
{}


turbulentATSMABLInletFvPatchVectorField::
turbulentATSMABLInletFvPatchVectorField
(
    const turbulentATSMABLInletFvPatchVectorField& ptf
)
:
    turbulentATSMInletFvPatchVectorField(ptf),
    profile_(ptf.profile_->clone())
{}

turbulentATSMABLInletFvPatchVectorField::
turbulentATSMABLInletFvPatchVectorField
(
    const turbulentATSMABLInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    turbulentATSMInletFvPatchVectorField(ptf, iF),
    profile_(ptf.profile_->clone())
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentATSMABLInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    turbulentATSMInletFvPatchVectorField::updateCoeffs();
}



void turbulentATSMABLInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);

    profile_->write(os);

    writeEntryIfDifferent<bool>(os, "periodicInY", false, periodicInY_);
    writeEntryIfDifferent<bool>(os, "periodicInZ", false, periodicInZ_);

    U_.writeEntry("U", os);
    R_.writeEntry("R", os);
    L_.writeEntry("Lxyz", os);

    writeEntryIfDifferent<scalar>(os, "density", 1.0, density_);
    writeEntryIfDifferent<scalar>(os, "perturb", 1e-5, perturb_);

    writeEntryIfDifferent<label>(os, "nVorton", 0, nVortonGlobal_);
    writeEntryIfDifferent<label>(os, "nVortonLocal", 0, nVortonLocal_);
    writeEntryIfDifferent<word>(os, "vortonType", "typeR", vortonType_);

    if (nVortonLocal_)
    {
        vortonLabel_.writeEntry("vortonLabel", os);
        vortonPosition_.writeEntry("vortonPosition", os);
        vortonDistance_.writeEntry("vortonDistance", os);
        vortonScale_.writeEntry("vortonScale", os);
        vortonIntensity_.writeEntry("vortonIntensity", os);
    }

    if (!mapMethod_.empty())
    {
        writeEntryIfDifferent<word>
        (
            os,
            "mapMethod",
            "nearestCell",
            mapMethod_
        );
    }

    if (nOutputFace_ > 0)
    {
        os.writeEntry("nOutputFace", nOutputFace_);
        outputFaceIndices_.writeEntry("outputFaceIndices", os);
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    turbulentATSMABLInletFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
