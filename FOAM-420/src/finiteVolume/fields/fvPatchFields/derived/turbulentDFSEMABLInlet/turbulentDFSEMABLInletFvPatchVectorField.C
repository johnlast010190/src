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

#include "fields/fvPatchFields/derived/turbulentDFSEMABLInlet/turbulentDFSEMABLInletFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void turbulentDFSEMABLInletFvPatchVectorField::initialiseBaseClassData
(
    const fvPatch& p,
    const dictionary& dict
)
{
    vectorField::operator=
    (
        vectorField("value", dict, p.size())
    );

    delta_ = readScalar(dict.lookup("delta"));

    d_ = dict.lookupOrDefault<scalar>("d", 1);

    kappa_ = dict.lookupOrDefault<scalar>("kappa", 0.41);

    perturb_ = dict.lookupOrDefault<scalar>("perturb", 1e-5);

    mapMethod_ = dict.lookupOrDefault<word>("mapMethod", "planarInterpolation");

    nCellPerEddy_ = dict.lookupOrDefault<label>("nCellPerEddy", 5);
    writeEddies_ = dict.lookupOrDefault<bool>("writeEddies", false);

    //set U,R,L fields in DFSEM base class
    profile_->setProfiles(U_, L_, R_);

    UMean_ = gSum(U_*patch().magSf())/gSum(patch().magSf());

    eddy::debug = debug;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentDFSEMABLInletFvPatchVectorField::
turbulentDFSEMABLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    turbulentDFSEMInletFvPatchVectorField(p, iF),
    profile_()
{}


turbulentDFSEMABLInletFvPatchVectorField::
turbulentDFSEMABLInletFvPatchVectorField
(
    const turbulentDFSEMABLInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentDFSEMInletFvPatchVectorField(ptf, p, iF, mapper),
    profile_(ptf.profile_->clone())
{}


turbulentDFSEMABLInletFvPatchVectorField::
turbulentDFSEMABLInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentDFSEMInletFvPatchVectorField(p, iF),
    profile_(ABLProfile::New(p, dict))
{
    initialiseBaseClassData(p, dict);
}


turbulentDFSEMABLInletFvPatchVectorField::
turbulentDFSEMABLInletFvPatchVectorField
(
    const turbulentDFSEMABLInletFvPatchVectorField& blpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    turbulentDFSEMInletFvPatchVectorField(blpvf, iF),
    profile_(blpvf.profile_->clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentDFSEMABLInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    turbulentDFSEMInletFvPatchVectorField::updateCoeffs();
}



void turbulentDFSEMABLInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeEntry("delta", delta_);
    writeEntryIfDifferent(os, "d", scalar(1.0), d_);
    writeEntryIfDifferent(os, "kappa", scalar(0.41), kappa_);
    writeEntryIfDifferent(os, "perturb", scalar(1e-5), perturb_);
    writeEntryIfDifferent(os, "nCellPerEddy", label(5), nCellPerEddy_);
    writeEntryIfDifferent(os, "writeEddies", false, writeEddies_);

    profile_->write(os);

    R_.writeEntry("R", os);
    L_.writeEntry("L", os);
    U_.writeEntry("U", os);

    if
    (
       !mapMethod_.empty()
     && mapMethod_ != "planarInterpolation"
    )
    {
        os.writeEntry("mapMethod", mapMethod_);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    turbulentDFSEMABLInletFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
