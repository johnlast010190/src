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
    (c) 2017-2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nutUSpaldingWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool nutUSpaldingWallFunctionFvPatchScalarField::isUpToDate() const
{
    scalar currentTime = db().time().value();
    scalar deltaT = db().time().deltaT().value();

    if
    (
        lastUpdateTimeValue_ < currentTime + 0.5*deltaT
     && lastUpdateTimeValue_ > currentTime - 0.5*deltaT
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

tmp<scalarField> nutUSpaldingWallFunctionFvPatchScalarField::calcNut() const
{
    if (isUpToDate())
    {
        return
        (
            tmp<scalarField>(new scalarField(*this))
        );
    }

    // get ref to circumvent the const
    scalar& lastTime = const_cast<scalar&>(lastUpdateTimeValue_);
    lastTime = db().time().value();

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


tmp<scalarField> nutUSpaldingWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.snGrad()/this->patch().deltaCoeffs()));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    forAll(uTau, facei)
    {
        if (nutw[facei] < 0 || nuw[facei] < 0 || magGradU[facei] < 0)
        {
            FatalErrorInFunction
                << "Negative input value for shear velocity calculation:" << nl
                << "    nuw[facei]: " << nuw[facei] << nl
                << "    nutw[facei]: " << nutw[facei] << nl
                << "    magGradU[facei]: " << magGradU[facei]
                << exit(FatalError);
        }

        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[facei]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[facei]/nuw[facei]
                    + magUp[facei]/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[facei]/nuw[facei]
                  + magUp[facei]/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;
            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau[facei] = max(0.0, ut);
        }
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    lastUpdateTimeValue_
    (
        db().time().value() - db().time().deltaT().value()
    )
{}


nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const nutUSpaldingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    lastUpdateTimeValue_
    (
        ptf.getLastUpdatedTime()
    )
{}


nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    lastUpdateTimeValue_
    (
        dict.lookupOrDefault<scalar>
        (
            "lastUpdated",
            db().time().value() - db().time().deltaT().value()
        )
    )
{}


nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const nutUSpaldingWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    lastUpdateTimeValue_
    (
        wfpsf.getLastUpdatedTime()
    )
{}


nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const nutUSpaldingWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    lastUpdateTimeValue_
    (
        wfpsf.getLastUpdatedTime()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUSpaldingWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutUSpaldingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeEntry("lastUpdated", lastUpdateTimeValue_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUSpaldingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
