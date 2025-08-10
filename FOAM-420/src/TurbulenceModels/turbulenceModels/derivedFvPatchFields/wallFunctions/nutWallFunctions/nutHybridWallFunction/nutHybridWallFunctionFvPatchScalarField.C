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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "nutHybridWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
nutHybridWallFunctionFvPatchScalarField::calcNut() const
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
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU
    (
        mag(Uw.snGrad() - patch().nf() * (Uw.snGrad() & patch().nf()))
    );
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(SMALL),
        sqr(calcUTau(magGradU))/max(magGradU, SMALL) - nuw
    );
}


tmp<scalarField> nutHybridWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patch().index()];

    const fvPatchVectorField& Uw =
        turbModel.U().boundaryField()[patch().index()];
    const scalarField magUp
    (
        mag((Uw.patchInternalField() - Uw)
        - patch().nf()*((Uw.patchInternalField() - Uw) & patch().nf()))
    );

    const tmp<scalarField> tnuw = turbModel.nu(this->patch().index());
    const scalarField& nuw = tnuw();
    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    label nErrorFaces = 0;
    label maxIter = 100;
    scalar maxErr = 0.001;
    scalar sqrtVSMALL = sqrt(VSMALL);

    forAll(uTau, facei)
    {
        scalar magUpara = magUp[facei];

        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);

        if (ut > VSMALL)
        {
            int iter = 0;
            scalar err = GREAT;
            scalar yPlus = ut*y[facei]/nuw[facei];
            scalar uTauNew = 0;

            do
            {
                if (yPlus < 60)
                {
                    scalar kUu = min(kappa_*magUpara/ut, 50);
                    scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu)
                        - 1.0/6.0*kUu*sqr(kUu);

                    scalar f =
                        - ut*y[facei]/nuw[facei] + magUpara/ut
                        + 1.0/E_*(fkUu - 1.0/24.0*sqr(sqr(kUu)));

                    scalar df =
                        - y[facei]/nuw[facei] - magUpara/sqr(ut)
                        - 1.0/E_*kUu*fkUu/ut;

                    uTauNew = ut - f/(sign(df)*max(mag(df), VSMALL));

                    err = mag(f);

                    if (false)
                    {
                        Pout
                            << "bl patch: " << patch().index()
                            << ", facei: "<< facei
                            << ", iter: " << iter
                            << ", y+i " << yPlus
                            << ", y+ " << uTauNew*y[facei]/nuw[facei]
                            << ", U+ " << magUpara/uTauNew
                            << ", ut " << ut
                            << ", uTauNew " << uTauNew
                            << ", f " << f
                            << ", df " << df
                            << endl;

                    }
                }
                else
                {
                    scalar f
                        = magUpara/ut
                        - 1.0/kappa_*log(E_*ut*y[facei]/nuw[facei]);
                    scalar df = - magUpara/sqr(ut)
                        - 1.0/kappa_/ut;

                    uTauNew = ut - f/(sign(df)*max(mag(df), VSMALL));

                    err = mag(f);

                    if (false)
                    {
                        Pout
                            << "ll patch: " << patch().index()
                            << ", facei: "<< facei
                            << ", iter: " << iter
                            << ", y+i " << yPlus
                            << ", y+ " << uTauNew*y[facei]/nuw[facei]
                            << ", U+ " << magUpara/uTauNew
                            << ", ut " << ut
                            << ", uTauNew " << uTauNew
                            << ", f " << f
                            << ", df " << df
                            << endl;

                    }

                    //stop the solver from jumping too far toward zero if
                    //initial guess is too high
                    if (yPlus > yPlusMax_ && uTauNew*y[facei]/nuw[facei] < 60)
                    {
                        uTauNew = SMALL + 60*nuw[facei]/y[facei];
                    }

                }

                err = max(err, 2*mag(ut - uTauNew)/(mag(ut)+mag(uTauNew)));

                ut = uTauNew;

                yPlus = ut*y[facei]/nuw[facei];

                /*
                if
                (
                    patch().index() == 0 && facei == 311 && Pstream::myProcNo() == 1
                )
                {
                    //nErrorFaces++;

                    Pout
                       // << " patch: " << patch().index()
                       // << ", facei: "<< facei
                       // << ", iter: " << iter
                        << ", y+ " << yPlus
                        << ", ut " << ut
                        << ", err " << err
                        << endl;
                }*/

                if
                (
                    (iter >= maxIter -1 && err > maxErr && yPlus > 70)
                )
                {
                    nErrorFaces++;
                    Pout
                        << " patch: " << patch().index()
                        << ", facei: "<< facei
                        << ", iter: " << iter
                        << ", y+ " << yPlus
                        << ", ut " << ut
                        << ", err " << err
                        << endl;
                }

            } while (mag(ut) > sqrtVSMALL && err > maxErr && ++iter < maxIter);

            if (yPlus > yPlusMax_)
            {
                ut = magUpara*kappa_/log(E_*yPlusMax_);
            }
            uTau[facei] = max(0.0, ut);
        }
    }

    reduce(nErrorFaces, sumOp<label>());

    if (nErrorFaces > 0)
    {
        WarningInFunction
            << "Patch " << patch().name() << ": max iter " << maxIter
            << " reached before wall function "
            << "convergence for " << nErrorFaces << " faces." << endl;
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutHybridWallFunctionFvPatchScalarField::
nutHybridWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF)
{}


nutHybridWallFunctionFvPatchScalarField::
nutHybridWallFunctionFvPatchScalarField
(
    const nutHybridWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutHybridWallFunctionFvPatchScalarField::
nutHybridWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutHybridWallFunctionFvPatchScalarField::
nutHybridWallFunctionFvPatchScalarField
(
    const nutHybridWallFunctionFvPatchScalarField& wfpsf
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf)
{}


nutHybridWallFunctionFvPatchScalarField::
nutHybridWallFunctionFvPatchScalarField
(
    const nutHybridWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
nutHybridWallFunctionFvPatchScalarField::yPlus() const
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
    tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutHybridWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutHybridWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
