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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nutkDemirelRoughWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar nutkDemirelRoughWallFunctionFvPatchScalarField::E
(
    const scalar KsPlus
) const
{
    if (KsPlus < 3.85)
    {
        return E_;
    }
    else if ((3.85 < KsPlus) && (KsPlus < 15.0))
    {
        return
            E_
            /exp
            (pow
            (
                log(0.26*KsPlus),
                sin(0.5*constant::mathematical::pi*log10(KsPlus/3.0)*1.43067655)
            ));
    }
    else
    {
        return E_
               /exp
               (log(0.26*KsPlus));
    }
}


tmp<scalarField> nutkDemirelRoughWallFunctionFvPatchScalarField::calcNut() const
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

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar uStar = Cmu25*sqrt(k[celli]);
        scalar Ks = roughWallCoeffs().roughnessHeight();
        scalar KsPlus = uStar*Ks/nuw[facei];
        scalar E = this->E(KsPlus);
        const scalar yPlusMin = constant::mathematical::e/E;
        const scalar yPlus = max(uStar*y[facei]/nuw[facei], yPlusMin);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        nutw[facei] =
            max
            (
                min
                (
                    nuw[facei]*max(yPlus*kappa_/log(E*yPlus) - 1, 0),
                    max(2*nutw[facei], nuw[facei])
                ),
                0.5*nutw[facei]
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", E = " << E
                << ", nutw = " << nutw[facei]
                << endl;
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkDemirelRoughWallFunctionFvPatchScalarField::nutkDemirelRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF)
{
    checkIfValidRoughnessCoeffs();
}


nutkDemirelRoughWallFunctionFvPatchScalarField::nutkDemirelRoughWallFunctionFvPatchScalarField
(
    const nutkDemirelRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{
    checkIfValidRoughnessCoeffs();
}


nutkDemirelRoughWallFunctionFvPatchScalarField::nutkDemirelRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict)
{
    checkIfValidRoughnessCoeffs();
}


nutkDemirelRoughWallFunctionFvPatchScalarField::nutkDemirelRoughWallFunctionFvPatchScalarField
(
    const nutkDemirelRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf)
{
    checkIfValidRoughnessCoeffs();
}


nutkDemirelRoughWallFunctionFvPatchScalarField::nutkDemirelRoughWallFunctionFvPatchScalarField
(
    const nutkDemirelRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF)
{
    checkIfValidRoughnessCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar nutkDemirelRoughWallFunctionFvPatchScalarField::getUPlus(const scalar& yPlus) const
{
    scalar UPlus = 0;
    return UPlus;
};

void nutkDemirelRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
}


void nutkDemirelRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);
}


void nutkDemirelRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkDemirelRoughWallFunctionFvPatchScalarField
)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
