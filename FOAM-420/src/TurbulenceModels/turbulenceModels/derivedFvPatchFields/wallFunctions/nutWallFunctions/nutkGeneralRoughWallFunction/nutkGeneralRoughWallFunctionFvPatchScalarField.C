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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nutkGeneralRoughWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar nutkGeneralRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs,
    const scalar B,
    const scalar KsPlusRough,
    const scalar KsPlusSmooth
) const
{
    // Return fn based on non-dimensional roughness height

    if (KsPlus < KsPlusRough)
    {
        return pow
        (
            B*(KsPlus - KsPlusSmooth)/(KsPlusRough-KsPlusSmooth) + Cs*KsPlus,
            sin(3.1428/2.0*(log(KsPlus/KsPlusSmooth)/log(KsPlusRough/KsPlusSmooth)))
            //sin(0.4258*(log(KsPlus) - 0.811))
        );
    }
    else
    {
        return (B + Cs*KsPlus);
    }
}


tmp<scalarField> nutkGeneralRoughWallFunctionFvPatchScalarField::calcNut() const
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

    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/nuw[faceI];
        scalar KsPlus = uStar*Ks_[faceI]/nuw[faceI];

        scalar Edash = E_;

        if (KsPlus > KsPlusSmooth[faceI])
        {
            Edash /= fnRough
            (
                KsPlus,
                Cs_[faceI],
                B[faceI],
                KsPlusRough[faceI],
                KsPlusSmooth[faceI]
            );
        }

        if (yPlus > yPlusLam_)
        {
            scalar limitingNutw = max(nutw[faceI], nuw[faceI]);

            // To avoid oscillations limit the change in the wall viscosity
            // which is particularly important if it temporarily becomes zero
            nutw[faceI] =
                max
                (
                    min
                    (
                        nuw[faceI]
                       *(yPlus*kappa_/log(max(Edash*yPlus, 1+1e-4)) - 1),
                        2*limitingNutw
                    ), 0.5*limitingNutw
                );
        }

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[faceI]
                << endl;
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkGeneralRoughWallFunctionFvPatchScalarField::
nutkGeneralRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0),
    B(p.size(), 0.0),
    KsPlusRough(p.size(), 0.0),
    KsPlusSmooth(p.size(), 0.0)
{}


nutkGeneralRoughWallFunctionFvPatchScalarField::
nutkGeneralRoughWallFunctionFvPatchScalarField
(
    const nutkGeneralRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(mapper(ptf.Ks_)),
    Cs_(mapper(ptf.Cs_)),
    B(mapper(ptf.B)),
    KsPlusRough(mapper(ptf.KsPlusRough)),
    KsPlusSmooth(mapper(ptf.KsPlusSmooth))
{}


nutkGeneralRoughWallFunctionFvPatchScalarField::
nutkGeneralRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size()),
    B("B", dict, p.size()),
    KsPlusRough("KsPlusRough", dict, p.size()),
    KsPlusSmooth("KsPlusSmooth", dict, p.size())
{}


nutkGeneralRoughWallFunctionFvPatchScalarField::
nutkGeneralRoughWallFunctionFvPatchScalarField
(
    const nutkGeneralRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_),
    B(rwfpsf.B),
    KsPlusRough(rwfpsf.KsPlusRough),
    KsPlusSmooth(rwfpsf.KsPlusSmooth)
{}


nutkGeneralRoughWallFunctionFvPatchScalarField::
nutkGeneralRoughWallFunctionFvPatchScalarField
(
    const nutkGeneralRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_),
    B(rwfpsf.B),
    KsPlusRough(rwfpsf.KsPlusRough),
    KsPlusSmooth(rwfpsf.KsPlusSmooth)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutkGeneralRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    m(Ks_, Ks_);
    m(Cs_, Cs_);
    m(B, B);
    m(KsPlusSmooth, KsPlusSmooth);
    m(KsPlusRough, KsPlusRough);
}


void nutkGeneralRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutkGeneralRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutkGeneralRoughWallFunctionFvPatchScalarField>(ptf);

    Ks_.rmap(nrwfpsf.Ks_, addr);
    Cs_.rmap(nrwfpsf.Cs_, addr);
    B.rmap(nrwfpsf.B, addr);
    KsPlusSmooth.rmap(nrwfpsf.KsPlusSmooth, addr);
    KsPlusRough.rmap(nrwfpsf.KsPlusRough, addr);
}


void nutkGeneralRoughWallFunctionFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    nutkWallFunctionFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(Ks_, scalar(0));
    mapper.map(Cs_, scalar(0));
    mapper.map(B, scalar(0));
    mapper.map(KsPlusSmooth, scalar(0));
    mapper.map(KsPlusRough, scalar(0));
}


void nutkGeneralRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    B.writeEntry("B", os);
    KsPlusRough.writeEntry("KsPlusRough", os);
    KsPlusSmooth.writeEntry("KsPlusSmooth", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(fvPatchScalarField, nutkGeneralRoughWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
