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

#include "derivedFvPatchFields/pressureVelocity/pressureVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "derivedFvPatchFields/basePressure/basePressureFvPatchScalarField.H"
#include "fields/Fields/symmTransformField/symmTransformField.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char*
NamedEnum
<
    pressureVelocityFvPatchVectorField::directionSpecification,
    4
>::names[] =
{
    "userSpecified",
    "zeroGradient",
    "patchNormal",
    "tangentialVelocity"
};

const NamedEnum
<
    pressureVelocityFvPatchVectorField::directionSpecification,
    4
> pressureVelocityFvPatchVectorField::directionSpecificationNames_;

}


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

bool Foam::pressureVelocityFvPatchVectorField::isCoupledSolver() const
{
    if (this->db().foundObject<foamCoupledControl>(solutionControl::typeName))
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::pressureVelocityFvPatchVectorField::setDirection()
{
    vectorField localCfs(coorFramePtr_->coorSys().localPosition(patch().Cf()));

    localCfs += velDir_;
    velDirField_ =
        coorFramePtr_->coorSys().globalPosition(localCfs)
      - patch().Cf();
}


void Foam::pressureVelocityFvPatchVectorField::setTangentialFraction
(
    const scalarField& phip
)
{
    valueFraction() = (I - sqr(patch().nf()))*neg(phip);

    // TODO Tangential velocity contribution has to be re-checked
    // if (directionType_ != zeroGradient && !isCoupledSolver())
    // {
    //     const label patchI = this->patch().index();
    //     const fvPatchScalarField& pp =
    //         db().lookupObject<volScalarField>(pName_).boundaryField()[patchI];
    //     if (isA<basePressureFvPatchScalarField>(pp))
    //     {
    //         const basePressureFvPatchScalarField& ppres =
    //             dynamic_cast<const basePressureFvPatchScalarField&>(pp);

    //         const scalarField& rhop =
    //             db().lookupObject<volScalarField>
    //             (
    //                 ppres.rhoName()
    //             ).boundaryField()[patchI];
    //         valueFraction() +=
    //             2.0*(ppres.C3Field()/rhop)*pos0(phip)*(I - sqr(patch().nf()));
    //         // Resistive pressure used to be:
    //         // valueFraction() +=
    //         //     ppres.Ct()*pos0(phip)*(I - sqr(patch().nf()));
    //     }
    // }
}


void
Foam::pressureVelocityFvPatchVectorField::transformAndSetTangentialVelocity()
{
    vectorField transformedField(velDirField_);

    tmp<vectorField> tlocCf
    (
        coorFramePtr_->coorSys().localPosition(patch().Cf())
    );
    vectorField& locCf = tlocCf.ref();

    locCf += transformedField;
    transformedField =
        coorFramePtr_->coorSys().globalPosition(locCf)
      - patch().Cf();

    const vectorField n(patch().nf());
    refValue() = transformedField - n*(n & transformedField);
}


Foam::tmp<Foam::vectorField>
Foam::pressureVelocityFvPatchVectorField::implicitCoeff
(
    const basePressureFvPatchScalarField& ppres
) const
{
    const vectorField nf(patch().nf());
    const vectorField& Up = *this;
    const vectorField Un(nf*(Up & nf));
    const vectorField Ut(Up - Un);

    return ppres.C1Field()*nf + ppres.C2Field()*Un + ppres.C3Field()*Ut;
}


void Foam::pressureVelocityFvPatchVectorField::reportBackFlow
(
    const fvsPatchField<scalar>& phip
) const
{
    scalar backFlowArea = 0.0;
    if
    (
        patch().type() == "outlet"
     || patch().patch().physicalType() == "outlet"
    )
    {
        const scalarField& magSf = patch().magSf();
        forAll(phip, facei)
        {
            if (phip[facei] < 0)
            {
                backFlowArea += magSf[facei];
            }
        }
        reduce(backFlowArea, plusOp<scalar>());
    }
    scalar areaPer = backFlowArea/totalArea_;
    if (areaPer > 1e-5) // report backflow for > 0.001%
    {
        Info<< "BackFlow of " << areaPer << "% "
            << "at patch " << this->patch().name()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureVelocityFvPatchVectorField::pressureVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    pName_("p"),
    directionType_(patchNormal),
    velDir_(Zero),
    velDirField_(this->size(), Zero),
    coorFramePtr_(nullptr),
    backFlowReport_(false),
    totalArea_(Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symmTensor::zero;
}


Foam::pressureVelocityFvPatchVectorField::pressureVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    directionType_
    (
        directionSpecificationNames_.readOrDefault
        (
            dict.lookupOrDefault<word>
            (
                "directionSpecification",
                "patchNormal"
            ),
            patchNormal
        )
    ),
    velDir_(Zero),
    velDirField_(this->size(), Zero),
    coorFramePtr_(nullptr),
    backFlowReport_(dict.lookupOrDefault<Switch>("backFlowReport", false)),
    totalArea_(Zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symmTensor::zero;

    if (directionType_ == userSpecified)
    {
        if (dict.found("velocityDirection"))
        {
            velDir_ = dict.lookupOrDefault<vector>("velocityDirection", Zero);
            velDir_.normalise();
        }
        else
        {
            FatalErrorInFunction
                << "velocityDirection is not specified for userSpecified "
                << "direction."
                << exit(FatalError);
        }

        if (dict.found("referenceFrame"))
        {
            coorFramePtr_ =
                coordinateFrame::lookupNew(this->internalField().mesh(), dict);
            setDirection();
        }
        else
        {
            //- velDir == absolute
            velDirField_ = velDir_;
        }
    }
    else if (directionType_ == tangentialVelocity)
    {
        if (dict.found("tangentialVelocity"))
        {
            velDirField_ =
                vectorField("tangentialVelocity", dict, p.size());

            if (dict.found("referenceFrame"))
            {
                coorFramePtr_ =
                    coordinateFrame::lookupNew
                    (
                        this->internalField().mesh(),
                        dict
                    );
                transformAndSetTangentialVelocity();
            }
            else
            {
                //- velDir == absolute
                refValue() = velDirField_;
            }
        }
        else
        {
            refValue() = Zero;
        }
    }
    if (backFlowReport_)
    {
        totalArea_ = gSum(patch().magSf());
    }
}


Foam::pressureVelocityFvPatchVectorField::pressureVelocityFvPatchVectorField
(
    const pressureVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    directionType_(ptf.directionType_),
    velDir_(ptf.velDir_),
    velDirField_(mapper(ptf.velDirField_)),
    coorFramePtr_(ptf.coorFramePtr_),
    backFlowReport_(ptf.backFlowReport_),
    totalArea_(ptf.totalArea_)
{
}


Foam::pressureVelocityFvPatchVectorField::pressureVelocityFvPatchVectorField
(
    const pressureVelocityFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    directionType_(ptf.directionType_),
    velDir_(ptf.velDir_),
    velDirField_(ptf.velDirField_),
    coorFramePtr_(ptf.coorFramePtr_),
    backFlowReport_(ptf.backFlowReport_),
    totalArea_(ptf.totalArea_)
{
}


Foam::pressureVelocityFvPatchVectorField::pressureVelocityFvPatchVectorField
(
    const pressureVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    directionType_(ptf.directionType_),
    velDir_(ptf.velDir_),
    velDirField_(ptf.velDirField_),
    coorFramePtr_(ptf.coorFramePtr_),
    backFlowReport_(ptf.backFlowReport_),
    totalArea_(ptf.totalArea_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
    m(velDirField_, velDirField_);
}


void Foam::pressureVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
    const pressureVelocityFvPatchVectorField& dmptf =
        refCast<const pressureVelocityFvPatchVectorField>(ptf);
    velDirField_.rmap(dmptf.velDirField_, addr);
}


void Foam::pressureVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    directionMixedFvPatchVectorField::autoMapGIB(mapper);
    mapper.map(velDirField_, vector::zero);
}


void Foam::pressureVelocityFvPatchVectorField::addMomentumGradPCoupledBC
(
    fvBlockMatrix<vector>& bUEq,
    const volScalarField& p
)
{
    if (isCoupledSolver())
    {
        const label patchi = this->patch().index();
        if (isA<basePressureFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            const basePressureFvPatchScalarField& ppres =
                dynamic_cast<const basePressureFvPatchScalarField&>
                (
                    p.boundaryField()[patchi]
                );
          // Switch off implicit treatment for the MRF for now
          // It has to be investigated how to properly treat the MRF
          if (!ppres.coorFramePtr())
          {
            const fvMesh& mesh = internalField().mesh();
            const vectorField& pSf = mesh.boundary()[patchi].Sf();
            Field<tensor>& diag = bUEq.diag().asSquare();
            Field<vector>& source = bUEq.source();

            tmp<vectorField> tcoeff(implicitCoeff(ppres));
            tmp<scalarField> tp0(ppres.C0Field());
            const scalarField& p0 = tp0();
            const vectorField& coeff = tcoeff();

            const labelList& faceCells = mesh.boundary()[patchi].faceCells();
            const vectorField Up(this->patchInternalField()());

            forAll(faceCells, facei)
            {
                const label faceCelli = faceCells[facei];
                const vector& coeffI = coeff[facei];
                /*
                // Consistent
                diag[faceCelli].xx() += coeffI.x()*pSf[facei].x();
                diag[faceCelli].xy() += coeffI.y()*pSf[facei].x();
                diag[faceCelli].xz() += coeffI.z()*pSf[facei].x();
                diag[faceCelli].yx() += coeffI.x()*pSf[facei].y();
                diag[faceCelli].yy() += coeffI.y()*pSf[facei].y();
                diag[faceCelli].yz() += coeffI.z()*pSf[facei].y();
                diag[faceCelli].zx() += coeffI.x()*pSf[facei].z();
                diag[faceCelli].zy() += coeffI.y()*pSf[facei].z();
                diag[faceCelli].zz() += coeffI.z()*pSf[facei].z();
                */
                tensor tCoeff = pSf[facei]*coeffI;
                if (sign(tCoeff.xx()) == sign(diag[faceCelli].xx()))
                {
                    diag[faceCelli].xx() += tCoeff.xx();
                }
                else
                {
                    source[faceCelli].x() -= tCoeff.xx()*Up[facei].x();
                }
                if (sign(tCoeff.xy()) == sign(diag[faceCelli].xy()))
                {
                    diag[faceCelli].xy() += tCoeff.xy();
                }
                else
                {
                    source[faceCelli].x() -= tCoeff.xy()*Up[facei].y();
                }
                if (sign(tCoeff.xz()) == sign(diag[faceCelli].xz()))
                {
                    diag[faceCelli].xz() += tCoeff.xz();
                }
                else
                {
                    source[faceCelli].x() -= tCoeff.xz()*Up[facei].z();
                }
                if (sign(tCoeff.yx()) == sign(diag[faceCelli].yx()))
                {
                    diag[faceCelli].yx() += tCoeff.yx();
                }
                else
                {
                    source[faceCelli].y() -= tCoeff.yx()*Up[facei].x();
                }
                if (sign(tCoeff.yy()) == sign(diag[faceCelli].yy()))
                {
                    diag[faceCelli].yy() += tCoeff.yy();
                }
                else
                {
                    source[faceCelli].y() -= tCoeff.yy()*Up[facei].y();
                }
                if (sign(tCoeff.yz()) == sign(diag[faceCelli].yz()))
                {
                    diag[faceCelli].yz() += tCoeff.yz();
                }
                else
                {
                    source[faceCelli].y() -= tCoeff.yz()*Up[facei].z();
                }
                if (sign(tCoeff.zx()) == sign(diag[faceCelli].zx()))
                {
                    diag[faceCelli].zx() += tCoeff.zx();
                }
                else
                {
                    source[faceCelli].z() -= tCoeff.zx()*Up[facei].x();
                }
                if (sign(tCoeff.zy()) == sign(diag[faceCelli].zy()))
                {
                    diag[faceCelli].zy() += tCoeff.zy();
                }
                else
                {
                    source[faceCelli].z() -= tCoeff.zy()*Up[facei].y();
                }
                if (sign(tCoeff.zz()) == sign(diag[faceCelli].zz()))
                {
                    diag[faceCelli].zz() += tCoeff.zz();
                }
                else
                {
                    source[faceCelli].z() -= tCoeff.zz()*Up[facei].z();
                }

                source[faceCelli] -= p0[facei]*pSf[facei];
            }
        }
      }
    }
}


void Foam::pressureVelocityFvPatchVectorField::addContinuityCoupledBC
(
    BlockLduSystem<vector, scalar>& bUEq,
    const volScalarField& p,
    const surfaceTensorField& rDf
) const
{
    const label patchi = this->patch().index();
    if (isA<basePressureFvPatchScalarField>(p.boundaryField()[patchi]))
    {
        const basePressureFvPatchScalarField& ppres =
            dynamic_cast<const basePressureFvPatchScalarField&>
            (
                p.boundaryField()[patchi]
            );
      // Switch off implicit treatment for the MRF for now
      // It has to be investigated how to properly treat the MRF
      if (!ppres.coorFramePtr())
      {
        const fvMesh& mesh = internalField().mesh();

        const vectorField& pSf = mesh.boundary()[patchi].Sf();
        const scalarField& dCoeff = this->patch().deltaCoeffs();
        const tensorField& rDfp = rDf.boundaryField()[patchi];
        tmp<vectorField> tnf(this->patch().nf());
        const vectorField& nf = tnf();

        vectorField& diag = bUEq.diag().asLinear();
        scalarField& source = bUEq.source();

        typename CoeffField<vector>::linearTypeField& blockI =
            bUEq.internalCoeffs()[patchi].asLinear();
        scalarField& blockB = bUEq.boundaryCoeffs()[patchi];

        const vectorField Up(this->patchInternalField()());
        tmp<vectorField> tcoeff(implicitCoeff(ppres));
        tmp<scalarField> tp0(ppres.C0Field());
        const scalarField& p0 = tp0();
        const vectorField& coeff = tcoeff();

        const labelList& faceCells = mesh.boundary()[patchi].faceCells();
        forAll(faceCells, facei)
        {
            const label faceCelli = faceCells[facei];

            const vector& coeffI = coeff[facei];
            const vector diffCoeff = (rDfp[facei]&pSf[facei])*dCoeff[facei];
            const scalar sdiffCoeff = diffCoeff&nf[facei];
            const vector coeffTerm = -(coeffI)*sdiffCoeff;
            const scalar sourceCoeff = -p0[facei]*sdiffCoeff;

            /*
            // Consistent
            diag[faceCelli].x() += coeffTerm.x();
            blockI[facei].x() += coeffTerm.x();
            diag[faceCelli].y() += coeffTerm.y();
            blockI[facei].y() += coeffTerm.y();
            diag[faceCelli].z() += coeffTerm.z();
            blockI[facei].z() += coeffTerm.z();
            */
            if (sign(coeffTerm.x()) == sign(diag[faceCelli].x()))
            {
                diag[faceCelli].x() += coeffTerm.x();
                blockI[facei].x() += coeffTerm.x();
            }
            else
            {
                scalar explicitTerm = -coeffTerm.x()*Up[facei].x();
                source[faceCelli] += explicitTerm;
                blockB[facei] += explicitTerm;
            }
            if (sign(coeffTerm.y()) == sign(diag[faceCelli].y()))
            {
                diag[faceCelli].y() += coeffTerm.y();
                blockI[facei].y() += coeffTerm.y();
            }
            else
            {
                scalar explicitTerm = -coeffTerm.y()*Up[facei].y();
                source[faceCelli] += explicitTerm;
                blockB[facei] += explicitTerm;
            }
            if (sign(coeffTerm.z()) == sign(diag[faceCelli].z()))
            {
                diag[faceCelli].z() += coeffTerm.z();
                blockI[facei].z() += coeffTerm.z();
            }
            else
            {
                scalar explicitTerm = -coeffTerm.z()*Up[facei].z();
                source[faceCelli] += explicitTerm;
                blockB[facei] += explicitTerm;
            }

            source[faceCelli] -= sourceCoeff;
            blockB[facei] -= sourceCoeff;
        }
      }
    }
}


void Foam::pressureVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            db(),
            phiName_
        );
    if (backFlowReport_ && totalArea_ > SMALL)
    {
        reportBackFlow(phip);
    }

    switch (directionType_)
    {
        case userSpecified:
        {
            vectorField nf(patch().nf());
            vectorField Un((*this&nf)*nf);
            vectorField Unt((Un&velDirField_)*velDirField_);
            scalarField magUn(mag(Un));
            scalarField magUnt(mag(Unt));

            scalarField magVector(magUn*magUn);
            forAll(magVector, facei)
            {
                if (magUnt[facei] > 0)
                {
                    magVector[facei] /= magUnt[facei];
                }
                else
                {
                    magVector[facei] = 0;
                }
            }
            refValue() = magVector*velDirField_;
            setTangentialFraction(phip);
            break;
        }
        case zeroGradient:
        {
            valueFraction() = Zero;
            break;
        }
        case patchNormal:
        {
            setTangentialFraction(phip);
            break;
        }

        case tangentialVelocity:
        {
            setTangentialFraction(phip);
            break;
        }
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
    //- reset update variable.
    directionMixedFvPatchVectorField::updateCoeffs();
}


void Foam::pressureVelocityFvPatchVectorField::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>& bEq
) const
{
    const fvPatchField<scalar>& ptf =
        patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), pName_);
    if
    (
        isA<basePressureFvPatchScalarField>(ptf)
    )
    {
        const basePressureFvPatchScalarField& bpp =
            refCast<const basePressureFvPatchScalarField>(ptf);
        bpp.boundaryRelaxMatrix(bEq);
    }
}


void Foam::pressureVelocityFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    os.writeEntry
    (
        "directionSpecification",
        directionSpecificationNames_[directionType_]
    );

    if (directionType_ == userSpecified)
    {
        writeEntryIfDifferent<vector>(os, "velocityDirection", Zero, velDir_);
        if (coorFramePtr_)
        {
            os.writeEntry("referenceFrame", coorFramePtr_->name());
        }
    }
    if (directionType_ == tangentialVelocity)
    {
        if (velDirField_.size())
        {
            velDirField_.writeEntry("tangentialVelocity", os);
            if (coorFramePtr_)
            {
                os.writeEntry("referenceFrame", coorFramePtr_->name());
            }
        }
    }
}


// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::pressureVelocityFvPatchVectorField::valueBoundaryCoeffs
// (
//     const tmp<scalarField>&
// ) const
// {
//     const fvsPatchField<scalar>& phip =
//         patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
//         (
//             db(),
//             phiName_
//         );
//     return
//         (
//             neg(phip)
//            *(
//                 *this
//               - cmptMultiply
//                 (
//                     valueInternalCoeffs(this->patch().weights()),
//                     this->patchInternalField()
//                 )
//             )
//         );
// }


void Foam::pressureVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    tmp<vectorField> normalValue(transform(valueFraction(), refValue()));
    tmp<vectorField> transformGradValue(transform(I - valueFraction(), pvf));
    fvPatchField<vector>::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
