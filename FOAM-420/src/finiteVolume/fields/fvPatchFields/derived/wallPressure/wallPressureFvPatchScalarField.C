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
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/wallPressure/wallPressureFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "memory/autoPtr/autoPtr.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * * Namespace Foam  * * * * * * * * * * * * * * * * //
namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

label wallPressureFvPatchScalarField::globInterpIndex = 0;
autoPtr<surfaceScalarField> wallPressureFvPatchScalarField::facePressure(nullptr);


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

void wallPressureFvPatchScalarField::calculateMultiFaceAdressing()
{
    multiFaceAddrPtr_.reset
    (
        new List<List<Tuple2<label, label>>>(this->size())
    );
    wallPressureFaceAddrPtr_.reset
    (
        new List<List<Tuple2<label, label>>>(this->size())
    );

    // reference to mesh
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    // reference parent field
    const volScalarField& p
    (
        this->db().lookupObject<volScalarField>
        (
            this->internalField().name()
        )
    );

    // check if faceCell has any additional boundary faces with type
    // wallPressure other than the one connecting to the cell
    // if so add them to a list

    const labelList& fcs = this->patch().faceCells();

    forAll(fcs, pfI)
    {
        const labelList& fCells = mesh.cells()[fcs[pfI]];
        label cfgi = pfI + this->patch().start();

        DynamicList<Tuple2<label, label>> addwPF(fCells.size());
        DynamicList<Tuple2<label, label>> addoPF(fCells.size());

        forAll(fCells, fci)
        {
            label cfaceI = fCells[fci];
            if (cfaceI != cfgi && cfaceI >= mesh.nInternalFaces())
            {
                label wP = mesh.boundaryMesh().whichPatch(cfaceI);

                if (!isA<emptyFvPatch>(mesh.boundary()[wP]))
                {

                    if
                    (
                        isA<wallPressureFvPatchScalarField>
                        (p.boundaryField()[wP])
                    )
                    {
                        addwPF.append
                        (
                            Tuple2<label,label>
                            (
                                wP, cfaceI - mesh.boundaryMesh()[wP].start()
                            )
                        );
                    }

                    addoPF.append
                    (
                        Tuple2<label,label>
                        (
                            wP, cfaceI - mesh.boundaryMesh()[wP].start()
                        )
                    );
                }
            }
        }


        addwPF.shrink();
        addoPF.shrink();

        multiFaceAddrPtr_->operator[](pfI) = addoPF;
        wallPressureFaceAddrPtr_->operator[](pfI) = addwPF;
    }
}


scalar wallPressureFvPatchScalarField::faceNormalGradient
(
    const fvPatchScalarField& cPatch,
    const label faceI,
    const vector& nf,
    const List<Tuple2<label, label>>& otherCellBoundaryFaces,
    label cellI
)
{
    // reference to mesh
    const fvMesh& mesh(cPatch.patch().boundaryMesh().mesh());

    // reference to field this belongs to
    const volScalarField& p
    (
        this->db().lookupObject<volScalarField>
        (
            this->internalField().name()
        )
    );

    // implicitly calculate the surface normal gradient in the faceCell
    // assuming boundary snGrad is the same as the cell center gradient

    const labelList& fCells = mesh.cells()[cellI];
    scalar gradSn = 0;

    //
    label cfgI = cPatch.patch().start() + faceI;
    const surfaceScalarField& ps
    (
        wallPressureFvPatchScalarField::facePressure()
    );
    const surfaceVectorField& Sf(mesh.Sf());

    forAll(fCells,fI)
    {
        label cfI = fCells[fI];

        if (cfI != cfgI)
        {
            if (cfI < mesh.nInternalFaces())
            {
                if (mesh.owner()[cfI] == cellI)
                {
                    gradSn += ps[cfI] * (Sf[cfI] & nf);
                }
                else
                {
                    gradSn -= ps[cfI] * (Sf[cfI] & nf);
                }
            }
        }
    }

    // other boundary faces contribution
    forAll(otherCellBoundaryFaces, ocfI)
    {
        label patchI = otherCellBoundaryFaces[ocfI].first();
        label faceII = otherCellBoundaryFaces[ocfI].second();

        gradSn += ps.boundaryField()[patchI][faceII]
            *(Sf.boundaryField()[patchI][faceII] & nf);
    }

    // add self contribution
    gradSn += p[cellI]*(cPatch.patch().magSf()[faceI]);


    vector deltaF = cPatch.patch().Cf()[faceI] - mesh.C()[cellI];
    scalar deltaFn = (nf & deltaF);

    //limit sub-volume for stability on tets/prisms
    scalar subV = cPatch.patch().magSf()[faceI] * deltaFn;
    subV = min(0.5*mesh.V()[cellI], subV);

    gradSn /= (mesh.V()[cellI] - subV);

    //damp gradient magnitude based on non-orthogonality
    scalar ninety = constant::mathematical::pi*0.5;

    if (startDampingAngle_ < -ninety)
    {
        gradSn = 0.0;
        alphaG_[faceI] = 0;
    }
    else
    {
        scalar dFndF = deltaFn/mag(deltaF);
        dFndF = min(1,dFndF);
        dFndF = max(-1,dFndF);
        scalar angle = Foam::acos(dFndF);

        scalar a = ninety/max(SMALL, (zeroGradientAngle_ - startDampingAngle_));
        scalar b = -a*startDampingAngle_;

        scalar northalpha(Foam::cos(min(max(0.0, (angle*a + b)), ninety)));
        gradSn *= northalpha;
        alphaG_[faceI] *= northalpha;
    }

    return gradSn;

}


vector wallPressureFvPatchScalarField::cellGradient
(
    label cellI,
    const Tuple2<label, label>& startCellBF,
    const List<Tuple2<label, label>>& wallPressureBFs,
    const List<Tuple2<label, label>>& otherCellBFs
) const
{
    // reference to mesh
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    // reference to field this belongs to
    const volScalarField& p
    (
        this->db().lookupObject<volScalarField>
        (
            this->internalField().name()
        )
    );

    const labelList& fCells = mesh.cells()[cellI];
    vector cGrad(vector::zero);

    //Add internal face contributions
    forAll(fCells,fI)
    {
        label cf = fCells[fI];

        if (cf < mesh.nInternalFaces())
        {
            if (mesh.owner()[cf] == cellI)
            {
                cGrad += wallPressureFvPatchScalarField
                    ::facePressure->operator[](cf) * mesh.Sf()[cf];
            }
            else
            {
                cGrad -= wallPressureFvPatchScalarField
                    ::facePressure->operator[](cf) * mesh.Sf()[cf];
            }
        }
    }

    //Add start face contributions
    label patchI = startCellBF.first();
    label faceI = startCellBF.second();
    cGrad
        += p.boundaryField()[patchI][faceI]
        *mesh.boundary()[patchI].Sf()[faceI];

    //Add other wallPressure boundary contributions
    forAll(wallPressureBFs, bfI)
    {
        patchI = wallPressureBFs[bfI].first();
        faceI = wallPressureBFs[bfI].second();
        cGrad
            += p.boundaryField()[patchI][faceI]
            *mesh.boundary()[patchI].Sf()[faceI];
    }

    //Add other boundary contributions
    forAll(otherCellBFs, bfI)
    {
        patchI = otherCellBFs[bfI].first();
        faceI = otherCellBFs[bfI].second();
        cGrad
            += p.boundaryField()[patchI][faceI]
            *mesh.boundary()[patchI].Sf()[faceI];
    }

    cGrad /= mesh.V()[cellI];

    return cGrad;
}


void wallPressureFvPatchScalarField::extrapolateGradient()
{

    // reference to mesh
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    // reference to field this belongs to
    const volScalarField& p
    (
        this->db().lookupObject<volScalarField>
        (
            this->internalField().name()
        )
    );

    //increment local index
    locInterpIndex_++;

    // check if this field surface interpolation exists/is up to date

    if
    (
        (
            locInterpIndex_ != wallPressureFvPatchScalarField::globInterpIndex
            || !wallPressureFvPatchScalarField::facePressure.valid()
        )
        ||
        (mesh.changing() && updateAddressing_)
    )
    {
        wallPressureFvPatchScalarField::facePressure.reset
        (
            fvc::interpolate(p).ptr()
        );
        wallPressureFvPatchScalarField::globInterpIndex = locInterpIndex_;
    }


    //face normals
    tmp<vectorField> nf(this->patch().nf());

    const Foam::labelUList& faceCells(this->patch().faceCells());

    forAll(*this, pfI)
    {
        label fcI = faceCells[pfI];

        // reset gradient relaxation factor
        alphaG_[pfI] = 1;

        scalar gradf = faceNormalGradient
        (
            *this, pfI, nf()[pfI],
            multiFaceAddrPtr_->operator[](pfI),
            fcI
        );

        if (wallPressureFaceAddrPtr_->operator[](pfI).size() > 0)
        {

            forAll(wallPressureFaceAddrPtr_->operator[](pfI), mfaI)
            {
                const Tuple2<label, label>& wallPressurePF
                (
                    wallPressureFaceAddrPtr_->operator[](pfI)[mfaI]
                );

                const fvPatch& nbrPatch
                (
                    mesh.boundary()[wallPressurePF.first()]
                );

                label nbrfI = wallPressurePF.second();

                scalar normalDot
                    = nf()[pfI]
                    & (nbrPatch.Sf()[nbrfI]/nbrPatch.magSf()[nbrfI]);

                scalar ffacing(1 - mag(normalDot));
                gradf *= ffacing;
                alphaG_[pfI] *= ffacing;
            }
        }

        gradient()[pfI] = (1-relaxG_)*gradient()[pfI] + relaxG_*gradf;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallPressureFvPatchScalarField::wallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    locInterpIndex_(-1),
    multiFaceAddrPtr_(),
    wallPressureFaceAddrPtr_(),
    startDampingAngle_(0.707),
    zeroGradientAngle_(0.0),
    alphaG_(p.size(), 1),
    relaxG_(0.5),
    updateAddressing_(false)
{}


wallPressureFvPatchScalarField::wallPressureFvPatchScalarField
(
    const wallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    locInterpIndex_(ptf.locInterpIndex_),
    multiFaceAddrPtr_(),
    wallPressureFaceAddrPtr_(),
    startDampingAngle_(ptf.startDampingAngle_),
    zeroGradientAngle_(ptf.zeroGradientAngle_),
    alphaG_(p.size(), 1),
    relaxG_(ptf.relaxG_),
    updateAddressing_(ptf.updateAddressing_)
{}


wallPressureFvPatchScalarField::wallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    locInterpIndex_(0),
    multiFaceAddrPtr_(),
    wallPressureFaceAddrPtr_(),
    startDampingAngle_
    (
        degToRad(dict.lookupOrDefault<scalar>("startDampingAngle", 30))
    ),
    zeroGradientAngle_
    (
        degToRad(dict.lookupOrDefault<scalar>("zeroGradientAngle", 60))
    ),
    alphaG_(p.size(), 1),
    relaxG_(dict.lookupOrDefault<scalar>("gradientRelaxationFactor", 0.5)),
    updateAddressing_
    (
        dict.lookupOrDefault<bool>("updateAddressing", false)
    )
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        gradient() = 0.0;
    }

    fvPatchField<scalar>::operator=
    (
        patchInternalField()+gradient()/patch().deltaCoeffs()
    );
}


wallPressureFvPatchScalarField::wallPressureFvPatchScalarField
(
    const wallPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    locInterpIndex_(wbppsf.locInterpIndex_),
    multiFaceAddrPtr_(),
    wallPressureFaceAddrPtr_(),
    startDampingAngle_(wbppsf.startDampingAngle_),
    zeroGradientAngle_(wbppsf.zeroGradientAngle_),
    alphaG_(wbppsf.size(), 1),
    relaxG_(wbppsf.relaxG_),
    updateAddressing_(wbppsf.updateAddressing_)
{}


wallPressureFvPatchScalarField::wallPressureFvPatchScalarField
(
    const wallPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    locInterpIndex_(wbppsf.locInterpIndex_),
    multiFaceAddrPtr_(),
    wallPressureFaceAddrPtr_(),
    startDampingAngle_(wbppsf.startDampingAngle_),
    zeroGradientAngle_(wbppsf.zeroGradientAngle_),
    alphaG_(wbppsf.size(), 1),
    relaxG_(wbppsf.relaxG_),
    updateAddressing_(wbppsf.updateAddressing_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallPressureFvPatchScalarField::~wallPressureFvPatchScalarField()
{
    if (facePressure.valid())
    {
        facePressure.clear();
    }

    //fixedGradientFvPatchScalarField::~fixedGradientFvPatchScalarField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if
    (
        !multiFaceAddrPtr_.valid()
        || !wallPressureFaceAddrPtr_.valid()
        ||
        (
            this->patch().boundaryMesh().mesh().changing()
            && updateAddressing_
        )
    )
    {
        calculateMultiFaceAdressing();
    }

    if (this->patch().boundaryMesh().mesh().changing())
    {
        alphaG_.resize(this->patch().size(), 1);
    }

    // extrapolate surface normal gradient to boundary
    extrapolateGradient();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void wallPressureFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<scalar>::operator=
    (
        this->patchInternalField()
        //+ gradient_/this->patch().deltaCoeffs()
        + this->gradient() * (this->patch().delta() & this->patch().nf())
    );

    fvPatchScalarField::evaluate();
}


tmp<scalarField> wallPressureFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


tmp<scalarField> wallPressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}

void wallPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    os.writeEntry("startDampingAngle", radToDeg(startDampingAngle_));
    os.writeEntry("zeroGradientAngle", radToDeg(zeroGradientAngle_));
    os.writeEntry("gradientRelaxationFactor", relaxG_);
    if (updateAddressing_)
    {
        os.writeEntry("updateAddressing", updateAddressing_);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
        fvPatchScalarField,
        wallPressureFvPatchScalarField
    );
}


// ************************************************************************* //
