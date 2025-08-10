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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/greyDiffusiveRadiation/greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "submodels/boundaryRadiationProperties/boundaryRadiationProperties.H"

#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "global/constants/constants.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    externalRadiationLoad_(false)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    externalRadiationLoad_(ptf.externalRadiationLoad_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    externalRadiationLoad_(dict.lookupOrDefault<bool>("externalRadiationLoad", true))
{
    if (dict.found("refValue"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    externalRadiationLoad_(ptf.externalRadiationLoad_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    externalRadiationLoad_(ptf.externalRadiationLoad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    const label patchi = patch().index();

    if (dom.nLambda() != 1)
    {
        FatalErrorInFunction
            << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;

    const vectorField n(patch().nf());

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    //we use a consistent direction for incident and emited radiation
    //otherwise flux imbalanced result in hot spots
    const scalarField nRayDave(-n & ray.dAve());

    //flux direction
    scalarField io(pos0(nRayDave));

    //reference boundary radiation container and get relevant properties
    const boundaryRadiationProperties& boundaryRadiation
        = dom.boundaryProperties();

    //emissivity
    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patchi)
    );
    const scalarField& emissivity = temissivity();

    //transmissivity - lookup
    const tmp<scalarField> ttransmissivity
    (
        boundaryRadiation.transmissivity(patchi)
    );
    const scalarField& transmissivity = ttransmissivity();

    //transmitted sources (coming from the other side of the patch)
    scalarField Itrans(size(), 0.0);

    //boundary sources - "diffuse" transmission
    Itrans += boundaryRadiation.radiantTransmissionSource(patchi);

    //coupled sources - "specular" transmission
    Itrans += coupledRadiantTransmission();

    //incident radiation that will contribute to reflection (not transmitted)
    scalarField qirradiance(this->size(), 0.0);

    //sum all incident ordinate contributions
    for (label rayi=0; rayi < dom.nRay(); rayi++)
    {
        const vector& dAve = dom.IRay(rayi).dAve();

        const scalarField& IFace =
            dom.IRay(rayi).ILambda(lambdaId).boundaryField()[patchi];

        forAll(Iw, facei)
        {

            const scalar nRayIDavef = n[facei] & dAve;

            if (nRayIDavef > 0.0)
            {
                // q into the wall
                qirradiance[facei] += IFace[facei]*nRayIDavef;
            }
        }
    }

    //solar or other external radiative heat load
    //note: the contribution to qabsorbed is added in top level radiation solver
    if
    (
        externalRadiationLoad_
     && db().foundObject<volScalarField>(radiation.externalRadHeatFieldName_)
    )
    {
        qirradiance += patch().lookupPatchField<volScalarField,scalar>
        (
            radiation.externalRadHeatFieldName_
        );
    }


    // Heat fluxes
    // net radiative fluxes of this ray
    //(does not include external sources, e.g. solar)
    scalarField& qr = ray.qr()[patchi];

    // energy absorbed from this ray
    scalarField& qabsorbed = ray.qin()[patchi];

    // irradiance at the patch
    scalarField& qirradianceI = ray.qg()[patchi];

    //update boundary conditions
    refGrad() = 0.0;
    valueFraction() = io;

    tmp<scalarField> qe
        = boundaryRadiation.emittedRadiantFlux(Tp, patchi);
    refValue()
        = io*
        (
            (qirradiance*(1-emissivity-transmissivity)       //reflection
            + qe())/pi                                       //emission
            + Itrans                                         //transmission
        );

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();

    //Update absorption heat flux
    //calculate incident heat flux in ray direction after value update
    //does not include solar or other external radiative fluxes
    qabsorbed += (1 - io)*emissivity*Iw*nRayDave;

    //Leave qr in until confirmation that it is not needed
    qr -= (1 - io)*(1-transmissivity)*Iw*nRayDave + io*(Iw-Itrans)*nRayDave;

    //Irradiance flux
    qirradianceI += (1 - io)*Iw*nRayDave;
}


Foam::tmp<Foam::scalarField> Foam::radiation::
greyDiffusiveRadiationMixedFvPatchScalarField::
coupledRadiantTransmission() const
{

    tmp<scalarField> Isource(new scalarField(patch().size(), 0.0));

    if (isA<mappedFvPatch>(this->patch()))
    {
        // Since we're inside initEvaluate/evaluate there might be processor
        // comms underway. Change the tag we use.
        int oldTag = UPstream::msgType();
        UPstream::msgType() = oldTag+1;

        const mappedFvPatch& mpp = refCast<const mappedFvPatch>
        (
            patch()
        );
        const fvPatch& nbrPatch = mpp.nbrPatch();
        const label nbrPatchI = nbrPatch.index();

        const radiationModel& radiation =
            nbrPatch.boundaryMesh().mesh().
            lookupObject<radiationModel>("radiationProperties");

        if (isA<fvDOM>(radiation))
        {
            const fvDOM& dom(refCast<const fvDOM>(radiation));

            const boundaryRadiationProperties& boundaryRadiation
                = dom.boundaryProperties();

            if
            (
                db().foundObject<volScalarField>
                (this->internalField().name())
            )
            {
                scalarField nbrFlux =
                    nbrPatch.lookupPatchField<volScalarField, scalar>
                    (this->internalField().name());
                Isource = boundaryRadiation.transmissivity(nbrPatchI)*nbrFlux;
                const scalarField defaultValues(patch().size(), 0.0);
                Isource =
                    mpp.interpolate
                    (
                        Isource,
                        defaultValues
                    );
            }
        }

    }
    else if (isA<cyclicAMIFvPatch>(this->patch()))
    {
        // Since we're inside initEvaluate/evaluate there might be processor
        // comms underway. Change the tag we use.
        int oldTag = UPstream::msgType();
        UPstream::msgType() = oldTag+1;

        const cyclicAMIFvPatch& mpp = refCast<const  cyclicAMIFvPatch>
        (
            patch()
        );
        const fvPatch& nbrPatch = mpp.nbrPatch();
        const label nbrPatchI = nbrPatch.index();

        const radiationModel& radiation =
            nbrPatch.boundaryMesh().mesh().
            lookupObject<radiationModel>("radiationProperties");

        if (isA<fvDOM>(radiation))
        {
            const fvDOM& dom(refCast<const fvDOM>(radiation));

            const boundaryRadiationProperties& boundaryRadiation
                = dom.boundaryProperties();

            if
            (
                db().foundObject<volScalarField>
                (this->internalField().name())
            )
            {
                scalarField nbrFlux =
                    nbrPatch.lookupPatchField<volScalarField, scalar>
                    (this->internalField().name());
                Isource = boundaryRadiation.transmissivity(nbrPatchI)*nbrFlux;

                const scalarField defaultValues(patch().size(), 0.0);
                Isource =
                    mpp.interpolate
                    (
                        Isource,
                        defaultValues
                    );
            }
        }
    }

    return Isource;
}


Foam::tmp<Foam::scalarField> Foam::radiation::
greyDiffusiveRadiationMixedFvPatchScalarField::emittedRadiantIntensity
(
    const scalarField& Tp
) const
{
    tmp<Foam::scalarField> efPtr(new scalarField(size(), 0.0));
    scalarField& ef = efPtr.ref();

    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    const boundaryRadiationProperties& boundaryRadiation
        = dom.boundaryProperties();

    tmp<scalarField> qe
        = boundaryRadiation.emittedRadiantFlux(Tp, patch().index());

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    if (dom.nLambda() != 1)
    {
        FatalErrorInFunction
            << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    const scalarField nRayDave( (-patch().nf() & ray.dAve()) );

    ef += pos0(nRayDave) * nRayDave * qe()/pi;

    return efPtr;
}


void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);

    if
    (
        this->patch().type() == word("cyclicAMI")
    && !this->overridesConstraint()
    )
    {
        os.writeEntry("patchType", "cyclicAMI");
    }

    os.writeEntry("externalRadiationLoad", externalRadiationLoad_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveRadiationMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
