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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/turbulentConvectiveTemperature/turbulentConvectiveTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "solidThermo/solidThermo.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void turbulentConvectiveTemperatureFvPatchScalarField
::calculateTotalConductivity() const
{
    at_.reset(new scalarField(this->size(), 0));
    scalarField& ta(at_());

    if (alphaWall_.valid())
    {
        ta += 1/max(VSMALL, alphaWall_());
    }

    forAll(t_, lI)
    {
        ta += t_[lI]/max(VSMALL, lambda_[lI]);
    }

    ta = 1/max(ta, VSMALL);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentConvectiveTemperatureFvPatchScalarField
::turbulentConvectiveTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    boundaryKappa
    (
        db(), patch(), "lookup", "undefined", "undefined", iF.group()
    ),
    Tinf_(p.size(), 0.0),
    alphaWall_(),
    relax_(1.0),
    layerNames_(0),
    t_(0),
    lambda_(0),
    at_(),
    qadd_()
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


turbulentConvectiveTemperatureFvPatchScalarField
::turbulentConvectiveTemperatureFvPatchScalarField
(
    const turbulentConvectiveTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    boundaryKappa(db(), patch(), ptf),
    Tinf_(mapper(ptf.Tinf_)),
    alphaWall_
    (
        ptf.alphaWall_.valid()
      ? mapper(ptf.alphaWall_()).ptr()
      : nullptr
    ),
    relax_(ptf.relax_),
    layerNames_(ptf.layerNames_),
    t_(ptf.t_.size()),
    lambda_(ptf.t_.size()),
    at_(),
    qadd_
    (
        ptf.qadd_.valid()
      ? mapper(ptf.qadd_()).ptr()
      : nullptr
    )
{
    forAll(ptf.t_, li)
    {
        t_.set(li, mapper(ptf.t_[li]).ptr());
        lambda_.set(li, mapper(ptf.lambda_[li]).ptr());
    }
}


turbulentConvectiveTemperatureFvPatchScalarField
::turbulentConvectiveTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    boundaryKappa(db(), patch(), dict, iF.group()),
    Tinf_("Tinf", dict, p.size()),
    alphaWall_
    (
        dict.found("alphaWall")
      ? new scalarField("alphaWall", dict, p.size())
      : nullptr
    ),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    layerNames_(0),
    t_(0),
    lambda_(0),
    at_(),
    qadd_
    (
        dict.found("qadd")
      ? new scalarField("qadd", dict, p.size())
      : nullptr
    )
{
    if (dict.found("layers"))
    {
        const dictionary& layerDicts(dict.subDict("layers"));
        layerNames_.setSize(layerDicts.size());
        t_.setSize(layerDicts.size());
        lambda_.setSize(layerDicts.size());

        label nLayers = 0;

        forAllConstIter(dictionary, layerDicts, iter)
        {
            // safety:
            if (!iter().isDict())
            {
                continue;
            }
            const word& key = iter().keyword();
            const dictionary& ldict = iter().dict();

            layerNames_[nLayers] = key;
            t_.set(nLayers, new scalarField("thickness", ldict, p.size()));
            lambda_.set(nLayers, new scalarField("lambda", ldict, p.size()));

            nLayers++;
        }
    }

    refValue() = Tinf_;
    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
    }

    refGrad() = 0.0;

    valueFraction() = 0.0;
    if (dict.found("valueFraction"))
    {
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


turbulentConvectiveTemperatureFvPatchScalarField
::turbulentConvectiveTemperatureFvPatchScalarField
(
    const turbulentConvectiveTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    boundaryKappa(db(), patch(), tppsf),
    Tinf_(tppsf.Tinf_),
    alphaWall_
    (
        tppsf.alphaWall_.valid()
      ? new scalarField(tppsf.alphaWall_())
      : nullptr
    ),
    relax_(tppsf.relax_),
    layerNames_(tppsf.layerNames_),
    t_(tppsf.t_),
    lambda_(tppsf.lambda_),
    at_(),
    qadd_
    (
        tppsf.qadd_.valid()
      ? new scalarField(tppsf.qadd_())
      : nullptr
    )
{}


turbulentConvectiveTemperatureFvPatchScalarField
::turbulentConvectiveTemperatureFvPatchScalarField
(
    const turbulentConvectiveTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    boundaryKappa(db(), patch(), tppsf),
    Tinf_(tppsf.Tinf_),
    alphaWall_
    (
        tppsf.alphaWall_.valid()
      ? new scalarField(tppsf.alphaWall_())
      : nullptr
    ),
    relax_(tppsf.relax_),
    layerNames_(tppsf.layerNames_),
    t_(tppsf.t_),
    lambda_(tppsf.lambda_),
    at_(),
    qadd_
    (
        tppsf.qadd_.valid()
      ? new scalarField(tppsf.qadd_())
      : nullptr
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& turbulentConvectiveTemperatureFvPatchScalarField
::alphaWall()
{
    at_.clear();
    if (!alphaWall_.valid())
    {
        alphaWall_.reset(new scalarField(this->size(), 0.0));
    }

    return alphaWall_();
}

List<word>& turbulentConvectiveTemperatureFvPatchScalarField
::layerNames()
{
    at_.clear();
    return layerNames_;
}

PtrList<scalarField>&
turbulentConvectiveTemperatureFvPatchScalarField::thickness()
{
    at_.clear();
    return t_;
}

PtrList<scalarField>&
turbulentConvectiveTemperatureFvPatchScalarField::lambda()
{
    at_.clear();
    return lambda_;
}

const scalarField&
turbulentConvectiveTemperatureFvPatchScalarField::conductivity() const
{
    if (!at_.valid())
    {
        calculateTotalConductivity();
    }

    return at_();
}

void turbulentConvectiveTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(Tinf_, Tinf_);

    if (alphaWall_.valid())
    {
        m(alphaWall_(), alphaWall_());
    }
    if (qadd_.valid())
    {
        m(qadd_(), qadd_());
    }

    forAll(t_, li)
    {
        m(t_[li], t_[li]);
        m(lambda_[li], lambda_[li]);
    }
}


void turbulentConvectiveTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const turbulentConvectiveTemperatureFvPatchScalarField& tiptf =
        refCast<const turbulentConvectiveTemperatureFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);

    if (tiptf.alphaWall_.valid())
    {
        alphaWall().rmap(tiptf.alphaWall_(), addr);
    }

    if (tiptf.qadd_.valid())
    {
        qadd_().rmap(tiptf.qadd_(), addr);
    }

    forAll(tiptf.t_, li)
    {
        t_[li].rmap(tiptf.t_[li], addr);
        lambda_[li].rmap(tiptf.lambda_[li], addr);
    }
}


void turbulentConvectiveTemperatureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(Tinf_, scalar(0));
    if (alphaWall_.valid())
    {
        mapper.map(alphaWall_(), scalar(0));
    }
    if (qadd_.valid())
    {
        mapper.map(qadd_(), scalar(0));
    }
    forAll(t_, li)
    {
        mapper.map(t_[li], scalar(0));
        mapper.map(lambda_[li], scalar(0));
    }
}


void turbulentConvectiveTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    vectorField nfKappaw;

    // non-wall temperature dependent heat flux
    scalarField Qfixed(this->size(), 0.0);

    nfKappaw = this->nfKappa();

    if (this->KMethod() != mtLookup)
    {
        const basicThermo& thermo = this->thermo();
        if
        (
            isA<solidThermo>(thermo)
         && !refCast<const solidThermo>(thermo).isotropic()
        )
        {
            word gradTName = "grad("+thermo.T().name()+")";
            autoPtr<volVectorField> pGrad;
            if (!db().foundObject<volVectorField>(gradTName))
            {
                serialThreads::pauseSwitching();
                pGrad.set
                (
                    fvc::grad(thermo.T()).ptr()
                );
                pGrad().rename(gradTName);
                serialThreads::resumeSwitching();
            }

            const volVectorField& gradT =
                db().lookupObject<volVectorField>(gradTName);
            vectorField nfKappaTang(this->nfKappa());
            nfKappaTang -= (nfKappaTang & patch().nf())*patch().nf();
            Qfixed +=
                (nfKappaTang & gradT.boundaryField()[patchi]);
        }
    }

    const scalarField& Tp = *this;
    const scalarField Tpi = this->patchInternalField()();

    const scalarField& alphaCond = conductivity();
    scalarField alphaConv
    (
        (nfKappaw & patch().nf()) * patch().deltaCoeffs()
    );

    if (qadd_.valid())
    {
        Qfixed += qadd_();
    }

    // Update option sources, using Newton's method to account for
    // possible patch value dependence

    scalarField C1( Qfixed + alphaConv*Tpi + alphaCond*Tinf_ );
    scalarField C2( alphaConv + alphaCond );

    label i = 0;
    scalar maxErr = 1e-5;
    label maxLoops = 100;

    scalarField Tw_new = Tp;
    scalarField Tw_old;

    scalarField Qopt;
    do
    {
        Tw_old = Tw_new;
        scalarField bSourceDeriv;
        Qopt = boundarySources(Tw_old, bSourceDeriv);

        scalarField f( C1 + Qopt - C2*Tw_old );
        scalarField df( bSourceDeriv - C2 );

        Tw_new = Tw_old - f/df;
        i++;
    }
    while (i < maxLoops && gMax(mag(Tw_new - Tw_old)) > maxErr);

    // Initialise on the first iteration
    if (!QoptOld_.valid())
    {
        QoptOld_.reset(new scalarField(Qopt));
    }

    if (i == maxLoops)
    {
        WarningInFunction
            << "Non-convergence in Newton's method for "
            << "temperature-dependent boundary source, patch: "
            << this->patch().name() << nl
            << "Using old timestep for Qopt." << endl;
        Qopt = QoptOld_();
    }

    if (relax_ < 1.0)
    {
        Qopt = QoptOld_() + (1.0 - relax_)*(Qopt - QoptOld_());
    }

    QoptOld_.reset(new scalarField(Qopt));

    valueFraction() =
        1.0/(1.0 + alphaConv/alphaCond);

    refValue() =
        1.0/C2 * (alphaCond*Tinf_ + Qfixed + Qopt)/valueFraction();

    mixedFvPatchScalarField::updateCoeffs();
}


void turbulentConvectiveTemperatureFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    boundaryKappa::write(os);

    Tinf_.writeEntry("Tinf", os);

    if (alphaWall_.valid())
    {
        alphaWall_().writeEntry("alphaWall", os);
    }

    if (qadd_.valid())
    {
        qadd_().writeEntry("qadd", os);
    }

    os.beginBlock("layers");

    forAll(t_, li)
    {
        os.beginBlock(layerNames_[li]);

        t_[li].writeEntry("thickness", os, false);
        lambda_[li].writeEntry("lambda", os, false);

        os.endBlock();
    }

    os.endBlock();
    refValue().writeEntry("refValue", os);
    valueFraction().writeEntry("valueFraction", os);
    writeEntry("value", os);

    if (relax_ < 1.0)
    {
        os.writeEntry("relax", relax_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentConvectiveTemperatureFvPatchScalarField
);

// New name
addNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    turbulentConvectiveTemperatureFvPatchScalarField,
    turbulentConvectiveTemperature
);


// ************************************************************************* //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //