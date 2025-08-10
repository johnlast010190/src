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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistribute.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    qs_(p.size()),
    solidDict_(),
    solidPtr_(nullptr),
    QPrevious_(p.size()),
    QRelaxation_(1)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase
    (
        p.patch(),
        ptf,
        isA<directFvPatchFieldMapper>(mapper)
      ? dynamic_cast<const directFvPatchFieldMapper&>(mapper).addressing()
      : dynamic_cast<const generalFvPatchFieldMapper&>(mapper)
       .directAddressing()
    ),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(mapper(ptf.thickness_)),
    qs_(mapper(ptf.qs_)),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    QPrevious_(mapper(ptf.QPrevious_)),
    QRelaxation_(ptf.QRelaxation_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), dict),
    mixedFvPatchScalarField(p, iF, dict, false),
    TName_("T"),
    baffleActivated_(dict.lookupOrDefault<bool>("baffleActivated", true)),
    thickness_(),
    qs_(p.size(), 0),
    solidDict_(dict),
    solidPtr_(),
    QPrevious_(p.size(), 0.0),
    QRelaxation_(dict.lookupOrDefault<scalar>("relaxation", 1))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dict, p.size());
    }

    if (dict.found("qs"))
    {
        qs_ = scalarField("qs", dict, p.size());
    }

    if (dict.found("QPrevious"))
    {
        QPrevious_ = scalarField("QPrevious", dict, p.size());
    }

    if (dict.found("refValue") && baffleActivated_)
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }

}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    QPrevious_(ptf.QPrevious_),
    QRelaxation_(ptf.QRelaxation_)
{}


template<class solidType>
thermalBaffle1DFvPatchScalarField<solidType>::
thermalBaffle1DFvPatchScalarField
(
    const thermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    qs_(ptf.qs_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_),
    QPrevious_(ptf.QPrevious_),
    QRelaxation_(ptf.QRelaxation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class solidType>
bool thermalBaffle1DFvPatchScalarField<solidType>::owner() const
{
    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    return (patchi < nbrPatchi);
}


template<class solidType>
const solidType& thermalBaffle1DFvPatchScalarField<solidType>::solid() const
{
    if (this->owner())
    {
        if (solidPtr_.empty())
        {
            solidPtr_.reset(new solidType(patch().boundaryMesh().mesh(), solidDict_));
        }
        return solidPtr_();
    }
    else
    {
        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        return nbrField.solid();
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::
baffleThickness() const
{
    if (this->owner())
    {
        if (thickness_.size() != patch().size())
        {
            FatalErrorInFunction
            << " Field thickness has not been specified "
            << " for patch " << this->patch().name()
            << " in dictionary " <<  solidDict_
            << abort(FatalError);
        }

        return thickness_;
    }
    else
    {
        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];
        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        // Using GREAT as a default value can be problematic due to roundoff
        // in the weights; rather interpolate 1/thickness
        tmp<scalarField> trThickness =
            scalar(1)/stabilise(nbrField.baffleThickness(), SMALL);
        scalarField& rThickness = trThickness.ref();
        scalarField defaultValues(patch().size(), SMALL);
        this->distribute
        (
            rThickness,
            static_cast<const UList<scalar>&>(defaultValues)
        );
        return 1/trThickness;
    }
}


template<class solidType>
tmp<scalarField> thermalBaffle1DFvPatchScalarField<solidType>::qs() const
{
    if (this->owner())
    {
         return qs_;
    }
    else
    {
        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        tmp<scalarField> tqs(new scalarField(nbrField.qs()));
        scalarField& qs = tqs.ref();

        const scalarField defaultValues(patch().size(), 0.0);
        this->distribute
        (
            qs,
            static_cast<const UList<scalar>&>(defaultValues)
        );

        return tqs;
    }
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mappedPatchBase::clearOut();

    mixedFvPatchScalarField::autoMap(m);

    if (this->owner())
    {
        m(thickness_, thickness_);
        m(qs_, qs_);
    }

    m(QPrevious_, QPrevious_);
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const thermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const thermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.rmap(tiptf.thickness_, addr);
        qs_.rmap(tiptf.qs_, addr);
    }
    QPrevious_.rmap(tiptf.QPrevious_, addr);
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(thickness_, scalar(0));
    mapper.map(qs_, scalar(0));
    mapper.map(QPrevious_, scalar(0));
}


template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    if (baffleActivated_)
    {
        const fvPatch& nbrPatch = patch().boundaryMesh()[nbrPatchi];

        const compressible::turbulenceModel& turbModel =
            db().template lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        // local properties
        const scalarField kappaw(turbModel.kappaEff(patchi));

        const fvPatchScalarField& Tp =
            patch().template lookupPatchField<volScalarField, scalar>(TName_);

        tmp<scalarField> Ti = patchInternalField();

        scalarField myh(patch().deltaCoeffs()*kappaw);

        // nbr properties
        const scalarField nbrKappaw(turbModel.kappaEff(nbrPatchi));


        const thermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const thermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        scalarField nbrTi(nbrField.patchInternalField());
        scalarField nbrTi_orig(nbrTi); //for NR solver
        this->distribute
        (
            nbrTi,
            static_cast<const UList<scalar>&>(Ti)
        );

        scalarField nbrTp =
           nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_);
        scalarField nbrTp_orig(nbrTp); //for NR solver
        this->distribute
        (
            nbrTp,
            static_cast<const UList<scalar>&>(Tp)
        );


        scalarField nbrh(nbrPatch.deltaCoeffs()*nbrKappaw);
        scalarField nbrh_orig(nbrh); //for NR solver
        this->distribute
        (
            nbrh,
            static_cast<const UList<scalar>&>(patch().deltaCoeffs()*kappaw)
        );

        // solid properties
        scalarField kappas(patch().size(), 0.0);
        forAll(kappas, i)
        {
            kappas[i] = solid().kappa(0.0, (Tp[i] + nbrTp[i])/2.0);
        }
        const scalarField KDeltaw(kappas/baffleThickness());

        scalarField Qfixed(0.5*qs_);

        if (false) //symmetric solution
        {
            // skip this for now
            // can be implemented if necessary
        }
        else //asymmetric solution, uses Tp from other side, not Ti
        {
            // Update option sources, using Newton's method to account for
            // possible face value dependence
            // We are solving
            // KDeltaw*(Twall-nbrTp) + myh*(Twall+Ti) = Qt + boundarySource(Twall)

            scalarField Tpi( this->patchInternalField() );

            scalarField C1( Qfixed + myh*Ti + KDeltaw*nbrTp );
            scalarField C2( myh + KDeltaw );

            label i = 0;
            scalar Tw_err = GREAT;
            scalar maxErr = 1e-5;
            label maxLoops = 100;

            scalarField Tw_new = *this;
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
                Tw_err = gMax(mag(Tw_new-Tw_old)/stabilise(Tw_old, SMALL));
                i++;
            }
            while (i < maxLoops && Tw_err > maxErr);

            if (i == maxLoops)
            {
                WarningInFunction
                    << "Non-convergence in Newton's method for "
                    << "temperature-dependent boundary source, patch: "
                    << this->patch().name() << nl << endl;
            }

            scalarField Qtot(Qfixed+Qopt);
            Qtot = QRelaxation_*Qtot + (1-QRelaxation_)*QPrevious_;
            QPrevious_ = Qtot;

            valueFraction() = 1.0/(1.0 + myh/KDeltaw);

            refValue() =
                1.0/C2 * (KDeltaw*nbrTp + Qtot)
                /valueFraction();


            if (debug)
            {
                const scalarField qDot(kappaw*Tp.snGrad());
                scalar Q = gSum(patch().magSf()*qDot);
                Info<< patch().boundaryMesh().mesh().name() << ':'
                    << patch().name() << ':'
                    << this->internalField().name() << " <- "
                    << nbrPatch.name() << ':'
                    << this->internalField().name() << " :"
                    << " heat[W]:" << Q
                    << " walltemperature "
                    << " min:" << gMin(*this)
                    << " max:" << gMax(*this)
                    << " avg:" << gAverage(*this)
                    << endl;
            }
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}

template<class solidType>
void thermalBaffle1DFvPatchScalarField<solidType>::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    mappedPatchBase::write(os);

    if (this->owner())
    {
        baffleThickness()().writeEntry("thickness", os);
        qs()().writeEntry("qs", os);
        solid().write(os);
    }

    QPrevious_.writeEntry("QPrevious", os);
    os.writeEntry("relaxation", QRelaxation_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
