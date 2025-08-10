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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "regionCoupledEnergyFluxFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "solidThermo/solidThermo.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/volFields/volFields.H"


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//

void regionCoupledEnergyFluxFvPatchScalarField::setThermo() const
{
    if (!thermoPtr_)
    {
        thermoPtr_ =
        (
            &this->db().lookupObject<basicThermo>
            (
                IOobject::groupName
                (
                    basicThermo::dictName,
                    internalField().group()
                )
            )
        );
    }
}


void regionCoupledEnergyFluxFvPatchScalarField::setNeighbourFieldName() const
{
    setThermo();
    const basicThermo& thermo = *thermoPtr_;
    const fvPatchField& Tpf = thermo.T().boundaryField()[patch().index()];
    const regionCoupledEnergyFluxFvPatchScalarField& Tepf =
        refCast<const regionCoupledEnergyFluxFvPatchScalarField>(Tpf);
    if (Tepf.neighbourFieldName_ != word::null)
    {
        // Neighbour field name manually specified in input - just copy it over
        neighbourFieldName_ = Tepf.neighbourFieldName_;
    }
    else
    {
        // Set coupled field name based on thermo types
        if (&Tepf == this)
        {
            // If this is the temperature patch, we have to bootstrap to find
            // the neigbouring temp field to prevent infinite recursion below
            neighbourFieldName_ = Tepf.internalField().name();
        }
        const regionCoupledEnergyFluxFvPatchScalarField& nbrTepf =
            refCast<const regionCoupledEnergyFluxFvPatchScalarField>
            (
                Tepf.neighbourFvPatchField()
            );
        nbrTepf.setThermo();
        if (nbrTepf.thermoPtr_->solvesTFromhe())
        {
            neighbourFieldName_ = nbrTepf.thermoPtr_->he().name();
        }
        else
        {
            neighbourFieldName_ = nbrTepf.thermoPtr_->T().name();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionCoupledEnergyFluxFvPatchScalarField::
regionCoupledEnergyFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFluxFvPatchScalarField(p, iF),
    thermoPtr_(nullptr)
{
    if (!regionCoupled())
    {
        return;
    }

    // The he field is constructed with this constructor, not read from
    // dictionary. Look up the dictionary settings from the corresponding
    // temperature field.

    if (iF.member() != "T")
    {
        setThermo();
        const basicThermo& thermo = *thermoPtr_;
        const fvPatchField& Tpf = thermo.T().boundaryField()[p.index()];
        const regionCoupledEnergyFluxFvPatchScalarField& Tepf =
            refCast<const regionCoupledEnergyFluxFvPatchScalarField>(Tpf);
        faceBCoeff_ = Tepf.faceBCoeff_;
        faceICoeff_ = Tepf.faceICoeff_;
        faceCorr_ = Tepf.faceCorr_;
        thicknessLayers_ = Tepf.thicknessLayers_;
        conductivityLayers_ = Tepf.conductivityLayers_;
        contactInsulance_ = Tepf.contactInsulance_;
        Qsum_ = Tepf.Qsum_;
        QsumNbr_ = Tepf.QsumNbr_;
        qURF_ = Tepf.qURF_;
        neighbourPhaseName_ = Tepf.neighbourPhaseName_;

        iF.mesh().schemes().setFluxRequired(iF.name());
    }

    // Delay working out the neighbourFieldName until it is needed,
    // as thermo is not fully constructed yet, and we need it available.
    // Also need other side patch field to be constructed.
    neighbourFieldName_ = word::null;
}


regionCoupledEnergyFluxFvPatchScalarField::
regionCoupledEnergyFluxFvPatchScalarField
(
    const regionCoupledEnergyFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCoupledFluxFvPatchScalarField(ptf, p, iF, mapper),
    thermoPtr_(nullptr)
{
    neighbourFieldName_ = word::null;
}


regionCoupledEnergyFluxFvPatchScalarField::
regionCoupledEnergyFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCoupledFluxFvPatchScalarField(p, iF, dict),
    thermoPtr_(nullptr)
{
    neighbourFieldName_ =
        dict.lookupOrDefault("neighbourFieldName", word::null);
}


regionCoupledEnergyFluxFvPatchScalarField::
regionCoupledEnergyFluxFvPatchScalarField
(
    const regionCoupledEnergyFluxFvPatchScalarField& ptf
)
:
    regionCoupledFluxFvPatchScalarField(ptf),
    thermoPtr_(nullptr)
{
    neighbourFieldName_ = word::null;
}


regionCoupledEnergyFluxFvPatchScalarField::
regionCoupledEnergyFluxFvPatchScalarField
(
    const regionCoupledEnergyFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFluxFvPatchScalarField(ptf, iF),
    thermoPtr_(nullptr)
{
    neighbourFieldName_ = word::null;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regionCoupledEnergyFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    else if (!regionCoupled())
    {
        regionCoupledFluxFvPatchScalarField::updateCoeffs();
        return;
    }

    setThermo();

    const label patchi = this->patch().index();
    const scalarField& pp = thermoPtr_->p().boundaryField()[patchi];
    const scalarField& pT = thermoPtr_->T().boundaryField()[patchi];
    if (ishe())
    {
        pMultiplier_ = thermoPtr_->Cpv(pp, pT, patchi);
        pOffset_ = thermoPtr_->he(pp, pT, patchi)/pMultiplier_-pT;
    }
    else
    {
        pMultiplier_ = scalarField(this->patch().size(), 1);
        pOffset_ = scalarField(this->patch().size(), 0);
    }
    regionCoupledFluxFvPatchScalarField::updateCoeffs();
}


tmp<scalarField> regionCoupledEnergyFluxFvPatchScalarField::boundarySources
(
    const scalarField& pf, scalarField& df
) const
{
    // Forward to the companion temperature boundary
    if (ishe())
    {
        const label patchi = this->patch().index();
        const fvPatchScalarField& pT0 =
            thermoPtr_->T().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            thermoPtr_->p().boundaryField()[patchi];
        scalarField pT( thermoPtr_->THE(pf, pp, pT0, patchi) );

        tmp<scalarField> f = pT0.boundarySources(pT, df);
        // Return derivative with respect to patch field value,
        // i.e. he not T
        scalarField dTdhe( 1.0/thermoPtr_->Cpv(pp, pT, patchi) );
        df *= dTdhe;
        return f;
    }
    else
    {
        return regionCoupledFluxFvPatchScalarField::boundarySources(pf, df);
    }
}


void regionCoupledEnergyFluxFvPatchScalarField::write(Ostream& os) const
{
    // Recover fields calculated in he
    scalarField QsumPrev(Qsum_), QsumNbrPrev(QsumNbr_);
    scalarField faceBCoeffPrev(faceBCoeff_), faceICoeffPrev(faceICoeff_);
    scalarField faceCorrPrev(faceCorr_);
    if
    (
        internalField().member() == "T"
     && this->db().foundObject<basicThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                internalField().group()
            )
        )
    )
    {
        setThermo();
        if (thermoPtr_->solvesTFromhe())
        {
            const fvPatchField& hep =
                thermoPtr_->he().boundaryField()[patch().index()];
            const regionCoupledEnergyFluxFvPatchScalarField& cwep =
                refCast<const regionCoupledEnergyFluxFvPatchScalarField>(hep);

            const_cast<scalarField&>(faceBCoeff_) = cwep.faceBCoeff_;
            const_cast<scalarField&>(faceICoeff_) = cwep.faceICoeff_;
            const_cast<scalarField&>(faceCorr_) = cwep.faceCorr_;
            const_cast<scalarField&>(Qsum_) = cwep.Qsum_;
            const_cast<scalarField&>(QsumNbr_) = cwep.QsumNbr_;
        }
    }

    regionCoupledFluxFvPatchScalarField::write(os);

    const_cast<scalarField&>(faceBCoeff_) = faceBCoeffPrev;
    const_cast<scalarField&>(faceICoeff_) = faceICoeffPrev;
    const_cast<scalarField&>(faceCorr_) = faceCorrPrev;
    const_cast<scalarField&>(Qsum_) = QsumPrev;
    const_cast<scalarField&>(QsumNbr_) = QsumNbrPrev;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchScalarField,
    regionCoupledEnergyFluxFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
