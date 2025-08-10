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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/turbulentTemperatureCoupledBaffleMixed/turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    boundaryKappa(db(), patch(), "undefined", "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),
    acht_(false),
    Cie_(0.02),
    switchTol_(0.5)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    boundaryKappa(db(), patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),
    acht_(ptf.acht_),
    Cie_(ptf.Cie_),
    switchTol_(ptf.switchTol_)
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    boundaryKappa(db(), patch(), dict),
    TnbrName_(dict.lookup("Tnbr")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),
    acht_(dict.lookupOrDefault<Switch>("acht", false)),
    Cie_(dict.lookupOrDefault<scalar>("Cie", 0.02)),
    switchTol_(dict.lookupOrDefault<scalar>("switchTol", 0.0))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("thicknessLayers"))
    {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    boundaryKappa(db(), patch(), wtcsf),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    kappaLayers_(wtcsf.kappaLayers_),
    contactRes_(wtcsf.contactRes_),
    acht_(wtcsf.acht_),
    Cie_(wtcsf.Cie_),
    switchTol_(wtcsf.switchTol_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    tmp<scalarField> intFld = patchInternalField();
    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& nbrField =
    refCast
    <
        const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            TnbrName_
        )
    );

    // Swap to obtain full local values of neighbour internal field
    tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
    tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

    if (contactRes_ == 0.0)
    {
        nbrIntFld.ref() = nbrField.patchInternalField();
        nbrKDelta.ref() = nbrField.kappa()*nbrPatch.deltaCoeffs();
    }
    else
    {
        nbrIntFld.ref() = nbrField;
        nbrKDelta.ref() = contactRes_;
    }

    mpp.distribute(nbrIntFld.ref());
    mpp.distribute(nbrKDelta.ref());

    tmp<scalarField> myKDelta = kappa()*patch().deltaCoeffs();


    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    this->refValue() = nbrIntFld();
    this->refGrad() = 0.0;
    this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

    if (acht_)
    {
        /*
            Description of acht idea:
            Combine harmonic and gradient boundary conditions
            using factor Cie

            Tboundary
                = Cie*(T_harmonic) + (1-Cie)*(T_gradient)

            So Cie = 1 is fully harmonic, while Cie = 0 is pure gradient

            Since we can only use the mixedFvPatchField framework
            we need to recalculate refValue, valueFraction and refGradient
            such that the correct boundary temperature is predicted from:

            Tboundary = valueFraction*refValue
                + (1-valueFraction)*(Tinternal + grad(T)sn*delta)

            We do this by assuming:
            valueFraction_new = Cie*valueFraction_old + Cie*(1-Cie)

            and then calculating refValue and refGradient based on the
            current internal value
        */


        //update wall current temperature
        scalarField Tp
        (
            valueFraction()*refValue()
            + (1.0 - valueFraction())*intFld()
        );

        scalarField deltaTib( Tp - intFld() );

        scalarField beta( Cie_*valueFraction() + Cie_*(1-Cie_) );
        scalar Ciemo = 1 - Cie_;
        scalar Ciemo2 = Ciemo*Ciemo;
        scalar CieCiemo = Cie_ * Ciemo;


        forAll(Tp, pfI)
        {

            if (intFld()[pfI] > nbrIntFld()[pfI] - switchTol_)
            {
                refValue()[pfI]
                    = (nbrIntFld()[pfI]*Cie_*valueFraction()[pfI]
                    + Tp[pfI] * CieCiemo)/beta[pfI];

                refGrad()[pfI]
                    = patch().deltaCoeffs()[pfI] *
                    (
                        (Cie_*(1-valueFraction()[pfI])*intFld()[pfI]
                        + Ciemo2*(intFld()[pfI] + deltaTib[pfI]))
                        / (1 - beta[pfI]) - intFld()[pfI]
                    );

                valueFraction()[pfI] = beta[pfI];
            }
        }
    }


    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa()*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    thicknessLayers_.writeEntry("thicknessLayers", os);
    kappaLayers_.writeEntry("kappaLayers", os);

    if (acht_)
    {
        os.writeEntry("acht", acht_);
        os.writeEntry("Cie", Cie_);
        os.writeEntry("switchTol", switchTol_);
    }

    boundaryKappa::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
);

// New name
addNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixed
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
