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
    (c) 2010-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "solidNodeHeatFluxTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidNodeHeatFluxTemperatureFvPatchScalarField::
solidNodeHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(),
    layersProp_()
{}


solidNodeHeatFluxTemperatureFvPatchScalarField::
solidNodeHeatFluxTemperatureFvPatchScalarField
(
    const solidNodeHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, false),
    layersProp_(ptf.layersProp_)

{}


solidNodeHeatFluxTemperatureFvPatchScalarField::
solidNodeHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(),
    layersProp_(dict)
{
    if (returnReduce(size(), sumOp<label>()))
    {
        if (dict.found("q"))
        {
            q_ = Function1<scalar>::New("q", dict);
        }
        else
        {
            FatalIOErrorInFunction
            (dict)
            << "Please supply 'q'"
            << exit(FatalIOError);
        }

        gradient() = 0.0;

        if (dict.found("gradient"))
        {
            gradient() = scalarField("gradient", dict, p.size());
        }

        Field<scalar>::operator=
        (
            this->patchInternalField() + gradient()/this->patch().deltaCoeffs()
        );
    }
}


solidNodeHeatFluxTemperatureFvPatchScalarField::
solidNodeHeatFluxTemperatureFvPatchScalarField
(
    const solidNodeHeatFluxTemperatureFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    q_(thftpsf.q_, false),
    layersProp_(thftpsf.layersProp_)
{}


solidNodeHeatFluxTemperatureFvPatchScalarField::
solidNodeHeatFluxTemperatureFvPatchScalarField
(
    const solidNodeHeatFluxTemperatureFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    q_(thftpsf.q_, false),
    layersProp_(thftpsf.layersProp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidNodeHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void solidNodeHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void solidNodeHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();
    scalarField Qt(patch().size());

    const scalar t = db().time().timeOutputValue();
    scalar qnode = q_->value(t);

    const incompressibleTurbulenceModel& turbModel =
        db().lookupObject<incompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const tmp<volScalarField> talphaEff
        = turbModel.alphaEff()*turbModel.rho();
    const volScalarField& alphaEff = talphaEff;
    const scalarField& alphaEffp = alphaEff.boundaryField()[patchi];

    const tmp<volScalarField> tCp = turbModel.Cp();
    const volScalarField& Cp = tCp;
    const scalarField& Cpp = Cp.boundaryField()[patchi];

    const scalarField& Tp = *this;

    scalarField alphaConv( Cpp*alphaEffp* patch().deltaCoeffs() ); // [W/m2 K]
    scalarField rCond // [m2 K/W]
    (
        alphaConv.size(),
        layersProp_.totalResistance()
    );

    scalarField alphaCond( 1.0/max(rCond, VSMALL) );  // [W/m2K]

    // fluid temperature next to the wall
    scalarField Tf( this->patchInternalField() );

    const scalarField& faceArea = patch().magSf();
    scalarField alphaCondA( alphaCond*faceArea ); // [W/K]
    scalarField alphaConvA( alphaConv*faceArea ); // [W/K]

    // solid node temperature
    scalar Tnode;

    // Add fvOption sources
    scalarField bSourceDeriv;
    scalarField Qfixed(boundarySources(Tp, bSourceDeriv));

    const radiation::fvDOM* domPtr(nullptr);

    // radiation handling
    if
    (
        db().foundObject<radiation::radiationModel>("radiationProperties")
     && db().foundObject<volScalarField>("qin")
    )
    {
        const radiation::radiationModel& radiation =
            db().lookupObject<radiation::radiationModel>("radiationProperties");

        // Remove -emittedRadiatnIntensity included in Qfixed from boundary
        // source
        if (domPtr) Qfixed += domPtr->emittedRadiantIntensity(patchi, Tp)();

        domPtr = &(refCast<const radiation::fvDOM>(radiation));
    }

    scalarField qabs( Qfixed*faceArea ); // [W]

    // init
    scalar TpAve = gAverage(Tp);

    scalar Tnode_old = TpAve;
    scalar Tnode_new = TpAve;

    scalarField oneAlpha( 1.0 + alphaConvA/max(alphaCondA, VSMALL) );

    scalarField Tw = Tp;

    scalarField pef(this->size(), 0.0);

    if (domPtr)
    {
        pef = domPtr->emittedRadiantIntensity(patchi, Tw);
    }

    label j = 0;
    scalar Tnode_err = GREAT;
    scalar maxOuterr = 1e-5;
    label maxOutloops = 100;

    do //outer loop for Tnode
    {
        Tnode_old = Tnode_new;

        // coeffs linearisation of radiative term wrt Tnode
        // considering that Twi = Tnode - qi/alphaCondAi
        scalarField f1( 4*pef*faceArea/Tw );
        scalarField f0( pef*faceArea - f1*Tnode_old );

        scalar sumNum = 0;
        scalar sumDen = 0;

        forAll(qabs, facei)
        {
            sumNum += (qabs[facei]/oneAlpha[facei])
                + (alphaConvA[facei]*Tf[facei]/oneAlpha[facei])
                - (f0[facei]/oneAlpha[facei]);

            sumDen += (alphaConvA[facei] + f1[facei])/oneAlpha[facei];
        }

        reduce(sumNum, sumOp<scalar>());
        reduce(sumDen, sumOp<scalar>());

        sumNum += qnode;

        Tnode_new = sumNum/sumDen;
        Tnode_err = mag(Tnode_new-Tnode_old)/Tnode_old;
        j++;

        //Now we know Tnode
        //so we can use Newton's method to solve for Twall (and get qi)
        Tnode = Tnode_new;

        scalarField Tw_old = Tp;
        scalarField Tw_new = Tp;

        scalarField C1( Qfixed + alphaConv*Tf + alphaCond*Tnode );
        scalarField C2( alphaConv + alphaCond );

        label i = 0;
        scalar Tp_err = GREAT;
        scalar maxerr = 1e-5;
        label maxloops = 100;

        // this update of pef slightly improves convergence
        if (domPtr)
        {
            pef = domPtr->emittedRadiantIntensity(patchi, Tw_old);
        }

        do // solve for Tw
        {
            Tw_old = Tw_new;

            scalarField f( C1 - C2*Tw_old - pef );

            scalarField df( -C2 - 4*pef/Tw_old );

            Tw_new = Tw_old - f/df;

            // update radiation using Tw_new
            if (domPtr)
            {
                pef = domPtr->emittedRadiantIntensity(patchi, Tw_new);
            }

            Tp_err = max(mag(Tw_new-Tw_old)/Tw_old);
            i++;
        }
        while (i < maxloops && Tp_err > maxerr);

        Tw = Tw_new; // update for outer loop
    }
    while (j < maxOutloops && Tnode_err > maxOuterr); // Tnode

    // Info<<"Tnode: "<<Tnode<<" Twall: "<<gAverage(Tw)<<endl;

    // we fix the convection in this bc
    Qt = alphaConv*(Tw - Tf);

    // fix gradient
    gradient() = Qt/(Cpp*alphaEffp);
    fixedGradientFvPatchScalarField::updateCoeffs();

}


void solidNodeHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    gradient().writeEntry("gradient", os);
    q_->writeData(os);
    writeEntry("value", os);

    layersProp_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    solidNodeHeatFluxTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
