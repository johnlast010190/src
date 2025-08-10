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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "environmentalTemperatureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void environmentalTemperatureFvPatchScalarField::calculateEnvRadiation()
{
    if (rSourcesPtr_.size())
    {
        Qenv_.reset(new scalarField(size(), 0.0));

        forAllIter(PtrList<entry>, rSourcesPtr_, iter)
        {
            // safety:
            if (!iter().isDict())
            {
                continue;
            }
            const dictionary& envRadDict = iter().dict();

            scalar radiationFlux(readScalar(envRadDict.lookup("intensity")));
            vector radiationDirection(envRadDict.lookup("direction"));
            radiationDirection /= mag(radiationDirection);

            scalarField projection( -(radiationDirection & this->patch().nf()) );
            projection *= pos0(projection);
            Qenv_() += radiationFlux * projection;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

environmentalTemperatureFvPatchScalarField
::environmentalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Trad_(),
    eps_(1.0),
    sigma_(5.6704e-8),
    Tconv_(),
    alphaOutside_(),
    alphaName_("lambda"),
    CpName_("Cp"),
    rSourcesPtr_(),
    Qenv_(),
    QPtr_()
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


environmentalTemperatureFvPatchScalarField
::environmentalTemperatureFvPatchScalarField
(
    const environmentalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Trad_(),
    eps_(ptf.eps_),
    sigma_(ptf.sigma_),
    Tconv_(),
    alphaOutside_(),
    alphaName_(ptf.alphaName_),
    CpName_(ptf.CpName_),
    rSourcesPtr_(),
    Qenv_(),
    QPtr_()
{
    if (ptf.QPtr_.valid())
    {
        QPtr_.reset(mapper(ptf.QPtr_()).ptr());
    }
    else if (!ptf.alphaOutside_.valid())
    {
        Tconv_.reset(new scalar(ptf.Tconv_));
    }
    else //full monty
    {
        // following 2 must be valid
        Tconv_.reset(new scalar(ptf.Tconv_));
        alphaOutside_.reset(new scalar(ptf.alphaOutside_));

        // optionals
        if (Trad_.valid()) Trad_.reset(new scalar(ptf.Trad_));

        if (ptf.Qenv_.valid())
        {
           Qenv_.reset(mapper(ptf.Qenv_()).ptr());
           PtrList<entry> tmpRS(ptf.rSourcesPtr_);
           rSourcesPtr_.transfer(tmpRS);
        }
    }
}


environmentalTemperatureFvPatchScalarField
::environmentalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    Trad_(),
    eps_(dict.lookupOrDefault<scalar>("emissivity", 1)),
    sigma_(dict.lookupOrDefault<scalar>("sigma", 5.6704e-8)),
    Tconv_(),
    alphaOutside_(),
    alphaName_(dict.lookupOrDefault<word>("alphaName", "lambda")),
    CpName_(dict.lookupOrDefault<word>("CpName", "Cp")),
    rSourcesPtr_(),
    Qenv_(),
    QPtr_()
{
    // logic
    // 1. If Q is specified ignore everything else and use fixed heat flux
    // 2. If Tconv is specified, but not alphaConv, use Tconv fixed temp
    // 3. If Tconv and alphaConv is specified use everything available

    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("q"))
    {
        QPtr_.reset(new scalarField(dict.lookup("q")));
        refValue() = 300;
    }
    else if (dict.found("Twall"))
    {
        Tconv_.reset(new scalar(readScalar(dict.lookup("Twall"))));
        refValue() = Tconv_;
        valueFraction() = 1.0;
    }
    else
    {
        Tconv_.reset(new scalar(readScalar(dict.lookup("Tenv"))));

        alphaOutside_.reset
        (
            new scalar(readScalar(dict.lookup("alphaConv")))
        );

        if (dict.found("Trad"))
        {
            Trad_.reset(new scalar(readScalar(dict.lookup("Trad"))));
        }

        if (dict.found("directedRadiation"))
        {
            PtrList<entry> tsource
            (
                (dict.lookupEntryPtr("directedRadiation",false,false))->stream()
            );

            rSourcesPtr_.transfer(tsource);

            calculateEnvRadiation();
        }
        refValue() = Tconv_;
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


environmentalTemperatureFvPatchScalarField
::environmentalTemperatureFvPatchScalarField
(
    const environmentalTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    Trad_(),
    eps_(ptf.eps_),
    sigma_(ptf.sigma_),
    Tconv_(),
    alphaOutside_(),
    alphaName_(ptf.alphaName_),
    CpName_(ptf.CpName_),
    rSourcesPtr_(),
    Qenv_(),
    QPtr_()
{
    if (ptf.QPtr_.valid())
    {
        QPtr_.reset(new scalarField(ptf.QPtr_()));
    }
    else if (!ptf.alphaOutside_.valid())
    {
        Tconv_.reset(new scalar(ptf.Tconv_));
    }
    else //full monty
    {
        // following 2 must be valid
        Tconv_.reset(new scalar(ptf.Tconv_));
        alphaOutside_.reset(new scalar(ptf.alphaOutside_));

        // optionals
        if (Trad_.valid()) Trad_.reset(new scalar(ptf.Trad_));

        if (ptf.Qenv_.valid())
        {
           Qenv_.reset(new scalarField(ptf.Qenv_()));
           PtrList<entry> tmpRS(ptf.rSourcesPtr_);
           rSourcesPtr_.transfer(tmpRS);
        }
    }
}


environmentalTemperatureFvPatchScalarField
::environmentalTemperatureFvPatchScalarField
(
    const environmentalTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    Trad_(),
    eps_(ptf.eps_),
    sigma_(ptf.sigma_),
    Tconv_(),
    alphaOutside_(),
    alphaName_(ptf.alphaName_),
    CpName_(ptf.CpName_),
    rSourcesPtr_(),
    Qenv_(),
    QPtr_()
{
    if (ptf.QPtr_.valid())
    {
        QPtr_.reset(new scalarField(ptf.QPtr_()));
    }
    else if (!ptf.alphaOutside_.valid())
    {
        Tconv_.reset(new scalar(ptf.Tconv_));
    }
    else //full monty
    {
        // following 2 must be valid
        Tconv_.reset(new scalar(ptf.Tconv_));
        alphaOutside_.reset(new scalar(ptf.alphaOutside_));

        // optionals
        if (Trad_.valid()) Trad_.reset(new scalar(ptf.Trad_));

        if (ptf.Qenv_.valid())
        {
           Qenv_.reset(new scalarField(ptf.Qenv_()));
           PtrList<entry> tmpRS(ptf.rSourcesPtr_);
           rSourcesPtr_.transfer(tmpRS);
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void environmentalTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    if (QPtr_.valid())
    {
        m(QPtr_(), QPtr_());
    }
    else if (alphaOutside_.valid())
    {
        if (Qenv_.valid())
        {
           m(Qenv_(), Qenv_());
        }
    }
}


void environmentalTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const environmentalTemperatureFvPatchScalarField& tiptf =
        refCast<const environmentalTemperatureFvPatchScalarField>(ptf);

    if (tiptf.QPtr_.valid())
    {
        QPtr_->rmap(tiptf.QPtr_, addr);
    }
    else if (tiptf.alphaOutside_.valid())
    {
        if (tiptf.Qenv_.valid())
        {
            Qenv_->rmap(tiptf.Qenv_, addr);
        }
    }
}


void environmentalTemperatureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMapGIB(mapper);
    if (QPtr_.valid())
    {
        mapper.map(QPtr_(), scalar(0));
    }
    if (alphaOutside_.valid() && Qenv_.valid())
    {
        mapper.map(Qenv_(), scalar(0));
    }
}


void environmentalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);
    const fvPatchField<scalar>& Cpp =
        patch().lookupPatchField<volScalarField, scalar>(CpName_);

    if (QPtr_.valid())
    {
        //fixed heat flux mode
        valueFraction() = 0.0;

        refGrad() = QPtr_()/(alphap * Cpp);
    }
    else if (!alphaOutside_.valid())
    {
        //fixed temperature mode
        valueFraction() = 1.0;

        refValue() = Tconv_();
    }
    else //full monty
    {
        const scalarField& Tp( *this );
        scalarField Tpi( this->patchInternalField() );

          //use Newton's method to solve for Twall

        scalarField Tw_old = Tp;
        scalarField Tw_new = Tp;

        scalarField alphaInside( alphap * patch().deltaCoeffs() );

        //get the right dimensions for alphaInside = W/m2/K
        if
        (
            alphap.internalField().dimensions()
            != dimensionSet(1,1,-3,1,0,0,0)
        )
        {
            FatalErrorInFunction
                << "This boundary does not yet support diffusion coefficient"
                << " with dimensions: "
                << alphap.internalField().dimensions()
                << nl << "Currently supported dimensions: "
                << dimensionSet(1,1,-3,1,0,0,0)
                << exit(FatalError);
        }

        scalar alphaRad = eps_*sigma_;

        // non-wall temperature dependent heat flux
        scalarField C1( alphaInside*Tpi + alphaOutside_ * Tconv_() );
        if (Qenv_.valid())
        {
            C1 += Qenv_();
        }
        if (Trad_.valid())
        {
            C1 += alphaRad*Trad_*Trad_*Trad_*Trad_;
        }

        scalarField C2( alphaInside + alphaOutside_ );

        label i = 0;
        scalar Tp_err = GREAT;
        scalar maxerr = 1e-5;
        label maxloops = 100;

        do
        {
            Tw_old = Tw_new;

            scalarField f( C1 - C2*Tw_old );

            scalarField df( -C2 );

            if (Trad_.valid())
            {
                f -= alphaRad*pow4(Tw_old);
                df -= 4*alphaRad*pow3(Tw_old);
            }

            Tw_new = Tw_old - f/df;

            Tp_err = max(mag(Tw_new-Tw_old)/Tw_old);
            i++;

        }
        while (i < maxloops && Tp_err > maxerr);

        //set the mixed boundary coeffs
        if (Trad_.valid())
        {
            C2 += alphaRad*pow3(Tw_new);
        }
        C1 -= alphaInside*Tpi;

        valueFraction() = 1 - alphaInside/C2;

        refValue() = C1/C2/valueFraction();

    }

    mixedFvPatchScalarField::updateCoeffs();
}


void environmentalTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    if (QPtr_.valid())
    {
        QPtr_->writeEntry("q", os);
    }
    else if (!alphaOutside_.valid())
    {
        os.writeEntry("Twall", Tconv_());
    }
    else //full monty
    {
        // following 2 must be valid
        os.writeEntry("Tenv", Tconv_());
        os.writeEntry("alphaConv", alphaOutside_());

        // optionals
        if (Trad_.valid())
        {
            os.writeEntry("Trad", Trad_());
        }

        if (Qenv_.valid())
        {
            os.beginBlock("directedRadiation");

            forAll(rSourcesPtr_, i)
            {
                // safety:
                if (!rSourcesPtr_(i)->isDict())
                {
                    continue;
                }
                const dictionaryEntry* rsEntry =
                    dynamic_cast<const dictionaryEntry*>(rSourcesPtr_(i));

                const dictionary& envRadDict = rsEntry->dict();
                const word& name = rsEntry->keyword();

                os.beginBlock(name);
                os.writeEntry
                (
                    "intensity",
                    envRadDict.lookup<scalar>("intensity")
                );
                os.writeEntry
                (
                    "direction",
                    envRadDict.lookup<vector>("direction")
                );
                os.endBlock();
            }
            os.endBlock();
        }
    }

    writeEntryIfDifferent<scalar>(os, "emissivity", 1.0, eps_);
    writeEntryIfDifferent<scalar>(os, "sigma", 5.6704e-8, sigma_);
    writeEntryIfDifferent<word>(os, "alphaName", "lambda", alphaName_);
    writeEntryIfDifferent<word>(os, "CpName", "Cp", CpName_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    environmentalTemperatureFvPatchScalarField
);

// ************************************************************************* //

} // End namespace incompressible
} // End namespace Foam


// ************************************************************************* //
