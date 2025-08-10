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

\*---------------------------------------------------------------------------*/

#include "SCOPEXiEq.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(SCOPEXiEq, 0);
    addToRunTimeSelectionTable(XiEqModel, SCOPEXiEq, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::SCOPEXiEq::SCOPEXiEq
(
    const dictionary& XiEqProperties,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),
    XiEqCoef_(readScalar(XiEqModelCoeffs_.lookup("XiEqCoef"))),
    XiEqExp_(readScalar(XiEqModelCoeffs_.lookup("XiEqExp"))),
    lCoef_(readScalar(XiEqModelCoeffs_.lookup("lCoef"))),
    SuMin_(0.01*Su.average()),
    uPrimeCoef_(readScalar(XiEqModelCoeffs_.lookup("uPrimeCoef"))),
    subGridSchelkin_
    (
        readBool(XiEqModelCoeffs_.lookup("subGridSchelkin"))
    ),
    MaModel
    (
        Su.mesh().lookupObject<IOdictionary>("combustionProperties"),
        thermo
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::SCOPEXiEq::~SCOPEXiEq()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::SCOPEXiEq::XiEq() const
{
    tmp<volScalarField> tk(turbulence_.k());
    const volScalarField& k = tk();
    tmp<volScalarField> tepsilon(turbulence_.epsilon());
    const volScalarField& epsilon = tepsilon();

    volScalarField up(sqrt((2.0/3.0)*k));
    if (subGridSchelkin_)
    {
        up.primitiveFieldRef() += calculateSchelkinEffect(uPrimeCoef_);
    }

    volScalarField l(lCoef_*sqrt(3.0/2.0)*up*k/epsilon);
    volScalarField Rl(up*l*thermo_.rhou()/thermo_.muu());

    volScalarField upBySu(up/(Su_ + SuMin_));
    volScalarField K(0.157*upBySu/sqrt(Rl));
    volScalarField Ma(MaModel.Ma());

    tmp<volScalarField> tXiEq
    (
        new volScalarField
        (
            IOobject
            (
                "XiEq",
                epsilon.time().timeName(),
                epsilon.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            epsilon.mesh(),
            dimensionedScalar("XiEq", dimless, 0.0)
        )
    );
    volScalarField& xieq = tXiEq.ref();

    forAll(xieq, celli)
    {
        if (Ma[celli] > 0.01)
        {
            xieq[celli] =
                XiEqCoef_*pow(K[celli]*Ma[celli], -XiEqExp_)*upBySu[celli];
        }
    }

    volScalarField::Boundary& xieqBf = xieq.boundaryFieldRef();

    forAll(xieq.boundaryField(), patchi)
    {
        scalarField& xieqp = xieqBf[patchi];
        const scalarField& Kp = K.boundaryField()[patchi];
        const scalarField& Map = Ma.boundaryField()[patchi];
        const scalarField& upBySup = upBySu.boundaryField()[patchi];

        forAll(xieqp, facei)
        {
            if (Ma[facei] > 0.01)
            {
                xieqp[facei] =
                    XiEqCoef_*pow(Kp[facei]*Map[facei], -XiEqExp_)
                   *upBySup[facei];
            }
        }
    }

    return tXiEq;
}


bool Foam::XiEqModels::SCOPEXiEq::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqModelCoeffs_.lookup("XiEqCoef") >> XiEqCoef_;
    XiEqModelCoeffs_.lookup("XiEqExp") >> XiEqExp_;
    XiEqModelCoeffs_.lookup("lCoef") >> lCoef_;
    XiEqModelCoeffs_.lookup("uPrimeCoef") >> uPrimeCoef_;
    XiEqModelCoeffs_.lookup("subGridSchelkin") >> subGridSchelkin_;

    return true;
}


// ************************************************************************* //
