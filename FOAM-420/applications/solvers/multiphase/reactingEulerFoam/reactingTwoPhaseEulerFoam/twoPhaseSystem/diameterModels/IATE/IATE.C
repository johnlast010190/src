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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "IATE.H"
#include "diameterModels/IATE/IATEsources/IATEsource/IATEsource.H"
#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcAverage.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "global/constants/fundamental/fundamentalConstants.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(IATE, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        IATE,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATE::IATE
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    kappai_
    (
        IOobject
        (
            IOobject::groupName("kappai", phase.name()),
            phase_.time().timeName(),
            phase_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.mesh()
    ),
    dMax_("dMax", dimLength, diameterProperties_),
    dMin_("dMin", dimLength, diameterProperties_),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        diameterProperties_
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.time().timeName(),
            phase_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dsm()
    ),
    sources_
    (
        diameterProperties_.lookup("sources"),
        IATEsource::iNew(*this)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::IATE::~IATE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::IATE::dsm() const
{
    return max(6/max(kappai_, 6/dMax_), dMin_);
}

void Foam::diameterModels::IATE::correct()
{
    volScalarField alphaAv
    (
        max
        (
            fvc::average(phase_ + phase_.oldTime()),
            residualAlpha_
        )
    );

    // Initialise the accumulated source term to the dilatation effect
    fvScalarMatrix R
    (
       -fvm::SuSp
        (
            ((1.0/3.0)/alphaAv)
           *(
                fvc::ddt(phase_) + fvc::div(phase_.alphaPhi())
              - phase_.continuityError()/phase_.rho()
            ),
            kappai_
        )
    );

    // Accumulate the run-time selectable sources
    forAll(sources_, j)
    {
        R += sources_[j].R(alphaAv, kappai_);
    }

    fv::options& fvOptions(fv::options::New(phase_.mesh()));

    // Construct the interfacial curvature equation
    fvScalarMatrix kappaiEqn
    (
        fvm::ddt(kappai_) + fvm::div(phase_.phi(), kappai_)
      - fvm::Sp(fvc::div(phase_.phi()), kappai_)
     ==
        R
      + fvOptions(kappai_)
    );

    kappaiEqn.relax();

    fvOptions.constrain(kappaiEqn);

    kappaiEqn.solve();

    // Update the Sauter-mean diameter
    d_ = dsm();
}


bool Foam::diameterModels::IATE::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    diameterProperties_.lookup("dMax") >> dMax_;
    diameterProperties_.lookup("dMin") >> dMin_;

    // Re-create all the sources updating number, type and coefficients
    PtrList<IATEsource>
    (
        diameterProperties_.lookup("sources"),
        IATEsource::iNew(*this)
    ).transfer(sources_);

    return true;
}


// ************************************************************************* //
