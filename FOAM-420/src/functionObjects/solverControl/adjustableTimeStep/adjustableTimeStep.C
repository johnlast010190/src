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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2012, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "adjustableTimeStep/adjustableTimeStep.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(adjustableTimeStep, 0);
    addToRunTimeSelectionTable(functionObject, adjustableTimeStep, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar Foam::functionObjects::adjustableTimeStep::getCoNum()
{
    scalar maxCoNum = 0;
    scalar meanCoNum = 0;


    //get phi
    const surfaceScalarField& phi
        = lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const volScalarField& rho = lookupObject<volScalarField>("rho");

        surfaceScalarField SfUfbyDelta
        (
            mesh_.surfaceInterpolation::deltaCoeffs()*mag(phi)
            /fvc::interpolate(rho)
        );

        maxCoNum = max(SfUfbyDelta/mesh_.magSf())
                    .value()*obr_.time().deltaT().value();

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh_.magSf()))
            .value()*obr_.time().deltaT().value();

    }
    else
    {
        surfaceScalarField SfUfbyDelta
        (
            mesh_.surfaceInterpolation::deltaCoeffs()*mag(phi)
        );

        maxCoNum = max(SfUfbyDelta/mesh_.magSf())
            .value()*obr_.time().deltaT().value();

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh_.magSf()))
            .value()*obr_.time().deltaT().value();

    }

    scalar CoNum = 0;
    if (CoType_ == "max")
    {
        CoNum = maxCoNum;
    }
    else if (CoType_ == "mean")
    {
        CoNum = meanCoNum;
    }
    else
    {
        FatalErrorInFunction
            << CoType_ << " from function: adjustableTimeStep: "
            << name()
            << ", is not a valid specification for "
            << "Courant number type." << nl
            << "Valied entries are: 'max' and 'mean'."
            << exit(FatalError);
    }


    return CoNum;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::adjustableTimeStep::adjustableTimeStep
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    targetCo_(1),
    maxDeltaT_(GREAT),
    minDeltaT_(SMALL),
    CoType_("max"),
    phiName_("phi"),
    deltaT_(0)
{
    read(dict);

    //set initial delta T

    scalar CoNum = getCoNum();

    if (CoNum > SMALL)
    {
        const_cast<Time&>(obr_.time()).setDeltaT
        (
            deltaT_ = min
            (
                targetCo_*obr_.time().deltaT().value()/CoNum,
                maxDeltaT_
            )
        );

        //exit Foam if deltaT is below user secified minDeltaT
        //if (obr_.time().deltaT().value() < minDeltaT_)
        if (deltaT_ < minDeltaT_)
        {
            FatalErrorInFunction
                << "deltaT = " <<  deltaT_ << nl
                << "minDeltaT = " << minDeltaT_ <<nl
                << "The computed initial deltaT is lower than minDeltaT"
                << exit(FatalError);
        }

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::adjustableTimeStep::~adjustableTimeStep()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::adjustableTimeStep::calculate()
{
    scalar CoNum = getCoNum();

    scalar maxDeltaTFact = targetCo_/(CoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    const_cast<Time&>(obr_.time()).setDeltaT
    (
        deltaT_ = min
        (
            deltaTFact*obr_.time().deltaT().value(),
            maxDeltaT_
        )
    );

    //exit Foam if deltaT is below user secified minDeltaT
    //if (obr_.time().deltaT().value() < minDeltaT_)
    if (deltaT_ < minDeltaT_)
    {
        FatalErrorInFunction
            << "From function: adjustableTimeStep: "
            << "The computed deltaT is lower than "
            << "minDeltaT = " << minDeltaT_
            << exit(FatalError);
    }

}


bool Foam::functionObjects::adjustableTimeStep::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    //calculate
    calculate();

    //write to screen
    Log << "deltaT = " <<  obr_.time().deltaT().value() << endl;
    Log << endl;

    return true;
}


bool Foam::functionObjects::adjustableTimeStep::write()
{
    return true;
}


bool Foam::functionObjects::adjustableTimeStep::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    //all entries are optional for this function
    if (dict.found("targetCFL"))
    {
        targetCo_ = readScalar(dict.lookup("targetCFL"));
    }

    if (dict.found("maxDeltaT"))
    {
        maxDeltaT_ = readScalar(dict.lookup("maxDeltaT"));
    }

    if (dict.found("minDeltaT"))
    {
        minDeltaT_ = readScalar(dict.lookup("minDeltaT"));
    }
    else
    {
        minDeltaT_ = scalar(SMALL);
        WarningInFunction
            << "minDeltaT set to default (" << minDeltaT_ << " s)." << endl;
    }

    if (dict.found("maxMean"))
    {
        CoType_ = word(dict.lookup("maxMean"));
    }

    if (dict.found("phiName"))
    {
        phiName_ = word(dict.lookup("phiName"));
    }

    Log << endl;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
