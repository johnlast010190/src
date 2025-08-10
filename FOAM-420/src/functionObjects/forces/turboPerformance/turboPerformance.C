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
    (c) D. Boger, B. Lewis M. Auvinen, H. Nilsson
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2021-2022 Esi

\*---------------------------------------------------------------------------*/

#include "turboPerformance/turboPerformance.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turboPerformance, 0);
    addToRunTimeSelectionTable(functionObject, turboPerformance, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turboPerformance::turboPerformance
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fluidPower(name, runTime, dict, true),
    omega_(vector::zero),
    movingWallSectors_(dict.lookupOrDefault<label>("movingWallSectors", 1))
{
    turboPerformance::read(dict);
    resetFile(turboPerformance::typeName);
    turboPerformance::writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::turboPerformance::~turboPerformance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::turboPerformance::writeFileHeader(Ostream& os) const
{
    const volScalarField& p
        = lookupObject<volScalarField>(pName_);
    os << "# Time" << tab;
    if (p.dimensions() == dimPressure)
    {
        os << "Flow (Kg/s)" << tab;
    }
    else
    {
        os << "Flow (m^3/s)" << tab;
    }
    os << "Density (kg/m^3)" << tab << "Gravity (m/s^2)" << tab
       << "Omega (rad^-1)" << tab << "Head (m)" <<tab
       << "Shaft Power (W)" << tab << "Efficiency (%)" << tab
       << "F_x (N)" << tab <<"F_y (N)" << tab <<"F_z (N)" << tab
       << "M_x (Nm)" << tab <<"M_y (Nm)" << tab <<"M_z (Nm)" << endl;
}


bool Foam::functionObjects::turboPerformance::read(const dictionary& dict)
{
    Log << turboPerformance::type() << " " << name() <<  " read:" << nl;

    fluidPower::read(dict);

    // For now omega (in rad/s) is the only additional info we need
    dict.lookup("omega") >> omega_;

    // Rotating patch
    patchSet_ =
        mesh_.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));

    return true;
}


bool Foam::functionObjects::turboPerformance::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    forces::calcForcesMoment();
    dEmHead dEmH  =  fluidPower::calcDEmHead();

    vector totForce = forces::forceEff();
    vector totMoment = forces::momentEff();

    if (movingWallSectors_>1)
    {
        //- Compute total (360) force and moments

        vector axis = omega_/mag(omega_);
        const tensor T(axis*axis);
        scalar rotationalAngle = degToRad(360/scalar(movingWallSectors_));

        const tensor S
        (
            0, -axis.z(), axis.y(),
            axis.z(), 0, -axis.x(),
            -axis.y(), axis.x(), 0
        );

        const tensor revTPos
        (
            T
          + cos(rotationalAngle)*(tensor::I - T)
          + sin(rotationalAngle)*S
        );

        vector forceRev(totForce);
        vector sumForce(forceRev);

        vector momRev(totMoment);
        vector sumMom(momRev);

        for (int iSec=2; iSec<=movingWallSectors_; iSec++)
        {
            forceRev = transform(revTPos, forceRev);
            sumForce += forceRev;

            momRev = transform(revTPos, momRev);
            sumMom += momRev;
        }
        totForce = sumForce;
        totMoment = sumMom;
    }

    // Shaft power (W)
    scalar TOmega = fabs( totMoment & omega_);

    // Pump Efficiency (%)
    scalar eff = ( dEmH.first() / TOmega ) * scalar(100);

    if (fluidPower::turbine_) // Turbine flag from fluidPower (Bryan)
    {
        // Turbine Efficiency (%)
        eff = ( TOmega / dEmH.first() ) * scalar(100);
    }

    scalar rhoave(rho()().weightedAverage(mesh_.V()).value());
    scalar grav(9.81);
    if (mesh_.foundObject<uniformDimensionedVectorField>("g"))
    {
        const uniformDimensionedVectorField& g =
            mesh_.lookupObject<uniformDimensionedVectorField>("g");

        grav = mag(g[0]);
    }

    // Tab separated output ... to avoid those irritating parenthesis. -- mikko
    file() << obr_.time().timeName() << tab
        << volflux_ << tab << rhoave << tab << grav << tab << mag(omega_) << tab
        << dEmH.second() <<tab << TOmega << tab << eff << tab
        << totForce[0] << tab << totForce[1] << tab <<totForce[2] << tab
        << totMoment[0] << tab << totMoment[1] << tab <<totMoment[2] << tab
        << endl;

    Log<< " performance data:" << nl;


    const volScalarField& p
        = lookupObject<volScalarField>(pName_);
    if (p.dimensions() == dimPressure)
    {
        Log<< "   Flow (Kg/s)     = " << volflux_ << nl;
    }
    else
    {
        Log<< "   Flow (m^3/s)     = " << volflux_ << nl;
    }

    Log<< "   Density (kg/m^3) = " << rhoave << nl
       << "   Gravity (m/s^2)  = " << grav << nl
       << "   Omega (rad^-1)   = " << mag(omega_) << nl
       << "   Head (m)     = " << dEmH.second()  << nl
       << "   TOmega (W)   = " << TOmega << nl
       << "   Eff (%)      = " << eff << nl
       << "   Forces (N)   = " << totForce << nl
       << "   Moments (Nm) = " << totMoment << nl
       << endl;

    // Write state/results information
    setResult("Flow", volflux_);
    setResult("Density", rhoave);
    setResult("Gravity", grav);
    setResult("Omega", mag(omega_));
    setResult("Head", dEmH.second());
    setResult("TOmega", TOmega);
    setResult("Efficiency", eff);
    setResult("Forces", totForce);
    setResult("Moments", totMoment);

    return true;
}


bool Foam::functionObjects::turboPerformance::write()
{
    return true;
}

// ************************************************************************* //
