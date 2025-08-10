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

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBodyMotionFunction> Foam::solidBodyMotionFunction::New
(
    const dictionary& SBMFCoeffs,
    const Time& runTime,
    const word& frameName
)
{
    const word motionType
    (
        SBMFCoeffs.lookup(solidBodyMotionFunction::typeName)
    );

    Info<< "Selecting motion function " << motionType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "interpolation type",
            dictionaryConstructorTable_(),
            motionType
        );

    return
        autoPtr<solidBodyMotionFunction>(ctor(SBMFCoeffs, runTime, frameName));
}


Foam::autoPtr<Foam::solidBodyMotionFunction> Foam::solidBodyMotionFunction::New
(
    const objectRegistry& obr,
    const dictionary& SBMFCoeffs,
    const word& frameName
)
{
    // Allow old keyword (solidBodyMotionFunction)
    word motionType;
    const word solidBodyName(solidBodyMotionFunction::typeName);
    if (SBMFCoeffs.found("type"))
    {
        motionType = SBMFCoeffs.lookup<word>("type");
    }
    else
    {
        motionType = SBMFCoeffs.lookup<word>(solidBodyName);
    }

    Info<< "Selecting motion function " << motionType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "motionFunctionType",
            registryConstructorTable_(),
            motionType
        );

    return autoPtr<solidBodyMotionFunction>(ctor(obr, SBMFCoeffs, frameName));
}

// ************************************************************************* //
