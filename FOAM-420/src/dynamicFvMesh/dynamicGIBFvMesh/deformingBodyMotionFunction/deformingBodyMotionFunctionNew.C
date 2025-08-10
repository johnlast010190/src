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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "deformingBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::deformingBodyMotionFunction>
Foam::deformingBodyMotionFunction::New
(
    const fvMesh& mesh,
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
{
    const word motionType(SBMFCoeffs.lookup("deformingBodyMotionFunction"));
    const dictionary& SBMFCoeffsi = SBMFCoeffs.subDict(motionType + "Coeffs");

    Info<< "Selecting solid-body motion function " << motionType << endl;

    const_cast<Time&>(mesh.time()).libs().open
    (
        "libdeformingBodyMotionFunctions.so"
    );
    const_cast<Time&>(mesh.time()).libs().open
    (
        SBMFCoeffs,
        "libs",
        dictionaryConstructorTable_()
    );

    const auto ctor = ctorTableLookup("deformingBodyMotionFunction", dictionaryConstructorTable_(), motionType);
    return autoPtr<deformingBodyMotionFunction>(ctor(mesh, SBMFCoeffsi, runTime));
}


// ************************************************************************* //
