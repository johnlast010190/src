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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSpecieMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSpecieMixture::basicSpecieMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const objectRegistry& obr,
    const word& phaseName
)
:
    speciesMassFractions(thermoDict, specieNames, obr, phaseName)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::basicSpecieMixture::W() const
{
    const PtrList<volScalarField>& Y(speciesMassFractions::Y());

    tmp<volScalarField> trW
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("W", Y[0].group()),
                Y[0].time().timeName(),
                Y[0].mesh()
            ),
            Y[0].mesh(),
            dimensionedScalar("zero", dimless, 0)
        )
    );

    volScalarField& rW = trW.ref();

    forAll(Y, i)
    {
        rW += Y[i]/W(i);
    }

    rW = 1.0/rW;

    return trW;
}


// ************************************************************************* //
