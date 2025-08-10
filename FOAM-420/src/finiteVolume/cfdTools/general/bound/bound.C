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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/bound/bound.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::volScalarField&
Foam::bound
(
    volScalarField& vsf,
    const dimensionedScalar& lowerBound,
    bool warn
)
{
    const scalar minVsf = min(vsf).value();

    if (minVsf < lowerBound.value())
    {
        if (warn)
        {
            Info<< "bounding " << vsf.name()
                << ", min: " << minVsf
                << " max: " << max(vsf).value()
                << " average: " << gAverage(vsf.primitiveField())
                << endl;
        }

        vsf.primitiveFieldRef() = max
        (
            max
            (
                vsf.primitiveField(),
                fvc::average(max(vsf, lowerBound))().primitiveField()
              * pos0(-vsf.primitiveField())
            ),
            lowerBound.value()
        );

        forAll(vsf.boundaryField(), pI)
        {
            vsf.boundaryFieldRef()[pI].fvPatchScalarField::operator=(
                max(vsf.boundaryField()[pI], lowerBound.value())
            );
        }
    }

    return vsf;
}


Foam::volScalarField&
Foam::bound
(
    volScalarField& vsf,
    const dimensionedScalar& lowerBound,
    const dimensionedScalar& upperBound,
    bool warn
)
{
    scalar minVsf = VGREAT;
    scalar maxVsf = -VGREAT;


    forAll(vsf.primitiveField(), i)
    {

        if (vsf[i] < minVsf) minVsf = vsf[i];
        if (vsf[i] > maxVsf) maxVsf = vsf[i];

        if (vsf[i] < lowerBound.value())
        {
            vsf[i] = lowerBound.value();
        }
        else if (vsf[i] > upperBound.value())
        {
            vsf[i] = upperBound.value();
        }

    }

    forAll(vsf.boundaryField(), pI)
    {
        vsf.boundaryFieldRef()[pI].fvPatchScalarField::operator=
        (
            max(vsf.boundaryField()[pI], lowerBound.value())
        );
        vsf.boundaryFieldRef()[pI].fvPatchScalarField::operator=
        (
            min(vsf.boundaryField()[pI], upperBound.value())
        );
    }

    if (warn)
    {
        reduce(std::tie(minVsf, maxVsf), ParallelOp<minOp<scalar>, maxOp<scalar>>{});
        // TODO: Could fuse `gAverage()` here, too...

        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << maxVsf
            << " average: " << gAverage(vsf.primitiveField())
            << endl;
    }

    return vsf;
}
// ************************************************************************* //
