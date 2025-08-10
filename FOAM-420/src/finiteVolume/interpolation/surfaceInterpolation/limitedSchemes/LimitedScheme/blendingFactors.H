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
    (c) 2019 Esi Ltd.

Namespace
    Foam::blendingFactors

Description
    Class to create blending factors for interfaceCapturing schemes.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches).

SourceFiles
    blendingFactors.C

\*---------------------------------------------------------------------------*/

#ifndef blendingFactors_H
#define blendingFactors_H

#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace blendingFactors
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline tmp<surfaceScalarField> calculateTheta
(
    const volScalarField&
);


class null
{
public:

    null()
    {}

    inline tmp<surfaceScalarField> operator()
    (
        const volScalarField& vsf
    ) const
    {
        const fvMesh& mesh = vsf.mesh();
        tmp<surfaceScalarField> tblendingFactor
        (
            new surfaceScalarField
            (
                IOobject
                (
                    vsf.name() + "BlendingFactor",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimless
            )
        );

        surfaceScalarField& blending = tblendingFactor.ref();
        blending *= 0.;

        return tblendingFactor;
    }
};


class cosine
{
public:

    cosine()
    {}

    inline tmp<surfaceScalarField> operator()
    (
        const volScalarField&
    ) const;
};


class root2
{
public:

    root2()
    {}

    inline tmp<surfaceScalarField> operator()
    (
        const volScalarField&
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingFactors
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/blendingFactors.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //




