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
    (c) 2016-2017 Wikki Ltd
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteArea/ddtSchemes/boundedBackwardFaDdtScheme/boundedBackwardFaDdtScheme.H"
#include "finiteArea/fac/facDiv.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar boundedBackwardFaDdtScheme::deltaT_() const
{
    return mesh().time().deltaT().value();
}


scalar boundedBackwardFaDdtScheme::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const dimensionedScalar dt
)
{
    // No change compared to backward

    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<areaScalarField> tdtdt
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                dimensionedScalar
                (
                    "0",
                    dt.dimensions()/dimTime,
                    0.0
                )
            )
        );

        tdtdt.ref().primitiveFieldRef() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().S0() - coefft00*mesh().S00())/mesh().S()
        );

        return tdtdt;
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                dimensionedScalar
                (
                    "0",
                    dt.dimensions()/dimTime,
                    0.0
                ),
                calculatedFaPatchScalarField::typeName
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const dimensionedScalar dt
)
{
    // No change compared to backward

    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    tmp<areaScalarField> tdtdt0
    (
        new areaScalarField
        (
            ddtIOobject,
            mesh(),
            -rDeltaT*(coefft0 - coefft00)*dt
        )
    );

    if (mesh().moving())
    {
        tdtdt0.ref().primitiveFieldRef() = (-rDeltaT.value()*dt.value())*
        (
            (coefft0*mesh().S0() - coefft00*mesh().S00())/mesh().S()
        );
    }

    return tdtdt0;
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos0(phict) - pos0(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*vf.primitiveField() -
                    (
                        coefft0.primitiveField()
                        *vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()
                        *vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                    coefft.boundaryField()*vf.boundaryField() -
                    (
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos0(phict) - pos0(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                  - (
                        coefft0.primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - (
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const dimensionedScalar& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos0(phict) - pos0(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.primitiveField() -
                    (
                        coefft0.primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft.boundaryField()*vf.boundaryField() -
                    (
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}

tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const dimensionedScalar& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            vf.oldTime().oldTime()
          - vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                vf.oldTime()
              - vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos0(phict) - pos0(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                   -(
                        coefft0.primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        vf.oldTime().oldTime().primitiveField()
                       *mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*rho.value()*
                (
                   -(
                        coefft0.boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt
(
    const areaScalarField& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            rho.oldTime().oldTime()*vf.oldTime().oldTime()
          - rho.oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", rho.dimensions()*vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos0(phict) - pos0(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*rho.primitiveField()*vf.primitiveField() -
                    (
                        coefft0.primitiveField()*
                        rho.oldTime().primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        rho.oldTime().oldTime().primitiveField()
                       *vf.oldTime().oldTime().primitiveField()*mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                    coefft.boundaryField()*vf.boundaryField() -
                    (
                        coefft0.boundaryField()*
                        rho.oldTime().boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        rho.oldTime().oldTime().boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<areaScalarField> boundedBackwardFaDdtScheme::facDdt0
(
    const areaScalarField& rho,
    const areaScalarField& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt0("+rho.name()+','+vf.name()+')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    areaScalarField phict
    (
        mag
        (
            rho.oldTime().oldTime()*vf.oldTime().oldTime()
          - rho.oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime()
        )/
        (
            mag
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            )
          + dimensionedScalar("small", rho.dimensions()*vf.dimensions(), SMALL)
        )
    );

    areaScalarField limiter(pos0(phict) - pos0(phict - scalar(1)));

    areaScalarField coefft(scalar(1) + limiter*deltaT/(deltaT + deltaT0));
    areaScalarField coefft00(limiter*sqr(deltaT)/(deltaT0*(deltaT + deltaT0)));
    areaScalarField coefft0(coefft + coefft00);

    if (mesh().moving())
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                  - (
                        coefft0.primitiveField()*
                        rho.oldTime().primitiveField()*
                        vf.oldTime().primitiveField()*mesh().S0()
                      - coefft00.primitiveField()*
                        rho.oldTime().oldTime().primitiveField()*
                        vf.oldTime().oldTime().primitiveField()*mesh().S00()
                    )/mesh().S()
                ),
                rDeltaT.value()*
                (
                  - (
                        coefft0.boundaryField()*
                        rho.oldTime().boundaryField()*
                        vf.oldTime().boundaryField()
                      - coefft00.boundaryField()*
                        rho.oldTime().oldTime().boundaryField()*
                        vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<areaScalarField>
        (
            new areaScalarField
            (
                ddtIOobject,
                rDeltaT*
                (
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


tmp<faScalarMatrix> boundedBackwardFaDdtScheme::famDdt
(
    const areaScalarField& vf
)
{
    tmp<faScalarMatrix> tfam
    (
        new faScalarMatrix
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    faScalarMatrix& fam = tfam.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    scalarField phict
    (
        mag
        (
            vf.oldTime().oldTime().internalField()
          - vf.oldTime().oldTime().oldTime().internalField()
        )/
        (
            mag
            (
                vf.oldTime().internalField()
              - vf.oldTime().oldTime().internalField()
            )
            + SMALL
        )
    );

    scalarField limiter(pos0(phict) - pos0(phict - 1.0));

    scalarField coefft(1.0 + limiter*deltaT/(deltaT + deltaT0));
    scalarField coefft00(limiter*deltaT*deltaT/(deltaT0*(deltaT + deltaT0)));
    scalarField coefft0(coefft + coefft00);

    fam.diag() = (coefft*rDeltaT)*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*
        (
            coefft0*vf.oldTime().primitiveField()*mesh().S0()
          - coefft00*vf.oldTime().oldTime().primitiveField()
           *mesh().S00()
        );
    }
    else
    {
        fam.source() = rDeltaT*mesh().S()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


tmp<faScalarMatrix> boundedBackwardFaDdtScheme::famDdt
(
    const dimensionedScalar& rho,
    const areaScalarField& vf
)
{
    tmp<faScalarMatrix> tfam
    (
        new faScalarMatrix
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faScalarMatrix& fam = tfam.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    scalarField phict
    (
        mag
        (
            vf.oldTime().oldTime().internalField()
          - vf.oldTime().oldTime().oldTime().internalField()
        )/
        (
            mag
            (
                vf.oldTime().internalField()
              - vf.oldTime().oldTime().internalField()
            )
            + SMALL
        )
    );

    scalarField limiter(pos0(phict) - pos0(phict - 1.0));

    scalarField coefft(1.0 + limiter*deltaT/(deltaT + deltaT0));
    scalarField coefft00(limiter*deltaT*deltaT/(deltaT0*(deltaT + deltaT0)));
    scalarField coefft0(coefft + coefft00);

    fam.diag() = (coefft*rDeltaT*rho.value())*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*rho.value()*
        (
            coefft0*vf.oldTime().primitiveField()*mesh().S0()
          - coefft00*vf.oldTime().oldTime().primitiveField()
           *mesh().S00()
        );
    }
    else
    {
        fam.source() = rDeltaT*mesh().S()*rho.value()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}


tmp<faScalarMatrix> boundedBackwardFaDdtScheme::famDdt
(
    const areaScalarField& rho,
    const areaScalarField& vf
)
{
    tmp<faScalarMatrix> tfam
    (
        new faScalarMatrix
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faScalarMatrix& fam = tfam.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Calculate unboundedness indicator
    // Note: all times moved by one because access to internal field
    // copies current field into the old-time level.
    scalarField phict
    (
        mag
        (
            rho.oldTime().oldTime().internalField()*
            vf.oldTime().oldTime().internalField()
          - rho.oldTime().oldTime().oldTime().internalField()*
            vf.oldTime().oldTime().oldTime().internalField()
        )/
        (
            mag
            (
                rho.oldTime().internalField()*
                vf.oldTime().internalField()
              - rho.oldTime().oldTime().internalField()*
                vf.oldTime().oldTime().internalField()
            )
            + SMALL
        )
    );

    scalarField limiter(pos0(phict) - pos0(phict - 1.0));

    scalarField coefft(1.0 + limiter*deltaT/(deltaT + deltaT0));
    scalarField coefft00(limiter*deltaT*deltaT/(deltaT0*(deltaT + deltaT0)));
    scalarField coefft0(coefft + coefft00);

    fam.diag() = (coefft*rDeltaT)*rho.primitiveField()*mesh().S();

    if (mesh().moving())
    {
        fam.source() = rDeltaT*
        (
            coefft0*rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().S0()
          - coefft00*rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()*mesh().S00()
        );
    }
    else
    {
        fam.source() = rDeltaT*mesh().S()*
        (
            coefft0*rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()
          - coefft00*rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfam;
}

tmp<edgeScalarField> boundedBackwardFaDdtScheme::meshPhi
(
    const areaScalarField& vf
)
{
    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    scalar coefftn_0 = 1 + coefft0_00;

    return tmp<edgeScalarField>
    (
        new edgeScalarField
        (
            IOobject
            (
                mesh().phi().name(),
                mesh()().time().timeName(),
                mesh()(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            coefftn_0*mesh().phi() - coefft0_00*mesh().phi().oldTime()
        )
    );
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundedBackwardFaDdtScheme, 0);

faDdtScheme<scalar>::addIstreamConstructorToTable<boundedBackwardFaDdtScheme>
    addboundedBackwardFaDdtSchemeIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
