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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "physicalProperties/physicalProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicalProperties::physicalProperties
(
    const dictionary& physicalPropertiesDict
)
:
    rho_(nullptr),
    Cp_(nullptr),
    lambda_(nullptr),
    Prt_(nullptr),
    Pr_(nullptr),
    physicalProperties_(physicalPropertiesDict)
{
    read(physicalPropertiesDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::physicalProperties::rho() const
{
    if (rho_.valid())
    {
        return rho_();
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Illegal call to undefined property, rho. "
            << "Access of this property requires that the quantaty be "
            << "defined in the constant/transportProperties dictionary."
            << exit(FatalError);
    }

    return rho_();
}

const Foam::dimensionedScalar& Foam::physicalProperties::Cp() const
{
    if (Cp_.valid())
    {
        return Cp_();
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Illegal call to undefined property, Cp. "
            << "Access of this property requires that the quantaty be "
            << "defined in the constant/transportProperties dictionary."
            << exit(FatalError);
    }

    return Cp_();
}

const Foam::dimensionedScalar& Foam::physicalProperties::lambda() const
{
    if (lambda_.valid())
    {
        return lambda_();
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Illegal call to undefined property, lambda. "
            << "Access of this property requires that the quantaty be "
            << "defined in the constant/transportProperties dictionary."
            << exit(FatalError);
    }

    return lambda_();
}

const Foam::dimensionedScalar& Foam::physicalProperties::Prt() const
{
    if (Prt_.valid())
    {
        return Prt_();
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Illegal call to undefined property, Prt. "
            << "Access of this property requires that the quantaty be "
            << "defined in the constant/transportProperties dictionary."
            << exit(FatalError);
    }

    return Prt_();
}

const Foam::dimensionedScalar& Foam::physicalProperties::Pr() const
{
    if (Pr_.valid())
    {
        return Pr_();
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Illegal call to undefined property, Prt. "
            << "Access of this property requires that the quantaty be "
            << "defined in the constant/transportProperties dictionary."
            << exit(FatalError);
    }

    return Pr_();
}


bool Foam::physicalProperties::read(const dictionary& physicalPropertiesDict)
{
    physicalProperties_ = physicalPropertiesDict;

    rho_.reset
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rho",
                physicalProperties_,
                dimDensity,
                1.0
            )
        )
    );

    if (physicalProperties_.found("Cp"))
    {
        Cp_.reset
        (
            new dimensionedScalar
            (
                "Cp",
                dimEnergy/dimMass/dimTemperature,
                physicalProperties_.lookup("Cp")
            )
        );
    }

    // usage of "k" as conductivity variable is deprecated
    if (physicalProperties_.found("k"))
    {
        lambda_.reset
        (
            new dimensionedScalar
            (
                "k",
                dimEnergy/dimTime/dimLength/dimTemperature,
                physicalProperties_.lookup("k")
            )
        );
    }
    //usage of lambda as conduction variable will be supersceded by kappa
    else if (physicalProperties_.found("lambda"))
    {
        lambda_.reset
        (
            new dimensionedScalar
            (
                "lambda",
                dimEnergy/dimTime/dimLength/dimTemperature,
                physicalProperties_.lookup("lambda")
            )
        );
    }


    Prt_.reset
    (
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "Prt",
                physicalProperties_,
                dimless,
                0.85
            )
        )
    );


    if (physicalProperties_.found("Pr"))
    {
        Pr_.reset
        (
            new dimensionedScalar
            (
                "Pr",
                dimless,
                physicalProperties_.lookup("Pr")
            )
        );
    }

    return true;
}


// ************************************************************************* //
