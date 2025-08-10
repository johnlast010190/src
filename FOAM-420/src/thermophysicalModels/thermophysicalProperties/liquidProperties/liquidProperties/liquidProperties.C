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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "liquidProperties/liquidProperties/liquidProperties.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidProperties, 0);
    defineRunTimeSelectionTable(liquidProperties,);
    defineRunTimeSelectionTable(liquidProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidProperties::liquidProperties
(
    scalar W,
    scalar Tc,
    scalar Pc,
    scalar Vc,
    scalar Zc,
    scalar Tt,
    scalar Pt,
    scalar Tb,
    scalar dipm,
    scalar omega,
    scalar delta
)
:
    thermophysicalProperties(W),
    Tc_(Tc),
    Pc_(Pc),
    Vc_(Vc),
    Zc_(Zc),
    Tt_(Tt),
    Pt_(Pt),
    Tb_(Tb),
    dipm_(dipm),
    omega_(omega),
    delta_(delta)
{}


Foam::liquidProperties::liquidProperties(const dictionary& dict)
:
    thermophysicalProperties(dict),
    Tc_(readScalar(dict.lookup("Tc"))),
    Pc_(readScalar(dict.lookup("Pc"))),
    Vc_(readScalar(dict.lookup("Vc"))),
    Zc_(readScalar(dict.lookup("Zc"))),
    Tt_(readScalar(dict.lookup("Tt"))),
    Pt_(readScalar(dict.lookup("Pt"))),
    Tb_(readScalar(dict.lookup("Tb"))),
    dipm_(readScalar(dict.lookup("dipm"))),
    omega_(readScalar(dict.lookup("omega"))),
    delta_(readScalar(dict.lookup("delta")))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    const auto ctor = ctorTableLookup("liquidProperties type", ConstructorTable_(), name);
    return autoPtr<liquidProperties>(ctor());
}


Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    const word& liquidPropertiesTypeName = dict.dictName();

    if (dict.found("defaultCoeffs"))
    {
        // Backward-compatibility

        if (Switch(dict.lookup("defaultCoeffs")))
        {
            return New(liquidPropertiesTypeName);
        }
        else
        {
            const auto ctor = ctorTableLookup("liquidProperties type", dictionaryConstructorTable_(), liquidPropertiesTypeName);
            return autoPtr<liquidProperties>
            (
                ctor
                (
                    dict.optionalSubDict(liquidPropertiesTypeName + "Coeffs")
                )
            );
        }
    }
    else
    {
        const auto ctor = ctorTableLookup("liquidProperties type", dictionaryConstructorTable_(), liquidPropertiesTypeName);
        return autoPtr<liquidProperties>(ctor(dict));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::liquidProperties::S(scalar p, scalar T) const
{
    NotImplemented;
}


Foam::scalar Foam::liquidProperties::pvInvert(scalar p) const
{
    // Check for critical and solid phase conditions
    if (p >= Pc_)
    {
        return Tc_;
    }
    else if (p < Pt_)
    {
        if (debug)
        {
            WarningInFunction
                << "Pressure below triple point pressure: "
                << "p = " << p << " < Pt = " << Pt_ <<  nl << endl;
        }
        return -1;
    }

    // Set initial upper and lower bounds
    scalar Thi = Tc_;
    scalar Tlo = Tt_;

    // Initialise T as boiling temperature under normal conditions
    scalar T = Tb_;

    while ((Thi - Tlo) > 1.0e-4)
    {
        if ((pv(p, T) - p) <= 0)
        {
            Tlo = T;
        }
        else
        {
            Thi = T;
        }

        T = (Thi + Tlo)*0.5;
    }

    return T;
}


void Foam::liquidProperties::readIfPresent(const dictionary &dict)
{
    thermophysicalProperties::readIfPresent(dict);
    dict.readIfPresent("Tc", Tc_);
    dict.readIfPresent("Pc", Pc_);
    dict.readIfPresent("Vc", Vc_);
    dict.readIfPresent("Zc", Zc_);
    dict.readIfPresent("Tt", Tt_);
    dict.readIfPresent("Pt", Pt_);
    dict.readIfPresent("Tb", Tb_);
    dict.readIfPresent("dipm", dipm_);
    dict.readIfPresent("omega", omega_);
    dict.readIfPresent("delta", delta_);
}


void Foam::liquidProperties::writeData(Ostream& os) const
{
    thermophysicalProperties::writeData(os);
    os  << token::SPACE
        << Tc_ << token::SPACE
        << Pc_ << token::SPACE
        << Vc_ << token::SPACE
        << Zc_ << token::SPACE
        << Tt_ << token::SPACE
        << Pt_ << token::SPACE
        << Tb_ << token::SPACE
        << dipm_ << token::SPACE
        << omega_<< token::SPACE
        << delta_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const liquidProperties& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
