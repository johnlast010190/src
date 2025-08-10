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
    (c) 2016-2023 Esi Ltd.
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "basicThermo/basicThermo.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchField.H"
#include "fields/fvPatchFields/basic/blended/blendedFvPatchField.H"
#include "derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/fixedEnergyZone/fixedEnergyZoneFvPatchScalarField.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/mixedEnergyCalculatedTemperature/mixedEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/fixedJump/fixedJumpFvPatchFields.H"
#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchFields.H"
#include "derivedFvPatchFields/energyJump/energyJump/energyJumpFvPatchScalarField.H"
#include "derivedFvPatchFields/energyJump/energyJumpAMI/energyJumpAMIFvPatchScalarField.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
#include "materialModels/materialTables/materialTables.H"
#include "materialModels/baseModels/baseModels.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
    defineRunTimeSelectionTable(basicThermo, objectRegistry);
}

Foam::word Foam::basicThermo::dictName("thermophysicalProperties");
const Foam::word Foam::basicThermo::matDictName("materialProperties");


// * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * //

Foam::word
Foam::basicThermo::heBoundaryBaseType(const fvPatchScalarField& tpf)
{
    word hbt = word::null;

    if (isA<fixedJumpFvPatchScalarField>(tpf))
    {
        const fixedJumpFvPatchScalarField& pf =
            dynamic_cast<const fixedJumpFvPatchScalarField&>(tpf);

        hbt = pf.interfaceFieldType();
    }
    else if (isA<fixedJumpAMIFvPatchScalarField>(tpf))
    {
        const fixedJumpAMIFvPatchScalarField& pf =
            dynamic_cast<const fixedJumpAMIFvPatchScalarField&>(tpf);

        hbt = pf.interfaceFieldType();
    }

    return hbt;
}


Foam::word
Foam::basicThermo::heBoundaryType(const fvPatchScalarField& tpf)
{
    word hbt = tpf.type();
    if (isA<fixedValueZoneFvPatchField<scalar>>(tpf))
    {
        hbt = fixedEnergyZoneFvPatchScalarField::typeName;
    }
    else if (isA<fixedValueFvPatchScalarField>(tpf))
    {
        hbt = fixedEnergyFvPatchScalarField::typeName;
    }
    else if
    (
        isA<zeroGradientFvPatchScalarField>(tpf)
     || isA<fixedGradientFvPatchScalarField>(tpf)
    )
    {
        hbt = gradientEnergyFvPatchScalarField::typeName;
    }
    else if
    (
        isA<mixedFvPatchScalarField>(tpf)
     || isA<mixedEnergyCalculatedTemperatureFvPatchScalarField>(tpf)
    )
    {
        hbt = mixedEnergyFvPatchScalarField::typeName;
    }
    else if (isA<fixedJumpFvPatchScalarField>(tpf))
    {
        hbt = energyJumpFvPatchScalarField::typeName;
    }
    else if (isA<fixedJumpAMIFvPatchScalarField>(tpf))
    {
        hbt = energyJumpAMIFvPatchScalarField::typeName;
    }
    else if (isA<blendedFvPatchField<scalar>>(tpf))
    {
        hbt = blendedEnergyFvPatchScalarField::typeName;
    }
    //else if (isA<genericFvPatchField<scalar>>(tpf))
    else if (tpf.type() == "generic")
    // so post-processing tools will work with unlinked boundaries
    {
        WarningInFunction
            << "Converting unknown fvPatchField type for patch "
            << tpf.patch().name()
            << " to `" << fixedEnergyFvPatchScalarField::typeName
            << "`" << nl
            << tab << "This is only a valid conversion for post-processing"
            << " tasks." << endl;
        hbt = fixedEnergyFvPatchScalarField::typeName;
    }
    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryBaseTypes()
{
    const volScalarField::Boundary& tbf =
        this->T_.boundaryField();

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        hbt[patchi] = heBoundaryBaseType(tbf[patchi]);
    }

    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryTypes()
{
    const volScalarField::Boundary& tbf =
        this->T_.boundaryField();

    wordList hbt(tbf.size());

    forAll(tbf, patchi)
    {
        hbt[patchi] = heBoundaryType(tbf[patchi]);
    }

    return hbt;
}


const Foam::fvMesh& Foam::basicThermo::mesh() const
{
    if (isA<fvSolutionRegistry>(this->db()))
    {
        return dynamic_cast<const fvSolutionRegistry&>(this->db()).mesh();
    }
    return refCast<const fvMesh>(this->db());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const objectRegistry& obr,
    const word& name,
    const dimensionSet& dims,
    const scalar& defaultValue
) const
{
    if (!obr.objectRegistry::foundObject<volScalarField>(name))
    {
        IOobject header
        (
            name,
            obr.time().timeName(),
            obr,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (header.typeHeaderOk<volScalarField>(true))
	    {
            volScalarField* fPtr
            (
                new volScalarField(header, mesh())
            );

            // Transfer ownership of this object to the objectRegistry
            fPtr->store(fPtr);
        }
        else
        {
            volScalarField* fPtr
            (
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        obr.time().timeName(),
                        obr,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar(name, dims, defaultValue),
                    zeroGradientFvPatchField<scalar>::typeName
                )
            );

            // Transfer ownership of this object to the objectRegistry
            fPtr->store(fPtr);
        }
    }

    return const_cast<volScalarField&>
    (
        obr.lookupObject<volScalarField>(name)
    );
}

Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const objectRegistry& mesh,
    const char* name,
    const dimensionSet& dims,
    const scalar& defaultValue
) const
{
    return lookupOrConstruct(mesh, word(name), dims, defaultValue);
}


Foam::volScalarField& Foam::basicThermo::lookupOrConstructPhasic
(
    const objectRegistry& obr,
    const word& name,
    const dimensionSet& dims,
    const scalar& defaultValue
) const
{
    const word phasicName = phasePropertyName(name, phaseName_);

    // Try to look up the phasic one
    volScalarField* Tphasic =
        obr.objectRegistry::lookupObjectRefPtr<volScalarField>(phasicName);
    if (Tphasic)
    {
        return *Tphasic;
    }

    // Try to read the phasic one
    IOobject phasicHeader
    (
        phasicName,
        obr.time().timeName(),
        obr,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (phasicHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField* fPtr
        (
            new volScalarField(phasicHeader, mesh())
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);

        return *fPtr;
    }

    if (phaseName_ != word::null)
    {
        // Try to read the global one and initialise phasic from that
        IOobject header
        (
            name,
            obr.time().timeName(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (header.typeHeaderOk<volScalarField>(true))
        {
            volScalarField f0(header, mesh());

            // Make a copy and turn off writing when the field did not exist in
            // setup
            phasicHeader.readOpt() = IOobject::NO_READ;
            phasicHeader.writeOpt() = IOobject::NO_WRITE;
            volScalarField* fPtr(new volScalarField(phasicHeader, f0));

            // Transfer ownership of this object to the objectRegistry
            fPtr->store(fPtr);

            return *fPtr;
        }
    }

    // Finally, create a new field; don't write
    volScalarField* fPtr
    (
        new volScalarField
        (
            IOobject
            (
                phasicName,
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(name, dims, defaultValue),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);

    return *fPtr;
}


void Foam::basicThermo::lookupAndCheckout(const char* name) const
{
    if (db().foundObject<volScalarField>(name))
    {
         db().checkOut(*db()[name]);
    }
}


Foam::IOobject Foam::basicThermo::loadDictionary
(
    const objectRegistry& obr,
    const word& phaseName,
    const bool isRegistred
)
{
    IOobject header
    (
        phasePropertyName(dictName, phaseName),
        obr.time().caseConstant(),
        obr,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );
    IOobject matHeader
    (
        matDictName,
        obr.time().caseSystem(),
        obr,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );
    if
    (
        header.typeHeaderOk<dictionary>(true)
     && matHeader.typeHeaderOk<dictionary>(true)
    )
    {
        FatalErrorInFunction
            << "It isn't allowed to use both \"" << matDictName << "\""
            << " and \"thermophysicalProperties\" dictionaries."
            << nl << exit(FatalError);
    }
    if (header.typeHeaderOk<dictionary>(true))
    {
        if (dictName == matDictName)
        {
            FatalErrorInFunction
                << "It isn't allowed to use both \"" << matDictName << "\""
                << " and \"thermophysicalProperties\" dictionaries."
                << nl << exit(FatalError);
        }
        return
            IOobject
            (
                phasePropertyName(dictName, phaseName),
                obr.time().constant(),
                obr,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                isRegistred
            );
    }
    else if (matHeader.typeHeaderOk<dictionary>(true))
    {
        dictName = matDictName;
        return
            IOobject
            (
                dictName,
                obr.time().system(),
                obr,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                isRegistred
            );
    }
    FatalErrorInFunction
        << "Either \"" << matDictName << "\""
        << " or \"thermophysicalProperties\" is required."
        << nl << exit(FatalError);

    return IOobject(obr, "dummy");
}


Foam::basicThermo::basicThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    IOdictionary(loadDictionary(obr, phaseName, true)),
    phaseName_(phaseName),
    p_(lookupOrConstruct(obr, "p", dimPressure)),
    pRef_(readPref()),
    T_(lookupOrConstructPhasic(obr, "T", dimTemperature, 300)),
    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionSet(1, -1, -1, 0, 0)
    ),
    dpdt_(lookupOrDefault<Switch>("dpdt", true))
{
    if (dictName == matDictName)
    {
        (*this).rename(IOobject::groupName(matDictName, phaseName));
    }
}


Foam::basicThermo::basicThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            obr.time().constant(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    phaseName_(phaseName),
    p_(lookupOrConstruct(obr, "p", dimPressure)),
    pRef_(readPref()),
    T_(lookupOrConstructPhasic(obr, "T", dimTemperature)),
    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionSet(1, -1, -1, 0, 0)
    ),
    dpdt_(lookupOrDefault<Switch>("dpdt", true))
{
    if (dictName == matDictName)
    {
        (*this).rename(IOobject::groupName(matDictName, phaseName));
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicThermo> Foam::basicThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return New<basicThermo>(obr, phaseName);
}


Foam::basicThermo* Foam::basicThermo::lookupPtr
(
    const objectRegistry& obr, const word& phaseName
)
{
    basicThermo* matPtr =
        obr.lookupObjectRefPtr<basicThermo>(matDictName);
    if (matPtr)
    {
        return matPtr;
    }
    else
    {
        return
            obr.lookupObjectRefPtr<basicThermo>
            (
                IOobject::groupName(dictName, phaseName)
            );
    }
}


Foam::basicThermo& Foam::basicThermo::lookupOrCreate
(
    const objectRegistry& obr, const word& phaseName
)
{
    basicThermo* thermoPtr = lookupPtr(obr, phaseName);
    if (!thermoPtr)
    {
        if (basicThermo::debug)
        {
            Pout<< "basicThermo::lookupOrCreate(const objectRegistry&, const word&) : "
                << "constructing thermophysical model for region "
                << obr.name() << endl;
        }
        thermoPtr = New(obr, phaseName).ptr();
        regIOobject::store(thermoPtr);
    }

    return *thermoPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{
    lookupAndCheckout("p");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::basicThermo::readPref() const
{
    const word pRefName("pRef");
    if (found("referenceFields") && isDict("referenceFields"))
    {
        return
            subDict("referenceFields").lookupOrDefault<dimensionedScalar>
            (
                "p",
                dimensionedScalar(pRefName, dimPressure, 0)
            ).value();
    }
    else if (found(pRefName))
    {
        if (!isDict("thermoType"))
        {
            Warning
                << "Loading old style " << pRefName << " "  << nl
                << "(not supported in new material library)"
                << nl << endl;
        }
        return lookup<scalar>(pRefName);
    }
    return 0.0;
}


const Foam::basicThermo& Foam::basicThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    word phaseDictName(phasePropertyName(dictName, pf.internalField().group()));

    if (pf.db().foundObject<basicThermo>(phaseDictName))
    {
        return pf.db().lookupObject<basicThermo>(phaseDictName);
    }
    else
    {
        HashTable<const basicThermo*> thermos =
            pf.db().lookupClass<basicThermo>();

        for
        (
            HashTable<const basicThermo*>::iterator iter = thermos.begin();
            iter != thermos.end();
            ++iter
        )
        {
            if
            (
                &(iter()->he().internalField())
              == &(pf.internalField())
            )
            {
                return *iter();
            }
        }
    }

    return pf.db().lookupObject<basicThermo>(phaseDictName);
}


void Foam::basicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorInFunction
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << " and " << phasePropertyName(c)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c,
    const word& d
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
         || he().name() == phasePropertyName(d)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << ", " << phasePropertyName(c)
            << " and " << phasePropertyName(d)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


Foam::wordList Foam::basicThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList::null();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList::null();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}


Foam::volScalarField& Foam::basicThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basicThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicThermo::alpha() const
{
    return alpha_;
}


const Foam::scalarField& Foam::basicThermo::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


bool Foam::basicThermo::read()
{
    if (regIOobject::read())
    {
        pRef_ = readPref();
        return true;
    }
    else
    {
        return false;
    }

}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::basicThermo::operator[](const word& modelName) const
{
    if (dictName != matDictName)
    {
        FatalErrorInFunction
            << "This operator requires material properties library"
            << exit(FatalError);
    }
    materialTables& mat =
        db().subRegistry("materialModels").lookupObjectRef<materialTables>
        (
            "materialTables"
        );

    return mat(modelName)();
}


// ************************************************************************* //
