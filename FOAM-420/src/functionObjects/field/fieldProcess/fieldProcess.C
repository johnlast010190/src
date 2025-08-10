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
    (c) 2008 Icon-CG Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldProcess/fieldProcess.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldProcess, 0);
    addToRunTimeSelectionTable(functionObject, fieldProcess, dictionary);
}
}

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::functionObjects::fieldProcess::calcTypes, 43>::names[] =
{
    "none",
    "ptot",
    "Cp",
    "CpComponents",
    "pGauge",
    "Cg",
    "lambda2",
    "Q",
    "Gcriterion",
    "vorticity",
    "enstrophy",
    "wallShear",
    "yPlus",
    "mag",
    "grad",
    "magGrad",
    "plus",
    "minus",
    "multiply",
    "divide",
    "max",
    "min",
    "bound",
    "pow",
    "average",
    "dot",
    "cross",
    "Co",
    "div",
    "Pe",
    "Unw",
    "wallHTC",
    "nearWallValue",
    "ddt",
    "patchWrite",
    "strainRate",
    "Urel",
    "Urtz",
    "wallHeatFlux",
    "Mach",
    "skinFriction",
    "cylindricalTransform",
    "proudmanSources"
};

const Foam::NamedEnum<Foam::functionObjects::fieldProcess::calcTypes, 43>
    Foam::functionObjects::fieldProcess::calcTypeNames_;

template<>
const char* Foam::NamedEnum<Foam::functionObjects::fieldProcess::opFieldTypes, 6>::names[] =
{
    "none",
    "scalar",
    "vector",
    "tensor",
    "symmTensor",
    "sphTensor"
};

const Foam::NamedEnum<Foam::functionObjects::fieldProcess::opFieldTypes, 6>
    Foam::functionObjects::fieldProcess::opFieldTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldProcess::writeField(word fieldName)
{
    switch(volFieldType(fieldName))
    {
        case ftnone:
        {
            WarningInFunction
                << "volField: " << fieldName
                << " not found in registry." << endl;
        } break;

        case ftscalar:
        {
            (lookupObject<volScalarField>(fieldName)).write();
        } break;

        case ftvector:
        {
            (lookupObject<volVectorField>(fieldName)).write();
        } break;

        case fttensor:
        {
            (lookupObject<volTensorField>(fieldName)).write();
        } break;

        case ftsymmTensor:
        {
            (lookupObject<volSymmTensorField>(fieldName)).write();
        } break;

        case ftsphTensor:
        {
            (lookupObject<volSphericalTensorField>(fieldName)).write();
        }
    }
}

void Foam::functionObjects::fieldProcess::clearField(word fieldName)
{
    switch(volFieldType(fieldName))
    {
        case ftnone:
        {
            //If field not found it has been deleted by object registry
            //WarningInFunction
            //    << "volField: " << fieldName
            //    << " not found in registry." << endl;
        } break;

        case ftscalar:
        {
            autoPtr<volScalarField> field
            (
                &const_cast<volScalarField&>
                (
                    lookupObject<volScalarField>(fieldName)
                )
            );

            field->release();


        } break;

        case ftvector:
        {
            autoPtr<volVectorField> field
            (
                &const_cast<volVectorField&>
                (
                    lookupObject<volVectorField>(fieldName)
                )
            );

            field->release();

        } break;

        case fttensor:
        {
            autoPtr<volTensorField> field
            (
                &const_cast<volTensorField&>
                (
                    lookupObject<volTensorField>(fieldName)
                )
            );

            field->release();

        } break;

        case ftsymmTensor:
        {
            autoPtr<volSymmTensorField> field
            (
                &const_cast<volSymmTensorField&>
                (
                    lookupObject<volSymmTensorField>(fieldName)
                )
            );

            field->release();

        } break;

        case ftsphTensor:
        {
            autoPtr<volSphericalTensorField> field
            (
                &const_cast<volSphericalTensorField&>
                (
                    lookupObject<volSphericalTensorField>(fieldName)
                )
            );

            field->release();

        }
    }

}

void Foam::functionObjects::fieldProcess::clearFields()
{
    forAll(resultFieldNames_, fnI)
    {
        clearField(resultFieldNames_[fnI]);
    }

    resultFieldNames_.clear();
    writeStatus_.clear();
}

Foam::functionObjects::fieldProcess::opFieldTypes
Foam::functionObjects::fieldProcess::volFieldType(word field)
{

    if (foundObject<volScalarField>(field))
    {
        return ftscalar;
    }
    else if (foundObject<volVectorField>(field))
    {
        return ftvector;
    }
    else if (foundObject<volTensorField>(field))
    {
        return fttensor;
    }
    else if (foundObject<volSymmTensorField>(field))
    {
        return ftsymmTensor;
    }
    else if (foundObject<volSphericalTensorField>(field))
    {
        return ftsphTensor;
    }
    else
    {
        return ftnone;
    }

    return ftnone;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldProcess::fieldProcess
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    operations_(dict.lookup("operations")),
    resultFieldNames_(100),
    writeStatus_(100),
    checkFields_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldProcess::~fieldProcess()
{
    clearFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldProcess::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    calculate();

    return true;
}


void Foam::functionObjects::fieldProcess::writeFields()
{
    forAll(resultFieldNames_, rfnI)
    {
        if (writeStatus_[rfnI])
        {
            writeField(resultFieldNames_[rfnI]);
        }
    }
}


void Foam::functionObjects::fieldProcess::calculate()
{
    forAll(operations_, opI)
    {
        const dictionary& opDict(operations_[opI]);
        word operation(opDict.lookup("operation"));
        switch(calcTypeNames_.read(IStringStream(operation)()))
        {
            case ctnone:
                none(opDict);
            break;
            case ctptot:
                ptot(opDict);
            break;
            case ctCp:
                Cp(opDict);
            break;
            case ctCpComponents:
                CpComponents(opDict);
            break;
            case ctpGauge:
                pGauge(opDict);
            break;
            case ctCg:
                Cg(opDict);
            break;
            case ctlambda2:
                lambda2(opDict);
            break;
            case ctQ:
                Q(opDict);
            break;
            case ctGcriterion:
                Gcriterion(opDict);
            break;
            case ctvorticity:
                vorticity(opDict);
            break;
            case ctenstrophy:
                enstrophy(opDict);
            break;
            case ctwallShear:
                wallShear(opDict);
            break;
            case ctyPlus:
                yPlus(opDict);
            break;
            case ctmag:
                magOp(opDict);
            break;
            case ctgrad:
                gradOp(opDict);
            break;
            case ctmaggrad:
                magGradOp(opDict);
            break;
            case ctplus:
                notImplementedYet(operation);
            break;
            case ctminus:
                notImplementedYet(operation);
            break;
            case ctmultiply:
                multiplyOp(opDict);
            break;
            case ctdivide:
                notImplementedYet(operation);
            break;
            case ctmax:
                maxOp(opDict);
            break;
            case ctmin:
                minOp(opDict);
            break;
            case ctbound:
                notImplementedYet(operation);
            break;
            case ctpow:
                powOp(opDict);
            break;
            case ctaverage:
                notImplementedYet(operation);
            break;
            case ctdot:
                dotOp(opDict);
            break;
            case ctcross:
                crossOp(opDict);
            break;
            case ctCo:
                Co(opDict);
            break;
            case ctdiv:
                divOp(opDict);
            break;
            case ctPe:
                Pe(opDict);
            break;
            case ctUnw:
                Unw(opDict);
            break;
            case ctHtc:
                wallHTC(opDict);
            break;
            case ctNwf:
                nearWallValue(opDict);
            break;
            case ctDdt:
                ddtOp(opDict);
            break;
            case ctPatchWrite:
                patchWrite(opDict);
            break;
            case ctStrainRate:
                strainRate(opDict);
            break;
            case ctUrel:
                Urel(opDict);
            break;
            case ctUrtz:
                Urtz(opDict);
            break;
            case ctHeatFlux:
                wallHeatFlux(opDict);
            break;
            case ctMach:
                Mach(opDict);
            break;
            case ctSkinFriction:
                skinFriction(opDict);
            break;
            case ctCylindricalTransform:
                cylindricalTransform(opDict);
            break;
            case ctProudmanSources:
                proudmanSources(opDict);
            break;
            default:
                WarningInFunction
                    << operation << " is an invalid operation." << endl;
            break;
        }
    }

    resultFieldNames_.shrink();
    writeStatus_.shrink();
}

bool Foam::functionObjects::fieldProcess::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    writeFields();

    Info<< endl;

    return true;
}


bool Foam::functionObjects::fieldProcess::end()
{
    return true;
}


bool Foam::functionObjects::fieldProcess::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    operations_.clear();
    operations_ = PtrList<dictionary>(dict.lookup("operations"));

    clearFields();

    //Report any duplicate fields
    checkFields_ = true;
    calculate();
    checkFields_ = false;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
