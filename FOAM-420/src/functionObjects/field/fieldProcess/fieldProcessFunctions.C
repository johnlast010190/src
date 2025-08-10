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
    (c) 2008 Icon-CG Ltd.
    (c) 2010-2018, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "fieldProcess/fieldProcess.H"
#include "fvMesh/wallDist/nearWallDist/nearWallDist.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "primitives/strings/wordRes/wordRes.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "coordinate/systems/cylindricalCS.H"
#include "sources/derived/MRFSource/MRFSource.H"
#include "solidThermo/solidThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
functionObjects::fieldProcess::createNearWallField
(
    const GeometricField<Type, fvPatchField, volMesh>& inField
)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    tmp<GeometricField<Type, fvPatchField, volMesh>> nwf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "nearWallValueField",
                obr_.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            ("nearWallValueField", inField.dimensions(), pTraits<Type>::zero)
        )
    );
    forAll(nwf->boundaryField(), i)
    {
        if (!nwf->boundaryField()[i].coupled())
        {
            nwf->boundaryFieldRef()[i] =
                inField.boundaryField()[i].patchInternalField()();
        }
    }

    return nwf;
}


template<class geoField>
void functionObjects::fieldProcess::setWallsToNearCell(geoField& gf)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];

        if (isType<wallFvPatch>(p))
        {
            gf.boundaryFieldRef()[patchi].forceAssign
            (
                gf.boundaryField()[patchi].patchInternalField()
            );
        }
    }
}


template<class geoField>
void functionObjects::fieldProcess::setWallsToNearCell(autoPtr<geoField>& gf)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];

        if (isType<wallFvPatch>(p))
        {
            gf->boundaryFieldRef()[patchi].forceAssign
            (
                gf->boundaryField()[patchi].patchInternalField()
            );
        }
    }
}


dimensionedScalar functionObjects::fieldProcess::readConstRho(word rhoName)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    IOdictionary tp
    (
        IOobject
        (
            "transportProperties",
            obr_.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    return dimensionedScalar("rho",tp.lookup(rhoName));
}


void functionObjects::fieldProcess::notImplementedYet(word nameOfFunction)
{
   WarningInFunction
        << "Operation " << nameOfFunction << " has not been implemented."
        << endl;

}

const tmp<volScalarField> functionObjects::fieldProcess::getMuEff()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);

        return model.muEff();
    }
    else if (foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<icoTurbModel>(turbulenceModel::propertiesName);

        return model.rho()*model.nuEff();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& model =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu",model.lookup("nu"));
        dimensionedScalar rho("rho",model.lookup("rho"));

        const auto& mesh = refCast<const fvMesh>(obr_);

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "muEff",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                mesh,
                rho*nu
             )
        );
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for getmuEff"
            << exit(FatalError);
        ::abort();
    }
}


const tmp<volScalarField> functionObjects::fieldProcess::getK()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);

        return model.k();
    }
    else if (foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<icoTurbModel>(turbulenceModel::propertiesName);

        return model.k();
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for k"
            << exit(FatalError);
        ::abort();
    }
}


const tmp<volScalarField> functionObjects::fieldProcess::getEpsilon()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);

        return model.epsilon();
    }
    else if (foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<icoTurbModel>(turbulenceModel::propertiesName);

        return model.epsilon();
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for epsilon"
            << exit(FatalError);
        ::abort();
    }
}


const tmp<volVectorField> functionObjects::fieldProcess::wallShear()
{
    const auto& mesh = refCast<const fvMesh>(obr_);

    tmp<volVectorField> wallShear
    (
        new volVectorField
        (
            IOobject
            (
                "wallShear",
                obr_.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "wallshear",
                dimensionSet(1, -1, -2, 0, 0, 0, 0),
                Zero
            ),
            fixedValueFvPatchScalarField::typeName
        )
    );

    const tmp<volScalarField> muEff = getMuEff();

    forAll(mesh.boundary(), pI)
    {
        if  (isA<wallFvPatch>(mesh.boundary()[pI]))
        {
            const fvPatch& cPatch = mesh.boundary()[pI];
            const fvPatchVectorField& Uboundary =
                cPatch.lookupPatchField<volVectorField, vector>("U");

            vectorField Up( Uboundary.patchInternalField() - Uboundary );
            Up -= cPatch.Sf()*(cPatch.Sf() & Up)/ magSqr(cPatch.Sf());

            const scalarField& ry = cPatch.deltaCoeffs();

            wallShear->boundaryFieldRef()[pI].forceAssign
            (
                muEff->boundaryField()[pI] * Up * ry
            );
        }
    }

    return wallShear;
}


template<class geoField>
void functionObjects::fieldProcess::setField
(
    word fieldName,
    geoField* inFieldPtr,
    const dictionary& dict
)
{
    tmp<geoField> rinField(inFieldPtr);
    Switch writeStatus = dict.lookupOrDefault<Switch>("write", true);

    Switch nearCellValue = false;
    if (dict.found("nearCellValue"))
    {
        nearCellValue = readBool(dict.lookup("nearCellValue"));
    }
    if (nearCellValue)
    {
        setWallsToNearCell(rinField.ref());
    }

    if (volFieldType(fieldName) == ftnone)
    {
        rinField->regIOobject::rename(fieldName);
        rinField->regIOobject::readOpt() = IOobject::NO_READ;
        rinField->regIOobject::writeOpt() = IOobject::NO_WRITE;

        store(fieldName, rinField);
        resultFieldNames_.append(fieldName);
        writeStatus_.append(bool(writeStatus));
    }
    else if (lookupObject<geoField>(fieldName).ownedByRegistry())
    {
        // The field already exists. Check type, if identical to the
        // current case, overwrite, if not provide warning and terminate op
        // remember to place field in resultFieldNames if not already present
        if
        (
            (
                rinField->type() == volScalarField::typeName
             && volFieldType(fieldName) == ftscalar
            )
         || (
                rinField->type() == volVectorField::typeName
             && volFieldType(fieldName) == ftvector
            )
         || (
                rinField->type() == volTensorField::typeName
             && volFieldType(fieldName) == fttensor
            )
         || (
                rinField->type() == volSymmTensorField::typeName
             && volFieldType(fieldName) == ftsymmTensor
            )
         || (
                rinField->type() == volSphericalTensorField::typeName
             && volFieldType(fieldName) == ftsphTensor
            )
        )
        {
            geoField& dbField =
                const_cast<geoField&>(lookupObject<geoField>(fieldName));

            dbField.dimensions().reset(rinField->dimensions());
            dbField = rinField();
            dbField.boundaryFieldRef().forceAssign(rinField().boundaryField());

            label resultIndex = -1;
            forAll(resultFieldNames_, rfI)
            {
                if (resultFieldNames_[rfI] == fieldName)
                {
                    resultIndex = rfI;
                    writeStatus_[rfI] = writeStatus;
                    break;
                }
            }

            if (resultIndex == -1)
            {
                resultFieldNames_.append(fieldName);
                writeStatus_.append(bool(writeStatus));
            }
            else if (checkFields_)
            {
                Info<< "Field " << fieldName << " already exists. "<< endl;
            }
        }
        else
        {
            WarningInFunction
                << fieldName
                << " already in registry and of different type to "
                << "operation result. "
                << fieldName << " will not be calculated."
                << endl;

            fieldName = "";
        }
    }
    else
    {
        FatalErrorInFunction
            << "ID CONFLICT:" << fieldName
            << " already in registry and owned by a non-registry object. "
            << "Change field name " << fieldName
            << " to a non-conflicting value."
            << exit(FatalError);

        fieldName = "";
    }

}


template<class geoField>
void functionObjects::fieldProcess::setField
(
    word fieldName,
    geoField& rinField,
    const dictionary& dict
)
{
    Switch writeStatus = true;
    if (dict.found("write"))
    {
        writeStatus = readBool(dict.lookup("write"));
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    autoPtr<geoField> inField
    (
        new geoField
        (
            IOobject
            (
                rinField.name(),
                obr_.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rinField
        )
    );

    Switch nearCellValue = false;
    if (dict.found("nearCellValue"))
    {
        nearCellValue = readBool(dict.lookup("nearCellValue"));
    }

    if (nearCellValue)
    {
        setWallsToNearCell(inField);
    }

    if (volFieldType(fieldName) == ftnone)
    {
        inField->regIOobject::rename(fieldName);

        inField->store(inField);
        resultFieldNames_.append(fieldName);
        writeStatus_.append(bool(writeStatus));
    }
    else if (lookupObject<geoField>(fieldName).ownedByRegistry())
    {
        // The field already exists. Check type, if identical to the
        // current case, overwrite, if not provide warning and terminate op
        // remember to place field in resultFieldNames if not already present
        if
        (
            (
                inField->type() == volScalarField::typeName
             && volFieldType(fieldName) == ftscalar
            )
         || (
                inField->type() == volVectorField::typeName
             && volFieldType(fieldName) == ftvector
            )
         || (
                inField->type() == volTensorField::typeName
             && volFieldType(fieldName) == fttensor
            )
         || (
                inField->type() == volSymmTensorField::typeName
             && volFieldType(fieldName) == ftsymmTensor
            )
         || (
                inField->type() == volSphericalTensorField::typeName
             && volFieldType(fieldName) == ftsphTensor
            )
        )
        {
            geoField& dbField =
                const_cast<geoField&>(lookupObject<geoField>(fieldName));

            dbField.dimensions().reset(inField->dimensions());
            dbField = inField();
            dbField.boundaryFieldRef().forceAssign(inField().boundaryField());

            label resultIndex = -1;
            forAll(resultFieldNames_, rfI)
            {
                if (resultFieldNames_[rfI] == fieldName)
                {
                    resultIndex = rfI;
                    writeStatus_[rfI] = writeStatus;
                    break;
                }
            }

            if (resultIndex == -1)
            {
                resultFieldNames_.append(fieldName);
                writeStatus_.append(bool(writeStatus));
            }
            else if (checkFields_)
            {
                Info<< "Field " << fieldName << " already exists. "<< endl;
            }
        }
        else
        {
            WarningInFunction
                << fieldName
                << " already in registry and of different type to "
                << "operation result. "
                << fieldName << " will not be calculated."
                << endl;
            fieldName = "";
        }
    }
    else
    {
        FatalErrorInFunction
            << "ID CONFLICT:" << fieldName
            << " already in registry and owned by a non-registry object. "
            << "Change field name " << fieldName
            << " to a non-conflicting value."
            << exit(FatalError);

        fieldName = "";
    }
}


void functionObjects::fieldProcess::none(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    switch(inFieldType)
    {
        case ftnone:
        {
            WarningInFunction
                << "Input field " << inField << " does not exist."
                << nl << "Skipping calculation of " << fieldName << endl;
        } break;

        case ftscalar:
        {
            setField
            (
                fieldName,
                lookupObjectRef<volScalarField>(inField),
                dict
            );
        } break;

        case ftvector:
        {
            setField
            (
                fieldName,
                lookupObjectRef<volVectorField>(inField),
                dict
            );
        } break;

        case fttensor:
        {
            setField
            (
                fieldName,
                lookupObjectRef<volTensorField>(inField),
                dict
            );
        } break;

        case ftsymmTensor:
        {
            setField
            (
                fieldName,
                lookupObjectRef<volSymmTensorField>(inField),
                dict
            );
        } break;

        case ftsphTensor:
        {
            setField
            (
                fieldName,
                lookupObjectRef<volSphericalTensorField>(inField),
                dict
            );
        }
    }
}


void functionObjects::fieldProcess::ptot(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");
    const word pName = dict.lookupOrDefault<word>("p", "p");

    if
    (
        foundObject<volScalarField>(pName)
     && foundObject<volVectorField>(Uname)
    )
    {

        const auto& pt = lookupObject<volScalarField>(pName);
        const auto& Ut = lookupObject<volVectorField>(Uname);

        word rhoName = "rho";
        if (dict.found("rho"))
        {
            rhoName = dict.lookup<word>("rho");
        }

        if (pt.dimensions() == dimensionSet(0,2,-2,0,0))
        {
            dimensionedScalar rhoIncomp = readConstRho(rhoName);

            setField
            (
                fieldName,
                (rhoIncomp*pt + 0.5*rhoIncomp*magSqr(Ut)).ptr(),
                dict
            );
        }
        else if (pt.dimensions() == dimensionSet(1,-1,-2,0,0))
        {
            const auto& rhoComp =
                lookupObject<volScalarField>(rhoName);

            setField
            (
                fieldName,
                (pt + 0.5*rhoComp*magSqr(Ut)).ptr(),
                dict
            );
        }
        else
        {
            WarningInFunction
                << pName << " does not have pressure dimensions."
                << nl << "Skipping calculation of ptot." << endl;

        }
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname << " and/or " << pName
            << " in the database." << nl
            << "Skipping calculation of ptot." << endl;
    }
}


void functionObjects::fieldProcess::Cp(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word pName = dict.lookupOrDefault<word>("p", "p");

    if (foundObject<volScalarField>(pName))
    {
        const auto& pt = lookupObject<volScalarField>(pName);
        scalar Uref = readScalar(dict.lookup("Uref"));
        dimensionedScalar U2("U2", pt.dimensions(), 0.5*magSqr(Uref));

        if (pt.dimensions() == dimensionSet(0, 2, -2, 0, 0))
        {
            dimensionedScalar PRef
            (
                "Pref",
                pt.dimensions(),
                dict.lookupOrDefault<scalar>("pRef", 0.0)
            );
            setField
            (
                fieldName,
                ((pt - PRef)/U2).ptr(),
                dict
            );
        }
        else if (pt.dimensions() == dimensionSet(1,-1,-2,0,0))
        {
            scalar pRefValue = 0;
            if (dict.isDict("referenceFields"))
            {
                const dictionary& refDict = dict.subDict("referenceFields");
                pRefValue =
                    refDict.found("p")
                  ? refDict.lookup<dimensionedScalar>("p").value()
                  : 0;
            }
            dimensionedScalar PRef
            (
                "Pref",
                pt.dimensions(),
                dict.isDict("referenceFields")
              ? pRefValue
              : readScalar(dict.lookup("pRef"))
            );
            dimensionedScalar rhoRef
            (
                "rhoRef",
                dimDensity,
                readScalar(dict.lookup("rhoRef"))
            );

            setField
            (
                fieldName,
                ((pt - PRef)/(rhoRef.value()*U2)).ptr(),
                dict
            );
        }
        else
        {
            WarningInFunction
                << pName << " does not have pressure dimensions."
                << nl << "Skipping calculation of Cp." << endl;

        }

    }
    else
    {
        WarningInFunction
            << "Could not find " << pName
            << " in the database." << nl
            << "Skipping calculation of Cp." << endl;
    }

}


void functionObjects::fieldProcess::CpComponents(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word pName = dict.lookupOrDefault<word>("p", "p");

    if (foundObject<volScalarField>(pName))
    {
        const auto& pt = lookupObject<volScalarField>(pName);
        scalar Uref = readScalar(dict.lookup("Uref"));
        dimensionedScalar U2("U2", pt.dimensions(), 0.5*magSqr(Uref));

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        tmp<volVectorField> patchNormals
        (
            new volVectorField
            (
                IOobject
                (
                    "patchNormals",
                    obr_.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("patchNormals", dimless, Zero)
            )
        );
        forAll(mesh.boundary(), pI)
        {
            if  (!mesh.boundary()[pI].coupled())
            {
                const fvPatch& cPatch = mesh.boundary()[pI];
                patchNormals->boundaryFieldRef()[pI].forceAssign(cPatch.nf());
            }
        }


        if (pt.dimensions() == dimensionSet(0,2,-2,0,0))
        {
            dimensionedScalar PRef
            (
                "Pref",
                pt.dimensions(),
                dict.lookupOrDefault<scalar>("pRef", 0.0)
            );
            setField
            (
                fieldName,
                ((pt - PRef)/U2*patchNormals).ptr(),
                dict
            );
        }
        else if (pt.dimensions() == dimensionSet(1,-1,-2,0,0))
        {
            scalar pRefValue = 0;
            if (dict.isDict("referenceFields"))
            {
                const dictionary& refDict = dict.subDict("referenceFields");
                pRefValue =
                    refDict.found("p")
                  ? refDict.lookup<dimensionedScalar>("p").value()
                  : 0;
            }
            dimensionedScalar PRef
            (
                "Pref",
                pt.dimensions(),
                dict.isDict("referenceFields")
              ? pRefValue
              : readScalar(dict.lookup("pRef"))
            );
            dimensionedScalar rhoRef
            (
                "rhoRef",
                dimDensity,
                readScalar(dict.lookup("rhoRef"))
            );

            setField
            (
                fieldName,
                ((pt - PRef)/(rhoRef*U2)*patchNormals).ptr(),
                dict
            );
        }
        else
        {
            WarningInFunction
                << pName << " does not have pressure dimensions."
                << nl << "Skipping calculation of Cp." << endl;
        }
    }
    else
    {
        WarningInFunction
            << "Could not find " << pName
            << " in the database." << nl
            << "Skipping calculation of Cp." << endl;
    }
}


void functionObjects::fieldProcess::pGauge(const dictionary& dict)
{
    const word fieldName = dict.lookupOrDefault<word>("fieldName", "pGauge");
    const word pName = dict.lookupOrDefault<word>("p", "p");

    if (foundObject<volScalarField>(pName))
    {
        const auto& pt = lookupObject<volScalarField>(pName);
        scalar pRefValue = 0;
        if (dict.isDict("referenceFields"))
        {
            const dictionary& refDict = dict.subDict("referenceFields");
            pRefValue =
                refDict.found("p")
              ? refDict.lookup<dimensionedScalar>("p").value()
              : 101325;
        }
        dimensionedScalar pRef
        (
            "pRef",
            pt.dimensions(),
            dict.isDict("referenceFields")
          ? pRefValue
          : dict.lookupOrDefault<scalar>("pRef", 101325)
        );

        setField(fieldName, (pt - pRef).ptr(), dict);
    }
    else
    {
        WarningInFunction
            << "Could not find " << pName
            << " in the database." << nl
            << "Skipping calculation of pGauge." << endl;
    }

}


void functionObjects::fieldProcess::Cg(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word pName = dict.lookupOrDefault<word>("p", "p");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if
    (
        foundObject<volScalarField>(pName)
     && foundObject<volVectorField>(Uname)
    )
    {
        const auto& pt = lookupObject<volScalarField>(pName);
        const auto& Ut = lookupObject<volVectorField>(Uname);
        scalar Uref = readScalar(dict.lookup("Uref"));
        dimensionedScalar U2("U2", dimVelocity*dimVelocity, 0.5*magSqr(Uref));


        if (pt.dimensions() == dimensionSet(0,2,-2,0,0))
        {
            setField
            (
                fieldName,
                ((pt + 0.5*magSqr(Ut))/U2).ptr(),
                dict
            );
        }
        else if (pt.dimensions() == dimensionSet(1,-1,-2,0,0))
        {
            word rhoRefName = "rho";
            if (dict.found("rho"))
            {
                rhoRefName = dict.lookup<word>("rho");
            }
            const auto& rhoComp = lookupObject<volScalarField>(rhoRefName);

            dimensionedScalar rhoRef
            (
                "rhoRef",
                dimDensity,
                readScalar(dict.lookup("rhoRef"))
            );
            scalar pRefValue = 0;
            if (dict.isDict("referenceFields"))
            {
                const dictionary& refDict = dict.subDict("referenceFields");
                pRefValue =
                    refDict.found("p")
                  ? refDict.lookup<dimensionedScalar>("p").value()
                  : 101325;
            }
            dimensionedScalar pRef
            (
                "pRef",
                pt.dimensions(),
                dict.isDict("referenceFields")
              ? pRefValue
              : readScalar(dict.lookup("pRef"))
            );

            setField
            (
                fieldName,
                ((pt - pRef + 0.5*rhoComp*magSqr(Ut))/(rhoRef*U2)).ptr(),
                dict
            );
        }
        else
        {
            WarningInFunction
                << pName << " does not have pressure dimensions."
                << nl << "Skipping calculation of ptot." << endl;

        }
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname << " and/or " << pName
            << " in the database." << nl
            << "Skipping calculation of Cg." << endl;
    }
}


void functionObjects::fieldProcess::lambda2(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (foundObject<volVectorField>(Uname))
    {
        const auto& Ut = lookupObject<volVectorField>(Uname);
        const volTensorField gradU(fvc::grad(Ut));

        volTensorField SSplusWW
        (
            (symm(gradU) & symm(gradU))
          + (skew(gradU) & skew(gradU))
        );

        setField
        (
            fieldName,
            (-eigenValues(SSplusWW)().component(vector::Y)).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of lambda2." << endl;
    }
}


void functionObjects::fieldProcess::Q(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (foundObject<volVectorField>(Uname))
    {
        const auto& Ut = lookupObject<volVectorField>(Uname);
        const volTensorField gradU(fvc::grad(Ut));

        setField
        (
            fieldName,
            (0.5*(sqr(tr(gradU)) - tr(((gradU)&(gradU))))).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of Q." << endl;
    }
}


void functionObjects::fieldProcess::Gcriterion(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (foundObject<volVectorField>(Uname))
    {
        const auto& Ut = lookupObject<volVectorField>(Uname);
        volSymmTensorField tsGradU(twoSymm(fvc::grad(Ut)));

        // Set diagonal to zero; Note: no usable operator currently exists
        forAll(tsGradU, cI)
        {
            tsGradU[cI].xx() = 0;
            tsGradU[cI].yy() = 0;
            tsGradU[cI].zz() = 0;
        }
        forAll(tsGradU.boundaryField(), pI)
        {
            if (!tsGradU.boundaryField()[pI].coupled())
            {
                forAll(tsGradU.boundaryField()[pI], fI)
                {
                    tsGradU.boundaryFieldRef()[pI][fI].xx() = 0;
                    tsGradU.boundaryFieldRef()[pI][fI].yy() = 0;
                    tsGradU.boundaryFieldRef()[pI][fI].zz() = 0;
                }
            }
        }

        setField(fieldName, (sqrt(0.5*(tsGradU && tsGradU))).ptr(), dict);
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of Gcriterion." << endl;
    }
}


void functionObjects::fieldProcess::vorticity(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (foundObject<volVectorField>(Uname))
    {
        const auto& Ut = lookupObject<volVectorField>(Uname);
        setField(fieldName, (fvc::curl(Ut)).ptr(), dict);
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of lambda2." << endl;
    }
}


void functionObjects::fieldProcess::enstrophy(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if
    (
        foundObject<volVectorField>(Uname)
    )
    {
        const auto& Ut = lookupObject<volVectorField>(Uname);

        setField
        (
            fieldName,
            (0.5*magSqr(fvc::curl(Ut))).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of lambda2." << endl;
    }
}


void functionObjects::fieldProcess::Co(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");

    word phiName = "phi";
    if (dict.found("phi"))
    {
        phiName = dict.lookup<word>("phi");
    }

    if (foundObject<surfaceScalarField>(phiName))
    {
        const auto& phit =
            lookupObject<surfaceScalarField>(phiName);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (phit.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {

            setField
            (
                fieldName,
                (
                    mag
                    (
                        fvc::reconstruct
                        (
                            phit *obr_.time().deltaT()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                        )
                    )
                ).ptr(),
                dict
            );
        }
        else if (phit.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            word rhoRefName = "rho";
            if (dict.found("rho"))
            {
                rhoRefName = dict.lookup<word>("rho");
            }
            const auto& rhoComp =
                lookupObject<volScalarField>(rhoRefName);

            setField
            (
                fieldName,
                (
                    mag
                    (
                        fvc::reconstruct
                        (
                            phit*obr_.time().deltaT()
                            /fvc::interpolate(rhoComp)
                          * mesh.surfaceInterpolation::deltaCoeffs()
                        )
                    )
                ).ptr(),
                dict
            );
        }
        else
        {
            WarningInFunction
                << "Incorrect dimensions of phi: " << phit.dimensions()
                << "Skipping calculation of Co." << endl;
        }


    }
    else
    {
        WarningInFunction
            << "Could not find " << phiName
            << " in the database." << nl
            << "Skipping calculation of Co." << endl;
    }
}


void functionObjects::fieldProcess::divOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    if (inFieldType != ftnone && inFieldType != ftscalar)
    {
        switch(inFieldType)
        {
            case ftvector:
            {
                setField
                (
                    fieldName,
                    (
                        fvc::div(lookupObject<volVectorField>(inField))
                    ).ptr(),
                    dict
                );
            } break;

            case fttensor:
            {
                setField
                (
                    fieldName,
                    (
                        fvc::div(lookupObject<volTensorField>(inField))
                    ).ptr(),
                    dict
                );
            } break;

            case ftsymmTensor:
            {
                setField
                (
                    fieldName,
                    (
                        fvc::div(lookupObject<volSymmTensorField>(inField))
                    ).ptr(),
                    dict
                );
            } break;

            case ftsphTensor:
            {
                setField
                (
                    fieldName,
                    (
                        fvc::div
                        (
                            lookupObject<volSphericalTensorField>(inField)
                        )
                    ).ptr(),
                    dict
                );
            } break;

            default:
                FatalErrorInFunction
                    << "Unsupported field type encountered."
                    << exit(FatalError);
            break;

        }
    }
    else
    {
        WarningInFunction
            << "Scalar or NULL input field for div operator with output  "
            << fieldName << nl
            << "Skipping calculation of div." << endl;
    }
}


void functionObjects::fieldProcess::magOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    if (inFieldType != ftnone)
    {
        switch(inFieldType)
        {
            case ftscalar:
            {
                setField
                (
                    fieldName,
                    (mag(lookupObject<volScalarField>(inField))).ptr(),
                    dict
                );
            } break;

            case ftvector:
            {
                setField
                (
                    fieldName,
                    (mag(lookupObject<volVectorField>(inField))).ptr(),
                    dict
                );
            } break;

            case fttensor:
            {
                setField
                (
                    fieldName,
                    (mag(lookupObject<volTensorField>(inField))).ptr(),
                    dict
                );
            } break;

            case ftsymmTensor:
            {
                setField
                (
                    fieldName,
                    (mag(lookupObject<volSymmTensorField>(inField))).ptr(),
                    dict
                );
            } break;

            case ftsphTensor:
            {
                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            lookupObject<volSphericalTensorField>(inField)
                        )
                    ).ptr(),
                    dict
                );
            } break;

            default:
            break;
        }
    }
    else
    {
        WarningInFunction
            << "NULL input field for mag operator with output  "
            << fieldName << nl
            << "Skipping calculation of mag." << endl;
    }
}


void functionObjects::fieldProcess::gradOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    if (inFieldType == ftscalar || inFieldType == ftvector)
    {
        switch(inFieldType)
        {
            case ftscalar:
            {
                setField
                (
                    fieldName,
                    (
                        fvc::grad(lookupObject<volScalarField>(inField))
                    ).ptr(),
                    dict
                );
            } break;

            case ftvector:
            {
                setField
                (
                    fieldName,
                    (
                        fvc::grad(lookupObject<volVectorField>(inField))
                    ).ptr(),
                    dict
                );
            } break;

            default:
            break;
        }
    }
    else
    {
        WarningInFunction
            << "NULL or tensor input field for grad operator with output  "
            << fieldName << nl
            << "Skipping calculation of grad." << endl;
    }
}


void functionObjects::fieldProcess::magGradOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    if (inFieldType == ftscalar || inFieldType == ftvector)
    {
        switch(inFieldType)
        {
            case ftscalar:
            {
                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            fvc::grad
                            (
                                lookupObject<volScalarField>(inField)
                            )
                        )
                    ).ptr(),
                    dict
                );
            } break;

            case ftvector:
            {
                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            fvc::grad
                            (
                                lookupObject<volVectorField>(inField)
                            )
                        )
                    ).ptr(),
                    dict
                );
            } break;

            default:
            break;
        }
    }
    else
    {
        WarningInFunction
            << "NULL or tensor input field for grad operator with output  "
            << fieldName << nl
            << "Skipping calculation of magGrad." << endl;
    }
}


void functionObjects::fieldProcess::multiplyOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField1 = dict.lookup<word>("inField1");
    word inField2 = dict.lookup<word>("inField2");
    opFieldTypes inField1Type = volFieldType(inField1);
    opFieldTypes inField2Type = volFieldType(inField2);

    if (inField1Type == ftscalar && inField2Type == ftscalar)
    {
        setField
        (
            fieldName,
            (
                lookupObject<volScalarField>(inField1)
              * lookupObject<volScalarField>(inField2)
            ).ptr(),
            dict
        );
    }
    else if (inField1Type == ftvector && inField2Type == ftvector)
    {
        setField
        (
            fieldName,
            (
                lookupObject<volVectorField>(inField1)
              & lookupObject<volVectorField>(inField2)
            ).ptr(),
            dict
        );
    }
    else if (inField1Type == ftscalar && inField2Type == ftvector)
    {
        setField
        (
            fieldName,
            (
                lookupObject<volScalarField>(inField1)
              * lookupObject<volVectorField>(inField2)
            ).ptr(),
            dict
        );
    }
    else if (inField1Type == ftvector && inField2Type == ftscalar)
    {
        setField
        (
            fieldName,
            (
                lookupObject<volVectorField>(inField1)
              * lookupObject<volScalarField>(inField2)
            ).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "NULL or tensor input field for multiply operator with output  "
            << fieldName << nl
            << "Skipping calculation of multiply." << endl;
    }
}


void functionObjects::fieldProcess::maxOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField1 = dict.lookup<word>("inField1");
    word inField2 = dict.lookup<word>("inField2");
    opFieldTypes inField1Type = volFieldType(inField1);
    opFieldTypes inField2Type = volFieldType(inField2);

    if (inField1Type == ftscalar && inField2Type == ftscalar)
    {
        setField
        (
            fieldName,
            (
                max
                (
                    lookupObject<volScalarField>(inField1),
                    lookupObject<volScalarField>(inField2)
                )
            ).ptr(),
            dict
        );
    }
    else if (inField1Type == ftvector && inField2Type == ftvector)
    {
        setField
        (
            fieldName,
            (
                max
                (
                    lookupObject<volVectorField>(inField1),
                    lookupObject<volVectorField>(inField2)
                )
            ).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Non-equal input field types for max operator with output  "
            << fieldName << nl
            << "Skipping calculation of max." << endl;
    }
}


void functionObjects::fieldProcess::minOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField1 = dict.lookup<word>("inField1");
    word inField2 = dict.lookup<word>("inField2");
    opFieldTypes inField1Type = volFieldType(inField1);
    opFieldTypes inField2Type = volFieldType(inField2);

    if (inField1Type == ftscalar && inField2Type == ftscalar)
    {
        setField
        (
            fieldName,
            (
                min
                (
                    lookupObject<volScalarField>(inField1),
                    lookupObject<volScalarField>(inField2)
                )
            ).ptr(),
            dict
        );
    }
    else if (inField1Type == ftvector && inField2Type == ftvector)
    {
        setField
        (
            fieldName,
            (
                min
                (
                    lookupObject<volVectorField>(inField1),
                    lookupObject<volVectorField>(inField2)
                )
            ).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Non-equal input field types for min operator with output  "
            << fieldName << nl
            << "Skipping calculation of max." << endl;
    }
}


void functionObjects::fieldProcess::powOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);
    scalar power = readScalar(dict.lookup("power"));

    if (inFieldType == ftscalar)
    {
        setField
        (
            fieldName,
            (
                Foam::pow
                (
                    mag(lookupObject<volScalarField>(inField)),
                    power
                )
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == ftvector)
    {
        setField
        (
            fieldName,
            (
                Foam::pow
                (
                    mag(lookupObject<volVectorField>(inField)),
                    power
                )
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == fttensor)
    {
        setField
        (
            fieldName,
            (
                Foam::pow
                (
                    mag(lookupObject<volTensorField>(inField)),
                    power
                )
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == ftsymmTensor)
    {
        setField
        (
            fieldName,
            (
                Foam::pow
                (
                    mag(lookupObject<volSymmTensorField>(inField)),
                    power
                )
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == ftsphTensor)
    {
        setField
        (
            fieldName,
            (
                Foam::pow
                (
                    mag(lookupObject<volSphericalTensorField>(inField)),
                    power
                )
            ).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Invalid input field for pow operator with output  "
            << fieldName << nl
            << "Skipping calculation of pow." << endl;
    }
}


void functionObjects::fieldProcess::dotOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField1 = dict.lookup<word>("inField1");
    word inField2 = dict.lookup<word>("inField2");
    opFieldTypes inField1Type = volFieldType(inField1);
    opFieldTypes inField2Type = volFieldType(inField2);

    if (inField1Type == ftscalar || inField2Type == ftscalar)
    {
        const volScalarField* sfPtr;
        opFieldTypes otherFieldType = inField2Type;
        word otherFieldName = inField2;

        if (inField1Type == ftscalar)
        {
            sfPtr = & (lookupObject<volScalarField>(inField1));
        }
        else
        {
            sfPtr = & (lookupObject<volScalarField>(inField2));
            otherFieldType = inField1Type;
            otherFieldName = inField1;
        }

        if (otherFieldType == ftscalar)
        {
            setField
            (
                fieldName,
                (
                    (*sfPtr)
                    * (lookupObject<volScalarField>(otherFieldName))
                ).ptr(),
                dict
            );
        }
        else if (otherFieldType == ftvector)
        {
            setField
            (
                fieldName,
                (
                    (*sfPtr)
                    * (lookupObject<volVectorField>(otherFieldName))
                ).ptr(),
                dict
            );

        }
        else if (otherFieldType == fttensor)
        {
            setField
            (
                fieldName,
                (
                    (*sfPtr)
                    * (lookupObject<volTensorField>(otherFieldName))
                ).ptr(),
                dict
            );

        }
        else if (otherFieldType == ftsymmTensor)
        {
            setField
            (
                fieldName,
                (
                    (*sfPtr)
                    * (lookupObject<volSymmTensorField>(otherFieldName))
                ).ptr(),
                dict
            );

        }
        else if (otherFieldType == ftsphTensor)
        {
            setField
            (
                fieldName,
                (
                    (*sfPtr)
                    * (lookupObject<volSphericalTensorField>
                    (otherFieldName))
                ).ptr(),
                dict
            );

        }

    }
    else if (inField1Type == ftvector && inField2Type == ftvector)
    {
        setField
        (
            fieldName,
            (
                (lookupObject<volVectorField>(inField1))
                & (lookupObject<volVectorField>(inField2))
            ).ptr(),
            dict
        );
    }
    else if (inField1Type == ftvector)
    {
        if (inField2Type == fttensor)
        {
            setField
            (
                fieldName,
                (
                    (lookupObject<volVectorField>(inField1))
                    & (lookupObject<volTensorField>(inField2))
                ).ptr(),
                dict
            );
        }
        else if (inField2Type == ftsymmTensor)
        {
            setField
            (
                fieldName,
                (
                    (lookupObject<volVectorField>(inField1))
                    & (lookupObject<volSymmTensorField>(inField2))
                ).ptr(),
                dict
            );
        }
        else if (inField2Type == ftsphTensor)
        {
            setField
            (
                fieldName,
                (
                    (lookupObject<volVectorField>(inField1))
                    & (lookupObject<volSphericalTensorField>(inField2))
                ).ptr(),
                dict
            );
        }
    }
    else if (inField2Type == ftvector)
    {
        if (inField1Type == fttensor)
        {
            setField
            (
                fieldName,
                (
                    (lookupObject<volTensorField>(inField1))
                    & (lookupObject<volVectorField>(inField2))
                ).ptr(),
                dict
            );
        }
        else if (inField1Type == ftsymmTensor)
        {
            setField
            (
                fieldName,
                (
                    (lookupObject<volSymmTensorField>(inField1))
                    & (lookupObject<volVectorField>(inField2))
                ).ptr(),
                dict
            );
        }
        else if (inField1Type == ftsphTensor)
        {
            setField
            (
                fieldName,
                (
                    (lookupObject<volSphericalTensorField>(inField1))
                    & (lookupObject<volVectorField>(inField2))
                ).ptr(),
                dict
            );
        }

    }
    else
    {
        WarningInFunction
            << "Illegal input field combination for dot operator with output  "
            << fieldName << nl
            << "Allowed operations: 'vector.vector', 'vector.tensor' "
            << "and 'tensor.vector'." << nl
            << "Skipping calculation of dot." << endl;
    }
}


void functionObjects::fieldProcess::crossOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField1 = dict.lookup<word>("inField1");
    word inField2 = dict.lookup<word>("inField2");
    opFieldTypes inField1Type = volFieldType(inField1);
    opFieldTypes inField2Type = volFieldType(inField2);

    //check both infields are vectors
    if (inField1Type == ftvector && inField2Type == ftvector)
    {
        setField
        (
            fieldName,
            (
                (lookupObject<volVectorField>(inField1))
                ^ (lookupObject<volVectorField>(inField2))
            ).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Non-vector input field for cross operator with output  "
            << fieldName << nl
            << "Skipping calculation of cross." << endl;
    }
}


void functionObjects::fieldProcess::wallShear(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");

    setField
    (
        fieldName,
        wallShear().ptr(),
        dict
    );
}


void functionObjects::fieldProcess::yPlus(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");

    word shearName = "tauw";

    if (dict.found("shearField"))
    {
        shearName = dict.lookup<word>("shearField");
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    volScalarField::Boundary d = nearWallDist(mesh).y();

    autoPtr<volVectorField> wallShearStress;

    if (foundObject<volVectorField>(shearName))
    {
        wallShearStress.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "wallShearStress",
                    obr_.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                lookupObject<volVectorField>(shearName)
            )
        );
    }
    else
    {
        wallShearStress.reset(wallShear().ptr());
    }

    volScalarField yPlus
    (
        IOobject
        (
            "dimensionlessY",
            obr_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("dimensionlessY", dimless, 0.0)
    );

    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const auto& model =
            lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);

        const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (typeid(currPatch) == typeid(wallFvPatch))
            {
                // Lookup Cmu corresponding to the turbulence model selected
                const scalar Cmu =
                    model.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);
                yPlus.boundaryFieldRef()[patchi].forceAssign
                (
                    model.yPlus(patchi, Cmu)
                );
            }
        }
    }
    else if (foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const icoTurbModel& model =
            lookupObject<icoTurbModel>(turbulenceModel::propertiesName);
        const fvPatchList& patches = mesh.boundary();

        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (typeid(currPatch) == typeid(wallFvPatch))
            {
                // Lookup Cmu corresponding to the turbulence model selected
                const scalar Cmu =
                    model.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);
                yPlus.boundaryFieldRef()[patchi].forceAssign
                (
                    model.yPlus(patchi, Cmu)
                );
            }
        }
    }

    setField(fieldName, yPlus, dict);
}


void functionObjects::fieldProcess::Pe(const dictionary& dict)
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    word fieldName = dict.lookup<word>("fieldName");

    word phiName = "phi";
    if (dict.found("phi"))
    {
        phiName = dict.lookup<word>("phi");
    }

    if
    (
        foundObject<surfaceScalarField>(phiName)
    )
    {
        const surfaceScalarField& phit =
            lookupObject<surfaceScalarField>(phiName);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (phit.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {
            if
            (
                foundObject<cmpTurbModel>(turbulenceModel::propertiesName)
            )
            {
                const cmpTurbModel& model =
                    lookupObject<cmpTurbModel>
                    (
                        turbulenceModel::propertiesName
                    );

                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            fvc::reconstruct
                            (
                                phit/mesh.surfaceInterpolation::deltaCoeffs()
                            )
                            /model.muEff()*model.rho()
                        )
                    ).ptr(),
                    dict
                );
            }
            else if
            (
                foundObject<icoTurbModel>(turbulenceModel::propertiesName)
            )
            {
                const icoTurbModel& model =
                    lookupObject<icoTurbModel>
                    (
                        turbulenceModel::propertiesName
                    );
                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            fvc::reconstruct
                            (
                                phit/mesh.surfaceInterpolation::deltaCoeffs()
                            )
                            /model.nuEff()
                        )
                    ).ptr(),
                    dict
                );
            }

        }
        else if (phit.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            if
            (
                foundObject<cmpTurbModel>(turbulenceModel::propertiesName)
            )
            {
                const cmpTurbModel& model =
                    lookupObject<cmpTurbModel>
                    (
                        turbulenceModel::propertiesName
                    );

                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            fvc::reconstruct
                            (
                                phit/mesh.surfaceInterpolation::deltaCoeffs()

                            )
                            /model.muEff()
                        )
                    ).ptr(),
                    dict
                );
            }
            else if
            (
                foundObject<icoTurbModel>(turbulenceModel::propertiesName)
            )
            {
                const icoTurbModel& model =
                    lookupObject<icoTurbModel>
                    (
                        turbulenceModel::propertiesName
                    );

                setField
                (
                    fieldName,
                    (
                        mag
                        (
                            fvc::reconstruct
                            (
                                phit/mesh.surfaceInterpolation::deltaCoeffs()
                            )
                            /model.nuEff()/model.rho()
                        )
                    ).ptr(),
                    dict
                );
            }

        }
        else
        {
            WarningInFunction
                << "Incorrect dimensions of phi: " << phit.dimensions()
                << "Skipping calculation of Pe." << endl;
        }
    }
    else
    {
        WarningInFunction
            << "Could not find " << phiName
            << " in the database." << nl
            << "Skipping calculation of Pe." << endl;
    }

}


void functionObjects::fieldProcess::Unw(const dictionary& dict)
{
    const word fieldName = dict.lookup<word>("fieldName");
    const word Uname = dict.lookupOrDefault<word>("U", "U");
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volVectorField Unw
    (
        IOobject
        (
            "nearWallVelocity",
            obr_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("Unw", dimensionSet(0, 1, -1, 0, 0, 0, 0), Zero),
        fixedValueFvPatchVectorField::typeName
    );

    if (foundObject<volVectorField>(Uname))
    {
        const volVectorField& Ut = lookupObject<volVectorField>(Uname);
        Unw.primitiveFieldRef() = Ut.primitiveField();
        setWallsToNearCell(Unw);
        setField(fieldName, Unw, dict);
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of Unw." << endl;
    }
}


void functionObjects::fieldProcess::wallHTC(const dictionary& dict)
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    word fieldName = dict.lookup<word>("fieldName");

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField alphaConv
    (
        IOobject
        (
            "wallheattransfercoefficient",
            obr_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "alphaConv",
            dimensionSet(1, 0, -3, -1, 0, 0, 0),
            0.0
        )
    );

    bool useTRef = dict.found("TRef");
    scalar TRef(0);
    const volScalarField* TPtr(nullptr);
    if (useTRef)
    {
        TRef = dict.lookup<scalar>("TRef");
        TPtr = &lookupObject<volScalarField>("T");
    }

    if (foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const cmpTurbModel& model =
            lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);

        const basicThermo& thermo = model.transport();

        forAll(alphaConv.boundaryField(), i)
        {
            if (!alphaConv.boundaryField()[i].coupled())
            {
                if (useTRef)
                {
                    alphaConv.boundaryFieldRef()[i] =
                        (
                            thermo.Cp()->boundaryField()[i]
                           *model.alphaEff()().boundaryField()[i]
                           *fvc::snGrad(*TPtr)->boundaryField()[i]
                        )
                      / stabilise(TPtr->boundaryField()[i] - TRef, VSMALL);
                }
                else
                {
                    alphaConv.boundaryFieldRef()[i] =
                        (
                            thermo.Cp()->boundaryField()[i]
                           *model.alphaEff()().boundaryField()[i]
                           *mesh.boundary()[i].deltaCoeffs()
                        );
                }
            }
        }

    }
    else if (foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const icoTurbModel& model =
            mesh.lookupObject<icoTurbModel>(turbulenceModel::propertiesName);

        forAll(alphaConv.boundaryField(), i)
        {
            if (!alphaConv.boundaryField()[i].coupled())
            {
                if (useTRef)
                {
                    alphaConv.boundaryFieldRef()[i] =
                        (
                            model.Cp()->boundaryField()[i]
                           *model.rho()->boundaryField()[i]
                           *model.alphaEff()().boundaryField()[i]
                           *fvc::snGrad(*TPtr)->boundaryField()[i]
                        )
                      / stabilise(TPtr->boundaryField()[i] - TRef, VSMALL);
                }
                else
                {
                    alphaConv.boundaryFieldRef()[i] =
                        model.Cp()->boundaryField()[i]
                       *model.rho()->boundaryField()[i]
                       *model.alphaEff()().boundaryField()[i]
                       *mesh.boundary()[i].deltaCoeffs();
                }
            }
        }

    }
    else
    {
        FatalErrorInFunction
            << "No valid source for viscosity for wall heat transfer "
            << "calculations. Aborting."
            << exit(FatalError);

    }

    setField(fieldName, alphaConv, dict);
}


void functionObjects::fieldProcess::nearWallValue(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    switch (inFieldType)
    {
        case ftnone:
        {
            WarningInFunction
                 << "Input field " << inField << " does not exist."
                 << nl << "Skipping calculation of " << fieldName << endl;
        } break;

        case ftscalar:
        {
            setField
            (
                fieldName,
                createNearWallField<scalar>
                (
                    lookupObject<volScalarField>(inField)
                ).ptr(),
                dict
            );
        } break;

        case ftvector:
        {
            setField
            (
                fieldName,
                createNearWallField<vector>
                (
                    lookupObject<volVectorField>(inField)
                ).ptr(),
                dict
            );
        } break;

        case fttensor:
        {
            setField
            (
                fieldName,
                createNearWallField<tensor>
                (
                    lookupObject<volTensorField>(inField)
                ).ptr(),
                dict
            );
        } break;

        case ftsymmTensor:
        {
            setField
            (
                fieldName,
                createNearWallField<symmTensor>
                (
                    lookupObject<volSymmTensorField>(inField)
                ).ptr(),
                dict
            );
        } break;

        case ftsphTensor:
        {
            setField
            (
                fieldName,
                createNearWallField<sphericalTensor>
                (
                    lookupObject<volSphericalTensorField>(inField)
                ).ptr(),
                dict
            );
        }
    }
}


void functionObjects::fieldProcess::ddtOp(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);

    // Check field type and create appropriate output
    if (inFieldType == ftscalar)
    {
        setField
        (
            fieldName,
            (
                fvc::ddt(lookupObject<volScalarField>(inField))
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == ftvector)
    {
        setField
        (
            fieldName,
            (
                fvc::ddt(lookupObject<volVectorField>(inField))
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == fttensor)
    {
        setField
        (
            fieldName,
            (
                fvc::ddt(lookupObject<volTensorField>(inField))
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == ftsymmTensor)
    {
        setField
        (
            fieldName,
            (
                fvc::ddt(lookupObject<volSymmTensorField>(inField))
            ).ptr(),
            dict
        );
    }
    else if (inFieldType == ftsphTensor)
    {
        setField
        (
            fieldName,
            (
                fvc::ddt(lookupObject<volSphericalTensorField>(inField))
            ).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Undefined input field type " << inFieldType
            << " for field " << fieldName << nl
            << "Skipping calculation of ddt." << endl;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
functionObjects::fieldProcess::createPatchOnlyField
(
    const GeometricField<Type, fvPatchField, volMesh>& inField,
    const List<wordRe>& patches
)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    tmp<GeometricField<Type, fvPatchField, volMesh>> gf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "nearWallValueField",
                obr_.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            ("nearWallValueField", inField.dimensions(), pTraits<Type>::zero)
        )
    );

    wordRes patchesToWrite(patches);

    forAll(gf->boundaryField(), i)
    {
        if (patchesToWrite.match(gf->boundaryField()[i].patch().name()))
        {
            gf->boundaryFieldRef()[i] = inField.boundaryField()[i];
        }
        else
        {
            gf->boundaryFieldRef()[i] = pTraits<Type>::zero;
        }
    }

    return gf;
}


void functionObjects::fieldProcess::patchWrite(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    word inField = dict.lookup<word>("inField");
    opFieldTypes inFieldType = volFieldType(inField);
    List<wordRe> patches(dict.lookup("patches"));

    switch(inFieldType)
    {
        case ftnone:
        {
            WarningInFunction
                 << "Input field " << inField << " does not exist."
                 << nl << "Skipping calculation of " << fieldName << endl;
        } break;

        case ftscalar:
        {
            setField
            (
                fieldName,
                createPatchOnlyField<scalar>
                (
                    lookupObject<volScalarField>(inField), patches
                ).ptr(),
                dict
            );
        } break;

        case ftvector:
        {
            setField
            (
                fieldName,
                createPatchOnlyField<vector>
                (
                    lookupObject<volVectorField>(inField), patches
                ).ptr(),
                dict
            );
        } break;

        case fttensor:
        {
            setField
            (
                fieldName,
                createPatchOnlyField<tensor>
                (
                    lookupObject<volTensorField>(inField), patches
                ).ptr(),
                dict
            );
        } break;

        case ftsymmTensor:
        {
            setField
            (
                fieldName,
                createPatchOnlyField<symmTensor>
                (
                    lookupObject<volSymmTensorField>(inField), patches
                ).ptr(),
                dict
            );
        } break;

        case ftsphTensor:
        {
            setField
            (
                fieldName,
                createPatchOnlyField<sphericalTensor>
                (
                    lookupObject<volSphericalTensorField>(inField),
                    patches
                ).ptr(),
                dict
            );
        }
    }
}


void functionObjects::fieldProcess::strainRate(const dictionary& dict)
{
    const word fieldName =
        dict.lookupOrDefault<word>("fieldName", "strainRate");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (foundObject<volVectorField>(Uname))
    {
        const volVectorField& Ut = lookupObject<volVectorField>(Uname);

        setField
        (
            fieldName,
            (sqrt(2.0)*mag(symm(fvc::grad(Ut)))).ptr(),
            dict
        );
    }
    else
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of strainRate." << endl;
    }
}


void functionObjects::fieldProcess::Urel(const dictionary& dict)
{
    const word fieldName = dict.lookupOrDefault<word>("fieldName", "Urel");
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (!foundObject<volVectorField>(Uname))
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of Urel." << endl;
        return;
    }

    const fv::options& fvOptionsObj =
        lookupObject<fv::options>(fv::options::typeName);

    if (fvOptionsObj.optionList::size() || mesh_.moving())
    {
        const volVectorField& U = lookupObject<volVectorField>(Uname);

        volVectorField Urel
        (
            IOobject
            (
                fieldName + "temp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U,
            "calculated"
        );

        if (fvOptionsObj.optionList::size())
        {
            for (int i=0; i < (fvOptionsObj.optionList::size()); ++i)
            {
                if (fvOptionsObj.optionList::operator[](i).MRF())
                {
                    fvOptionsObj.optionList::operator[](i).makeRelative(Urel);
                }
            }
        }

        if (mesh_.moving())
        {
            Urel -= fvc::reconstruct(mesh_.phi());
        }

        setField(fieldName, Urel, dict);
    }
    else
    {
        WarningInFunction
            << "Could not find MRF's"
            << " in the database or mesh is stationary." << endl;
    }
}


void functionObjects::fieldProcess::Urtz(const dictionary& dict)
{
    // Load velocity field name
    const word Uname = dict.lookupOrDefault<word>("U", "U");

    if (!foundObject<volVectorField>(Uname))
    {
        WarningInFunction
            << "Could not find " << Uname
            << " in the database." << nl
            << "Skipping calculation of Urtz." << endl;
        return;
    }

    // Load optional specification of the field name
    const word fieldName =
        dict.lookupOrDefault<word>("fieldName", Uname + word("rtz"));

    // Name of the field to be transformed to cylindrical coordinates
    coordinateFrame* framePtr = nullptr;
    if (dict.found("referenceFrame"))
    {
        framePtr = coordinateFrame::lookupNew(mesh_, dict);
    }
    else
    {
        WarningInFunction
            << "Could not find referenceFrame in the dict." << nl
            << "Skipping calculation of Urtz." << endl;
        return;
    }

    const fv::options& fvOptionsObj =
        lookupObject<fv::options>(fv::options::typeName);

    if (fvOptionsObj.optionList::size() || mesh_.moving())
    {
        const volVectorField& U = lookupObject<volVectorField>(Uname);

        volVectorField Urtz
        (
            IOobject
            (
                fieldName + "temp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U,
            "calculated"
        );

        if (fvOptionsObj.optionList::size())
        {
            for (int i=0; i < (fvOptionsObj.optionList::size()); ++i)
            {
                if (fvOptionsObj.optionList::operator[](i).MRF())
                {
                    fvOptionsObj.optionList::operator[](i).makeRelative(Urtz);
                }
            }
        }

        if (mesh_.moving())
        {
            Urtz -= fvc::reconstruct(mesh_.phi());
        }

        // Always does cylindical transformation regardless on original
        // system type of coordinates
        coordSystem::cylindrical csys
        (
            framePtr->coorSys().origin(),
            framePtr->axis()
        );

        // Internal field
        Urtz.primitiveFieldRef() =
            csys.invTransform(mesh_.cellCentres(), Urtz.primitiveField());

        forAll(Urtz.boundaryField(), patchi)
        {
            vectorField& vfp = Urtz.boundaryFieldRef()[patchi];
            vfp = csys.invTransform(mesh_.boundary()[patchi].Cf(), vfp);
        }

        Urtz.correctBoundaryConditions();

        setField(fieldName, Urtz, dict);
    }
    else
    {
        WarningInFunction
            << "Could not find MRF's"
            << " in the database or mesh is stationary." << endl;
    }
}


void functionObjects::fieldProcess::wallHeatFlux(const dictionary& dict)
{
    // Calculate/decompose total heat flux into conduction,
    // convection and radiation components
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    word phaseName = dict.lookupOrDefault("phaseName", word::null);
    const word thermoName =
            IOobject::groupName(basicThermo::dictName, phaseName);

    word convectiveFluxName =
        dict.lookupOrDefault<word>("convectiveFluxName", "convectiveHeatFlux");
    word conductiveFluxName =
        dict.lookupOrDefault<word>("conductiveFluxName", "conductiveHeatFlux");
    word radiativeFluxName =
        dict.lookupOrDefault<word>("radiativeFluxName", "radiativeHeatFlux");
    word totalFluxName =
        dict.lookupOrDefault<word>("totalFluxName","totalHeatFlux");

    volScalarField convectiveHeatFlux
    (
        IOobject
        (
            convectiveFluxName + "temp",
            obr().time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, 0, -3, 0, 0, 0, 0), 0.0)
    );

    volScalarField conductiveHeatFlux
    (
        IOobject
        (
            conductiveFluxName + "temp",
            obr().time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        ("zero", dimensionSet(1, 0, -3, 0, 0, 0, 0), 0.0)
    );

    volScalarField radiativeHeatFlux
    (
        IOobject
        (
            radiativeFluxName + "temp",
            obr().time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, 0, -3, 0, 0, 0, 0), 0.0)
    );

    if (obr().foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const volScalarField& T = obr().lookupObject<volScalarField>("T");
        const cmpTurbModel& model =
            obr().lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);

        const basicThermo& thermo = model.transport();

        forAll(convectiveHeatFlux.boundaryField(), i)
        {
            if (!convectiveHeatFlux.boundaryField()[i].coupled())
            {
                convectiveHeatFlux.boundaryFieldRef()[i] =
                    thermo.Cp()->boundaryField()[i]
                   *model.alphaEff()().boundaryField()[i]
                   *fvc::snGrad(T)->boundaryField()[i];
            }
        }
    }
    else if (obr().foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const volScalarField& T = obr().lookupObject<volScalarField>("T");
        const icoTurbModel& model =
            mesh_.lookupObject<icoTurbModel>(turbulenceModel::propertiesName);

        forAll(convectiveHeatFlux.boundaryField(), i)
        {
            if (!convectiveHeatFlux.boundaryField()[i].coupled())
            {
                convectiveHeatFlux.boundaryFieldRef()[i] =
                    model.Cp()->boundaryField()[i]
                   *model.rho()->boundaryField()[i]
                   *model.alphaEff()->boundaryField()[i]
                   *fvc::snGrad(T)->boundaryField()[i];
            }
        }
    }
    else if (obr().foundObject<solidThermo>(thermoName))
    {
        const volScalarField& T = obr().lookupObject<volScalarField>("T");

        const solidThermo& thermo =
            obr().lookupObject<solidThermo>(thermoName);

        autoPtr<volSymmTensorField> aniAlpha(nullptr);
        if (!thermo.isotropic())
        {
            aniAlpha.reset
            (
                new volSymmTensorField
                (
                    obr().lookupObject<volSymmTensorField>
                    (
                        IOobject::groupName("Anialpha", phaseName)
                    )
                )
            );
        }

        forAll(conductiveHeatFlux.boundaryField(), i)
        {
            if (!conductiveHeatFlux.boundaryField()[i].coupled())
            {
                conductiveHeatFlux.boundaryFieldRef()[i] =
                    thermo.Cp()->boundaryField()[i];

                if (thermo.isotropic())
                {
                    conductiveHeatFlux.boundaryFieldRef()[i] *=
                        thermo.alpha().boundaryField()[i]
                       *fvc::snGrad(T)->boundaryField()[i];
                }
                else
                {
                    conductiveHeatFlux.boundaryFieldRef()[i] *=
                        (
                            conductiveHeatFlux.boundaryField()[i].patch().nf()
                          & aniAlpha->boundaryField()[i]
                        )
                      & fvc::grad(T)->boundaryField()[i];
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Failed to locate fluid turbulence model or solid thermo "
            << "for wall heat transfer calculations. Aborting."
            << exit(FatalError);
    }

    // Radiative heat flux if available
    if (obr().foundObject<volScalarField>("Qr"))
    {
        const volScalarField& Qr = obr().lookupObject<volScalarField>("Qr");
        forAll(radiativeHeatFlux.boundaryField(), i)
        {
            radiativeHeatFlux.boundaryFieldRef()[i] = -Qr.boundaryField()[i];
        }
    }

    setField(convectiveFluxName, convectiveHeatFlux, dict);
    setField(conductiveFluxName, conductiveHeatFlux, dict);
    setField(radiativeFluxName, radiativeHeatFlux, dict);
    setField
    (
        totalFluxName,
        (convectiveHeatFlux + conductiveHeatFlux + radiativeHeatFlux).ptr(),
        dict
    );
}


void functionObjects::fieldProcess::Mach(const dictionary& dict)
{
    // fieldProcess calculating Mach number distributions
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    word fieldName (dict.lookupOrDefault<word>("fieldName","Mach"));

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const volVectorField& U = lookupObject<volVectorField>("U");

    if (foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
    {
        const volScalarField& T = lookupObject<volScalarField>("T");

        const cmpTurbModel& model =
            lookupObject<cmpTurbModel>(turbulenceModel::propertiesName);
        const basicThermo& thermo = model.transport();

        setField
        (
            fieldName,
            (
                mag(U)
              / sqrt((thermo.Cp()/thermo.Cv())*(thermo.Cp()-thermo.Cv())*T)
            ).ptr(),
            dict
        );
    }
    else if (foundObject<icoTurbModel>(turbulenceModel::propertiesName))
    {
        const icoTurbModel& model =
            mesh.lookupObject<icoTurbModel>(turbulenceModel::propertiesName);

        dimensionedScalar bulkMod
        (
            "bulkModulus",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            dict.lookup("bulkModulus")
        );

        setField(fieldName, (mag(U)/sqrt(bulkMod/model.rho())).ptr(), dict);

        Info<< "Warning - the Mach number calculation is "
            << "based on the bulk modulus specified as an input." << nl
            << "This will affect results in case of multiphase flows." << endl;

    }
    else
    {
        FatalErrorInFunction
            << "Cannot compute Mach number for this case. "
            << "Aborting."
            << exit(FatalError);
    }

}


void functionObjects::fieldProcess::skinFriction(const dictionary& dict)
{
    // Name for the field
    word fieldName = dict.lookup<word>("fieldName");
    word shearName = dict.lookupOrDefault<word>("shearField", "tauw");

    // Get Uref and rhoRef from dictionary
    scalar Uref = readScalar(dict.lookup("Uref"));
    dimensionedScalar U2("U2", dimVelocity*dimVelocity, 0.5*magSqr(Uref));
    scalar rho = readScalar(dict.lookup("rhoRef"));
    dimensionedScalar rhoRef("rhoRef", dimDensity, rho);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Get wall shear stress
    autoPtr<volVectorField> wallShearStress;
    if (foundObject<volVectorField>(shearName))
    {
        wallShearStress.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "wallShearStress",
                    obr_.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lookupObject<volVectorField>(shearName)
            )
        );
    }
    else
    {
        wallShearStress.reset(wallShear().ptr());
    }

    volScalarField Cf
    (
        IOobject
        (
            "skinFrictionCoeff",
            obr_.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Cf", dimless, 0.0)
    );

    // init with magnitude of the wallShearStress
    forAll(mesh.boundary(),patchi)
    {
        if  (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            Cf.boundaryFieldRef()[patchi].forceAssign
            (
                mag(wallShearStress->boundaryField()[patchi])
            );
        }
    }

    setField(fieldName, (Cf/(rhoRef*U2)).ptr(), dict);
}


void functionObjects::fieldProcess::cylindricalTransform
(
    const dictionary& dict
)
{
    // Name of the field to be transformed to cylindrical coordinates
    const word inFieldName = dict.lookup<word>("inField");
    const word fieldName =
        dict.lookupOrDefault<word>("fieldName", inFieldName + word("_rtz"));

    if (foundObject<volVectorField>(inFieldName))
    {
        coordSystem::cylindrical csys
        (
            dict.lookupOrDefault<vector>("origin", Zero),
            dict.lookup<vector>("axis")
        );

        // Copy input field
        tmp<volVectorField> outFieldPtr
        (
            new volVectorField
            (
                word("tmp"),
                lookupObject<volVectorField>(inFieldName)
            )
        );

        // Internal field
        outFieldPtr->primitiveFieldRef() =
            csys.invTransform
            (
                mesh_.cellCentres(),
                outFieldPtr->primitiveFieldRef()
            );

        forAll(outFieldPtr->boundaryField(), patchi)
        {
            vectorField& vfp = outFieldPtr->boundaryFieldRef()[patchi];
            vfp = csys.invTransform(mesh_.boundary()[patchi].Cf(), vfp);
        }

        outFieldPtr->correctBoundaryConditions();

        setField(fieldName, outFieldPtr.ptr(), dict);
    }
}


void functionObjects::fieldProcess::proudmanSources(const dictionary& dict)
{
    word fieldName = dict.lookup<word>("fieldName");
    scalar rhoInf_ = dict.lookupOrDefault<scalar>("rhoInf", 1.0);
    scalar alpha0_ = dict.lookupOrDefault<scalar>("alpha0", 343);
    scalar alphae_ = dict.lookupOrDefault<scalar>("alphae", 0.1);
    scalar refAcPower_ = dict.lookupOrDefault<scalar>("refAcPower", 1e-12);

    tmp<volScalarField> acPowerdB
    (
        new volScalarField
        (
            IOobject
            (
                "acPowerdB",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    const tmp<volScalarField> kTurb = getK();
    const tmp<volScalarField> eTurb = getEpsilon();

    tmp<volScalarField> Mt
    (
        sqrt(2*kTurb())
       /dimensionedScalar("vel", dimLength/dimTime, alpha0_)
    );

    tmp<volScalarField> acPower
    (
        alphae_
       *eTurb()
       *pow(Mt, 5)
       *dimensionedScalar( "rho", dimMass/sqr(dimLength)/dimLength, rhoInf_)
    );

    acPowerdB.ref() =
        10*log10
        (
            acPower()
           /dimensionedScalar
            (
                "pv",
                dimMass/dimLength/sqr(dimTime)/dimTime,
                refAcPower_
            )
        );

    setField(fieldName, acPowerdB.ptr(), dict);
}

}
