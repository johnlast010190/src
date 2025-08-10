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

\*---------------------------------------------------------------------------*/

#include "regionState.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "fvMesh/fvMeshStitchers/stationary/fvMeshStitchersStationary.H"
#include "nonConformal/nonConformalFuncs/nonConformalFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace stateFunctions
{
    defineTypeNameAndDebug(regionState, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void regionState::clearNonConformalSetup(const bool initFields)
{
    // Check for the existence of the non-conformal finite-volume addressing
    IOobject polyFacesBfIO
    (
        "polyFaces",
        mesh().pointsInstance(),
        fvMesh::typeName,
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!initFields)
    {
        // Return early if the mesh is conformal or if there are no initialised
        // fields in the time directory
        if
        (
            !fileHandler().isFile(polyFacesBfIO.objectPath())
        || !isDir(mesh().time().timePath())
        ) return;
    }
    else
    {
        // Return early for conformal meshes. Needed by the FOAM GUI.
        if (!fileHandler().isFile(polyFacesBfIO.objectPath())) return;
    }

    // Remove the directory containing the finite-volume addressing field
    rmDir(polyFacesBfIO.path());

    // Delete the nonConformalCouplesDict dictionary, if it wasn't written
    // by the current caseSetup run.
    if (!initFields)
    {
        IOobject nccDict
        (
            "nonConformalCouplesDict",
            mesh().time().caseSystem(),
            mesh().dbDir(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if
        (
            fileHandler().isFile(nccDict.objectPath())
        && !system().found("nonConformalCouplesDict")
        )
        {
            rm(nccDict.objectPath());
        }
    }

    // Remove the non-conformal patches from the boundary mesh
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    List<polyPatch*> newPolyPatches;

    // Clone all patches, excluding the non-conformal ones.
    label newIndexi = 0;
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (isA<nonConformalPolyPatch>(pp)) continue;

        const directPolyPatch& dpp =
            refCast<const directPolyPatch>(patches[patchi]);

        newPolyPatches.append
        (
            dpp.clone(patches, newIndexi++, dpp.size(), dpp.start()).ptr()
        );
    }

    // Re-patch the mesh
    meshPtr_->removeFvBoundary();
    meshPtr_->addFvPatches(newPolyPatches);

    // Exclude CBC patch field entries for removed non-conformal patches
    const dictionary& fieldDef(fieldDefinitions());
    forAllConstIter(dictionary, fieldDef, iter)
    {
        const word fieldName(iter().keyword());

        autoPtr<IOdictionary> cbcPtr
        (
            new IOdictionary
            (
                IOobject
                (
                    fieldName + ".reference",
                    mesh().time().timeName(),
                    "uniform/boundaryConditions",
                    localDb().registry(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
        IOdictionary& cbcDictRef
        (
            localDb().registry().template lookupObjectRef<IOdictionary>
            (
                fieldName + ".reference"
            )
        );

        forAllIter(dictionary, cbcDictRef, iter)
        {
            const entry& e = iter();
            const word pfType(e.dict().lookup("type"));

            if
            (
                pfType == nonConformalCyclicFvPatch::typeName
             || pfType == nonConformalErrorFvPatch::typeName
             || pfType == nonConformalProcessorCyclicFvPatch::typeName
            )
            {
                cbcDictRef.remove(e.keyword());
            }
        }

        cbcDictRef.regIOobject::write();
    }
}


Xfer<dictionary> regionState::initialPatchDict(const label& patchI) const
{
    const polyBoundaryMesh& polyPatches = meshPtr_->boundaryMesh();

    dictionary tpd;

    //type entry
    tpd.add
    (
        "type",
        polyPatches[patchI].type(),
        true
    );

    //size entry
    tpd.add
    (
        "nFaces",
        polyPatches[patchI].size(),
        true
    );
    //start entry
    tpd.add
    (
        "startFace",
        polyPatches[patchI].start(),
        true
    );

    return tpd.xfer();
}

word regionState::correctPatchPhysicalType(dictionary& patchDict)
{
    word bType = word(patchDict.lookup("type"));
    word pType = patchDict.lookupOrDefault<word>("physicalType", "");

    //translate to or add physical type to patch
    if (bType == "patch" && pType == "")
    {
        pType = "opening";
    }

    if (bType == "inlet" || bType == "outlet")
    {
        pType = bType;
        bType = "patch";
    }

    //physicalType checking
    if (pType != "" && bType == "patch")
    {
        if
        (
            pType == "inlet"
            || pType == "outlet"
            || pType == "opening"
            || pType == "symmetry"
        )
        {
            patchDict.add<word>
            (
                "physicalType",
                pType,
                true
            );
            patchDict.add<word>
            (
                "type",
                bType,
                true
            );
        }
        else
        {
            FatalError << pType
               << " is a currently unsupported"
               << " physicalType for patch "
               << "boundaries. Currently "
               << "supported physicalTypes are: "
               << "inlet, outlet, opening and symmetry."
               << exit(FatalError);

        }
    }

    return bType;
}

bool regionState::validMultiPatchModify
(
    const dictionary& defaults,
    const dictionary& patchDict
)
{

    if (patchDict.found("name"))
    {
        FatalErrorInFunction
            << "'name' entry in multi-patch modification entry: "
            << patchDict
            << exit(FatalError);

        return false;
    }

    //only allow non-mapping, non-coupled types
    word bType(patchDict.lookup("type"));

    wordList allowedMultiSetPatchTypes
    (
        defaults.subDict("mesh").lookup("allowedMultiSetPatchTypes")
    );

    bool typeAllowed = false;

    forAll(allowedMultiSetPatchTypes, i)
    {
        if (bType == allowedMultiSetPatchTypes[i])
        {
            typeAllowed = true;
            break;
        }
    }

    if (!typeAllowed)
    {
        FatalErrorInFunction
            << "Inelligible type entry in multi-patch modification entry: "
            << bType << nl
            << "Allowed patch type entries are: " << allowedMultiSetPatchTypes
            << exit(FatalError);

        return false;
    }

    return true;
}

void regionState::modifyBoundaryPatch
(
    const fvMesh& mesh,
    const dictionary& dict,
    label patchI
)
{

    polyBoundaryMesh& polyPatches
       = const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    Info<< "Modifying boundary patch: "
         << polyPatches[patchI].name()
         << endl;

    fvBoundaryMesh& fvPatches
        = const_cast<fvBoundaryMesh&>(mesh.boundary());

    //renaming
    dictionary patchDict(dict);
    word patchName = patchDict.lookupOrDefault<word>
    (
        "name",
        polyPatches[patchI].name()
    );
    patchDict.remove("name");


    polyPatches.set
    (
        patchI,
        polyPatch::New
        (
            patchName,
            patchDict,
            patchI,
            polyPatches
        )
    );

    fvPatches.set
    (
        patchI,
        fvPatch::New
        (
            polyPatches[patchI],
            // needs const reference
            mesh.boundary()
        )
    );
}


void regionState::modifyBoundaryMesh
(
    const dictionary& defaults,
    List<bool>& modifiedPatch
)
{
    const dictionary& patchDefaults
    (
        defaults.subDict("mesh").subDict("patchDefaults")
    );

    // No modification if 'boundaryMesh' entry is not present
    if (!input().found("boundaryMesh")) return;

    // Clear any residual non-conformal setup
    clearNonConformalSetup();

    const dictionary& bmDict(input().subDict("boundaryMesh"));

    //fail if any entry fails to find at least one hit
    const bool strictCheck(switches().strictPatchNameChecking());

    //ref to boundary mesh
    const fvBoundaryMesh& bmesh(meshPtr_->boundary());

    //We can overwrite the patch multiple times, make sure it sticks in the
    //following order of importance
    //1. Explicit
    //2. Group
    //3. RegExp
    //Where the same interface gets multiple hits, the last
    //entry will win.

    //Construct dictionary to record hits for each entry
    dictionary hitDict;

    // 1. Handle explicit patch names.
    forAllConstIter(dictionary, bmDict, iter)
    {
        if (iter().isDict() && !iter().keyword().isPattern())
        {
            word patchName = iter().keyword();
            label patchi = bmesh.findPatchID(patchName);
            const dictionary& indict(iter().dict());

            //check if renaming has already occured
            label newNamePatchi = -1;

            if (indict.found("name"))
            {
                newNamePatchi = bmesh.findPatchID
                (
                    word(indict.lookup("name"))
                );

                if (patchi == -1)
                {
                    patchi = newNamePatchi;
                }
            }

            if (patchi != -1)
            {

                dictionary patchDict(initialPatchDict(patchi));
                patchDict.merge(indict);

                //remove old style inlet/outlet patches
                word boundaryType = correctPatchPhysicalType(patchDict);

                //add boiler plate from coupling patches
                if (patchDefaults.found(boundaryType))
                {
                    dictionary mergedDict
                    (
                        patchDefaults.subDict(boundaryType)
                    );
                    mergedDict.merge(patchDict); // overwrites existing
                    patchDict.transfer(mergedDict);
                }

                modifyBoundaryPatch
                (
                    mesh(),
                    patchDict,
                    patchi
                );

                modifiedPatch[patchi] = true;

                if (!hitDict.found(patchName))
                {
                    hitDict.add(patchName, "");
                }
            }
            /*
            else if (strictCheck)
            {
                FatalErrorInFunction
                    << "No patch matching boundary "
                    << "name: " << iter().keyword()
                    << ". Disable strictPatchNameChecking to prevent"
                    << " this exception (1)."
                    << exit(FatalError);
            }*/
        }
    }

    // 2. Patch-groups.
    if (bmDict.size())
    {
        for
        (
            IDLList<entry>::const_reverse_iterator iter = bmDict.rbegin();
            iter != bmDict.rend();
            ++iter
        )
        {
            const entry& e = iter();

            if (e.isDict() && !e.keyword().isPattern())
            {
                const labelList patchIDs = bmesh.findIndices
                (
                    e.keyword(),
                    true                    // use patchGroups
                );

                /*
                if (!patchIDs.size() && strictCheck)
                {
                    FatalErrorInFunction
                        << "No patch matching boundary "
                        << "name: " << iter().keyword()
                        << ". Disable strictPatchNameChecking to prevent"
                        << " this exception (2)."
                        << exit(FatalError);
                }
                */

                forAll(patchIDs, i)
                {
                    label patchi = patchIDs[i];

                    if (!modifiedPatch[patchi])
                    {
                        dictionary patchDict(initialPatchDict(patchi));
                        patchDict.merge(e.dict());
                        //remove old style inlet/outlet patches
                        correctPatchPhysicalType(patchDict);

                        if (validMultiPatchModify(defaults, patchDict))
                        {
                            modifyBoundaryPatch
                            (
                                *meshPtr_,
                                patchDict,
                                patchi
                            );

                            modifiedPatch[patchi] = true;

                            if (!hitDict.found(e.keyword()))
                            {
                                hitDict.add(e.keyword(), "");
                            }
                        }

                    }
                }
            }
        }
    }

    //check if all non-regexp entries have been used
    forAllConstIter(dictionary, bmDict, iter)
    {
        if (iter().isDict() && !iter().keyword().isPattern())
        {
            if (!hitDict.found(iter().keyword()) && strictCheck)
            {
                FatalErrorInFunction
                    << "No patch matching boundary "
                    << "name/group: " << iter().keyword()
                    << ". Disable strictPatchNameChecking to prevent"
                    << " this exception."
                    << exit(FatalError);
            }
        }
    }


    // 3. Wildcard patch overrides
    if (bmDict.size())
    {
        forAllReverse(bmesh, patchi)
        {
            if (!modifiedPatch[patchi])
            {
                if
                (
                    bmDict.found(bmesh[patchi].name())
                    && bmesh[patchi].type() != "processor"
                )
                {
                    dictionary patchDict(initialPatchDict(patchi));
                    patchDict.merge(bmDict.subDict(bmesh[patchi].name()));

                    //remove old style inlet/outlet patches
                    correctPatchPhysicalType(patchDict);

                    if (validMultiPatchModify(defaults, patchDict))
                    {
                        modifyBoundaryPatch
                        (
                            *meshPtr_,
                            patchDict,
                            patchi
                        );

                        modifiedPatch[patchi] = true;
                    }
                }
            }
        }
    }

    forAll(modifiedPatch, i)
    {
        if (modifiedPatch[i])
        {
            fvMesh& pm(const_cast<fvMesh&>(*meshPtr_));
            pm.recalculatePatches();

            //find pointMesh and delete if present
            if (meshPtr_->foundObject<pointMesh>(pointMesh::typeName))
            {
                autoPtr<pointMesh> pntMesh
                (
                    &const_cast<pointMesh&>
                    (
                        meshPtr_->lookupObject<pointMesh>
                        (pointMesh::typeName)
                    )
                );

                pntMesh->release();
            }

            break;
        }
    }
}

template<>
void regionState::storeInternalField
(
    const word& name, const Field<scalar>& v
)
{
    scalarFields_.append
    (
        new Tuple2<word, Field<scalar>>
        (
            name,
            v
        )
    );
}

template<>
void regionState::storeInternalField(const word& name, const Field<vector>& v)
{
    vectorFields_.append
    (
        new Tuple2<word, Field<vector>>
        (
            name,
            v
        )
    );
}
template<>
void regionState::storeInternalField(const word& name, const Field<tensor>& v)
{
    tensorFields_.append
    (
        new Tuple2<word, Field<tensor>>
        (
            name,
            v
        )
    );
}
template<>
void regionState::storeInternalField
(
    const word& name,
    const Field<symmTensor>& v
)
{
    symmTensorFields_.append
    (
        new Tuple2<word, Field<symmTensor>>
        (
            name,
            v
        )
    );
}
template<>
void regionState::storeInternalField
(
    const word& name,
    const Field<sphericalTensor>& v
)
{
    sphericalTensorFields_.append
    (
        new Tuple2<word, Field<sphericalTensor>>
        (
            name,
            v
        )
    );
}

template<>
bool regionState::setInternalValue<scalar>(const word& name)
{
    typedef GeometricField<scalar, fvPatchField, volMesh> GeoField;

    if (localDb().registry().foundObject<GeoField>(name))
    {
        GeoField& f
            = const_cast<GeoField&>(localDb().registry().lookupObject<GeoField>(name));

        forAll(scalarFields_, i)
        {
            if (scalarFields_[i].first() == name)
            {
                Info<< "   Initialising " << f.name()
                     << " with existing values" << endl;

                f.primitiveFieldRef() = scalarFields_[i].second();
                return true;
            }
        }
    }

    return false;
}

template<>
bool regionState::setInternalValue<vector>(const word& name)
{
    typedef GeometricField<vector, fvPatchField, volMesh> GeoField;

    if (localDb().registry().foundObject<GeoField>(name))
    {
        GeoField& f
            = const_cast<GeoField&>(localDb().registry().lookupObject<GeoField>(name));

        forAll(vectorFields_, i)
        {
            if (vectorFields_[i].first() == name)
            {
                Info<< "   Initialising " << f.name()
                     << " with existing values" << endl;

                f.primitiveFieldRef() = vectorFields_[i].second();
                return true;
            }
        }
    }

    return false;
}

template<>
bool regionState::setInternalValue<tensor>(const word& name)
{
    typedef GeometricField<tensor, fvPatchField, volMesh> GeoField;

    if (localDb().registry().foundObject<GeoField>(name))
    {
        GeoField& f
            = const_cast<GeoField&>(localDb().registry().lookupObject<GeoField>(name));

        forAll(tensorFields_, i)
        {
            if (tensorFields_[i].first() == name)
            {
                Info<< "   Initialising " << f.name()
                     << " with existing values" << endl;

                f.primitiveFieldRef() = tensorFields_[i].second();
                return true;
            }
        }
    }

    return false;
}

template<>
bool regionState::setInternalValue<symmTensor>(const word& name)
{
    typedef GeometricField<symmTensor, fvPatchField, volMesh> GeoField;

    if (localDb().registry().foundObject<GeoField>(name))
    {
        GeoField& f
            = const_cast<GeoField&>(localDb().registry().lookupObject<GeoField>(name));

        forAll(symmTensorFields_, i)
        {
            if (symmTensorFields_[i].first() == name)
            {
                Info<< "   Initialising " << f.name()
                     << " with existing values" << endl;

                f.primitiveFieldRef() = symmTensorFields_[i].second();
                return true;
            }
        }
    }

    return false;
}

template<>
bool regionState::setInternalValue<sphericalTensor>(const word& name)
{
    typedef GeometricField<sphericalTensor, fvPatchField, volMesh> GeoField;

    if (localDb().registry().foundObject<GeoField>(name))
    {
        GeoField& f
            = const_cast<GeoField&>(localDb().registry().lookupObject<GeoField>(name));

        forAll(sphericalTensorFields_, i)
        {
            if (sphericalTensorFields_[i].first() == name)
            {
                Info<< "   Initialising " << f.name()
                     << " with existing values" << endl;

                f.primitiveFieldRef() = sphericalTensorFields_[i].second();
                return true;
            }
        }
    }

    return false;
}



void regionState::readInitialFields()
{
    //read fields and store components
    const wordList& fnList(this->fieldNames());
    const dictionary& fieldDefDicts(this->fieldDefinitions());

    forAll(fnList, i)
    {
        word field(fnList[i]);

        const dictionary& fd
        (
            fieldDefDicts.subDict(field).subDict("fieldDefinition")
        );

        fieldType ft
        (
            stateFunction::fieldTypeNames_[word(fd.lookup("type"))]
        );


        switch (ft)
        {

            case ftScalar :
            {
                readInitialTypeField<scalar>(fd, field);
            } break;

            case ftVector :
            {
                readInitialTypeField<vector>(fd, field);
            } break;

            case ftTensor :
            {
                readInitialTypeField<tensor>(fd, field);
            } break;

            case ftSymmTensor :
            {
                readInitialTypeField<symmTensor>(fd, field);
            } break;

            case ftSphericalTensor :
            {
                readInitialTypeField<sphericalTensor>(fd, field);
            } break;

            default:
            {
                FatalErrorInFunction
                    << "Invalid field type: '"
                    << word(fd.lookup("type"))
                    << "'" << exit(FatalError);
            }
        }
    }

}

void regionState::instantiateFields()
{
    const dictionary& fieldDef(fieldDefinitions());

    //loop through fields dict, create fields and store in database
    forAllConstIter(dictionary, fieldDef, iter)
    {
        dictionary fd(iter().dict().subDict("fieldDefinition"));
        word fieldName(iter().keyword());

        //if boundary is not reset, then reuse old boundaries
        if (!switches().resetBoundaryFields())
        {
            if (initialFieldBoundaryDictsPtr_->found(fieldName))
            {
                fd.subDict("boundaryConditions").merge
                (
                    initialFieldBoundaryDictsPtr_->subDict(fieldName)
                );
            }
        }
        //if boundary is reset, but internal field is not
        //try to use existing boundaries for coupled patches
        else if (!switches().resetInternalFields())
        {

            if (initialFieldBoundaryDictsPtr_->found(fieldName))
            {

                const dictionary& oldFieldBCs
                (
                    initialFieldBoundaryDictsPtr_->subDict(fieldName)
                );

                dictionary& newFieldBCs
                (
                    fd.subDict("boundaryConditions")
                );

                forAll(meshPtr_->boundary(), pi)
                {
                    if (meshPtr_->boundary()[pi].coupled())
                    {
                        const fvPatch& cPatch(meshPtr_->boundary()[pi]);

                        if
                        (
                           !newFieldBCs.found(cPatch.name())
                         && oldFieldBCs.found(cPatch.name())
                        )
                        {
                            newFieldBCs.add
                            (
                                cPatch.name(),
                                oldFieldBCs.subDict(cPatch.name())
                            );
                        }
                    }
                }
            }
        }

        fieldType ft
        (
            stateFunction::fieldTypeNames_[fd.lookup("type")]
        );

        switch (ft)
        {

            case ftScalar :
            {
                createTypeField<scalar>(fd, fieldName);
            } break;

            case ftVector :
            {
                createTypeField<vector>(fd, fieldName);
            } break;

            case ftTensor :
            {
                createTypeField<tensor>(fd, fieldName);
            } break;

            case ftSymmTensor :
            {
                createTypeField<symmTensor>(fd, fieldName);
            } break;

            case ftSphericalTensor :
            {
                createTypeField<sphericalTensor>(fd, fieldName);
            } break;

            default:
            {
                FatalErrorInFunction
                    << "Invalid field type: '"
                    << word(fd.lookup("type"))
                    << "'" << exit(FatalError);
            }
        }
    }
}


void regionState::updateInitialisation(fieldInit& fIS)
{
    //logic

    //if resetInternalFields == true
    //then intialialise
    //else
    //  if stored internalfield exists - revert
    //  else initialise

    if (switches().resetInternalFields())
    {
        fIS.correct();
        // Write and re-construct boundary with updated internal field
        reInitBoundary(fIS);
    }
    else
    {
        //search for stored internal field

        if (setInternalValue<scalar>(fIS.name()));
        else if (setInternalValue<vector>(fIS.name()));
        else if (setInternalValue<tensor>(fIS.name()));
        else if (setInternalValue<symmTensor>(fIS.name()));
        else if (setInternalValue<sphericalTensor>(fIS.name()));
        else
        {
            Info<< "   No existing field found for " << fIS.name()
                 << ". Re-initialising" << endl;

            fIS.correct();
            // Write and re-construct boundary with updated internal field
            reInitBoundary(fIS);
        }
    }
}

void regionState::reInitBoundary(const fieldInit& fIS)
{
    const dictionary& fd = fIS.fieldDict();
    meshType met
    (
        stateFunction::meshTypeNames_
        [
            fd.lookupOrDefault<word>("meshType", "volume")
        ]
    );

    switch (met)
    {
        case meVolume:
        {
            reInitBoundaryTypes<fvPatchField, volMesh>(fIS);
        } break;
        case meSurface:
        {
            reInitBoundaryTypes<fvsPatchField, surfaceMesh>(fIS);
        } break;
        case mePoint:
        {
            reInitBoundaryTypes<pointPatchField, pointMesh>(fIS);
        } break;
    }
}

void regionState::initialiseFields()
{

    const dictionary& fieldDicts(fieldDefinitions());

    // create field initialisation list
    PtrList<fieldInit> fIS(0);

    // append new initializations
    forAllConstIter(dictionary, fieldDicts, iter)
    {
        fIS.append
        (
            fieldInit::New
            (
                *meshPtr_,
                localDb(),
                iter().dict(),
                iter().keyword()
            )
        );
    }


    //use priority data to schedule initialisation
    //velocity is priority 2, temperature is priority 4, everything else is 10
    for (label pr = 1; pr <= 100; pr++)
    {
        forAll(fIS, i)
        {
            if (fIS[i].initDict().lookupOrDefault<scalar>("priority", 10) == pr)
            {
                updateInitialisation(fIS[i]);
            }
        }
    }
}

void regionState::clearFieldFiles()
{
    word timeInstance = mesh().time().timeName();

    if
    (
        timeInstance != "constant"
     && isDir(mesh().objectRegistry::path(timeInstance))
    )
    {
        fileName regionDir = mesh().objectRegistry::path(timeInstance);

        if (regionName() != polyMesh::defaultRegion)
        {
            regionDir = fileName(regionDir/regionName());
        }

        fileNameList ObjectNames = readDir(regionDir);

        forAll(ObjectNames, objI)
        {
            if (!fieldDefinitions().found(ObjectNames[objI]))
            {
                rm(regionDir/ObjectNames[objI]);
            }
        }
    }
}

void regionState::clearInstanceFields()
{
    DynamicList<word> namesToCheckout;
    forAll(localDb().registry().names(), nameI)
    {
        if (localDb().registry().foundObject<regIOobject>(localDb().registry().names()[nameI]))
        {
            word name = localDb().registry().names()[nameI];
            namesToCheckout.append(name);
        }
    }
    for (auto const& i : namesToCheckout)
    {
        Info<<"Checkout field "<< i <<endl;
        localDb().registry().checkOut(localDb().registry().lookupObjectRef<regIOobject>(i));
    }
    scalarFields_.clear();
    vectorFields_.clear();
    tensorFields_.clear();
    symmTensorFields_.clear();
    symmTensorFields_.clear();
    solutionRegPtr_.release();
}

void regionState::sync() const
{
    Pstream::barrier();
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionState::regionState
(
    const word regionName,
    const dictionary& regionInput,
    const dictionary& defaults,
    const stateFunction& global,
    const stateIndex& index,
    word meshName
)
:
    stateFunction
    (
        regionName,
        regionInput,
        defaults,
        global,
        index,
        meshName
    ),
    meshPtr_(),
    scalarFields_(),
    vectorFields_(),
    tensorFields_(),
    symmTensorFields_(),
    sphericalTensorFields_(),
    initialFieldBoundaryDictsPtr_(new dictionary())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

regionState::~regionState()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void regionState::createObjects(const Time& runTime)
{
    //read mesh for field read
    if (runTime.foundObject<fvMesh>(meshName()))
    {
        meshPtr_ = &runTime.lookupObjectRef<fvMesh>(meshName());
    }
    else
    {
        meshPtr_ = new fvMesh
            (
                IOobject
                (
                    meshName(),
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            );
            meshPtr_->objectRegistry::store();
       //runTime.db().store(meshPtr_);
      // store(meshPtr_);
    }

    solutionRegPtr_.set
    (
        new fvSolutionRegistry
        (
            runTime,
            *meshPtr_,
            regionName()
        )
    );

    // Check for the existence of the non-conformal finite-volume addressing.
    // Needed by the FOAM GUI.
    IOobject polyFacesBfIO
    (
        "polyFaces",
        mesh().pointsInstance(),
        fvMesh::typeName,
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!fileHandler().isFile(polyFacesBfIO.objectPath()))
    {
        // Read old fields and store components
        if
        (
            !switches().resetInternalFields()
         || !switches().resetBoundaryFields()
        )
        {
            Info<< nl << "Storing existing fields" << endl;
            Info<<       "=======================" << endl;
            readInitialFields();
            Info<< endl;
        }
    }
}


void regionState::modifyMesh(scalar writePause)
{
    //modify mesh
    if (switches().resetBoundaryMesh())
    {
        List<bool> modifiedPatch(mesh().boundaryMesh().size(), false);
        modifyBoundaryMesh(defaults(), modifiedPatch);

        if (distributed_ || collated_)
        {
            if (!mesh().boundaryMesh().write())
            {
                FatalError << "Boundary mesh write failure. "
                           << exit(FatalError);
           }
        }
        else
        {

            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {

                sync();
                if (Pstream::myProcNo() == procI)
                {
                    if (!meshPtr_->boundaryMesh().write())
                    {
                        FatalError << "Boundary mesh write failure. "
                                   << exit(FatalError);
                    }
                }
                Foam::nanosleep(writePause);
            }
            sync();
        }


        //remove stored boundary entries for modified patches
        //will use defaults instead
        forAllIter(dictionary, initialFieldBoundaryDictsPtr_(), iter)
        {
            dictionary& fd(iter().dict());

            forAll(modifiedPatch, i)
            {
                if (modifiedPatch[i])
                {
                    word removePatch(meshPtr_->boundaryMesh()[i].name());

                    fd.remove(removePatch);
                }
            }
        }
    }
}


void regionState::createFields()
{
    Info<< nl << regionName() << endl;
    Info<< "---------------------------------------------------"
         << "-------------------------------------" << endl;

    // Clear eventual residual files from a previous setup with
    // non-conformal meshes. This is a helper workaround function for the
    // FOAM GUI to be able to clean previous NCC setups, but this workflow
    // should be refactored in the future.
    clearNonConformalSetup(true);

    if (mesh().needsStitching())
    {
        // On-demand construction of non-conformal couples
        Info<< nl;
        nonConformalFuncs::createNonConformalCouples(*meshPtr_);
    }
    else if (mesh().stitcher().stitches())
    {
        // Update non-conformal mesh data
        Info<< nl;

        // Clear mesh geometry and addressing
        meshPtr_->conform();

        // Re-stitch the mesh
        fvMeshStitchers::stationary stitcher(*meshPtr_);
        stitcher.connect(false, true, true);

        Info<< nl;
    }

    //create fields
    {
        if (switches().printBoundaries())
        {
            Info<< nl << nl << "Boundary conditions" << endl;
            Info<<             "-------------------" << endl;
        }

        instantiateFields();
    }

    //initialise fields
    {
        Info<< nl << "Initialisation" << endl;
        Info<<       "--------------" << nl << endl;
        initialiseFields();
    }

    //write fields
    {
        localDb().objectRegistry::write();
        localDb().mesh().objectRegistry::write();
    }

    //delete unused fields if "deleteUnusedFields==true"
    if (switches().deleteUnusedFields())
    {
        clearFieldFiles();
    }
    if (switches().cleanInstanceFields())
    {
        clearInstanceFields();
    }

    Info<< "---------------------------------------------------"
         << "-------------------------------------" << endl;
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
