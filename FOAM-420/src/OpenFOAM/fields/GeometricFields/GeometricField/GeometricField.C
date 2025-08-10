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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "db/Time/Time.H"
#include "db/IOobjects/IOdictionary/localIOdictionary.H"
#include "meshes/data/data.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define checkField(gf1, gf2, op)                                               \
if ((gf1).mesh() != (gf2).mesh())                                              \
{                                                                              \
    FatalErrorInFunction                                                       \
        << "different mesh for fields "                                        \
        << (gf1).name() << " and " << (gf2).name()                             \
        << " during operation " <<  op                                         \
        << abort(FatalError);                                                  \
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::
applyCentralizedBoundaryConditions
(
    dictionary& fieldDict
)
{
    //Function overview:
    // Update the on-file boundaryField definition using central boundary
    // conditions.
    // Initially only support replacement mode
    // i.e. entire patch boundary definition is replaced by CBC
    // In future consider adding "merge/modify" mode so micro adjustments are
    // allowed

    // Algorithm:
    // 1. Loop over CBC entries
    // 2. Check that they are valid (have not been applied before)
    // 3. Remove all entries in boundaryField matching CBC entry
    // (suports regular expressions and patchGroups)
    // 4. In a second loop add all valid CBC entries to boundaryField

    //modification indicator for CBC reporting
    bool modified = false;


    //create central boundary condition dictionary
    //read if present
    IOdictionary cbcDict
    (
        IOobject
        (
            this->name(),
            this->time().caseSystem(),
            "boundaryConditions",
            this->db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    // check for common formatting mistake - illegal use of "boundaryField"
    // sub-dictionary
    if (cbcDict.found("boundaryField"))
    {
        FatalErrorInFunction
            << "Illegal keyword `boundaryField` found in the boundary "
            << "conditions dictionary of field " << this->name() << " in file "
             << nl << cbcDict.filePath() << "." << nl
            << abort(FatalError);
    }

    // There is no central BC definition
    // CBC system disabled, no generation of boundary reference
    if (!cbcDict.typeHeaderOk<IOdictionary>())
    {
        return modified;
    }

    //CBC exists

    //Load(if present) and store the central boundary condition reference
    {
        // First remove existing one if present
        IOdictionary* dictPtr =
            this->db().template lookupObjectRefPtr<IOdictionary>
            (
                this->name() + ".reference"
            );
        if (dictPtr)
        {
            dictPtr->checkOut();
        }

        autoPtr<IOdictionary> cbcRefPtr
        (
            //new localIOdictionary
            new IOdictionary
            (
                IOobject
                (
                    //the "reference" postfix is so that this object doesnt
                    //clash with the field of the same name
                    this->name() + ".reference",
                    this->time().timeName(),
                    "uniform/boundaryConditions",
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        cbcRefPtr.ptr()->store();
    }

    //lookup the central boundary reference as non-const reference
    IOdictionary& cbcRefDict
    (
        this->db().template lookupObjectRef<IOdictionary>
        //this->mesh().time().template lookupObjectRef<IOdictionary>
        (
            this->name() + ".reference"
        )
    );

    const typename GeoMesh::BoundaryMesh& bmesh(this->mesh().boundary());

    //Loop through CBC and remove boundaryField entries that match
    //valid CBC entries

    //make sure the field dictionary contains a boundaryField sub-dictionary
    if (!fieldDict.found("boundaryField"))
    {
        fieldDict.add("boundaryField", dictionary());
        modified = true;
    }
    dictionary& boundaryFieldDict(fieldDict.subDict("boundaryField"));

    //record entry validity
    HashTable<bool> activeEntry(cbcDict.size());

    forAllConstIter(dictionary, cbcDict, iter)
    {
        const entry& e = iter();
        word ename(e.keyword());

        activeEntry.set(ename, true);

        if (!e.isDict())
        {
            activeEntry.set(ename, false);
        }
        else // SHA1 check
        {
            SHA1Digest referenceDigest;
            if (cbcRefDict.found(e.keyword()))
            {
                referenceDigest = cbcRefDict.subDict(e.keyword()).digest();
            }

            if (referenceDigest == e.dict().digest())
            {
                activeEntry.set(ename, false);
            }
        }

        //check for matches
        if (activeEntry[ename])
        {
            const labelList patchIDs = bmesh.findIndices
            (
                e.keyword(),
                true                    // use patchGroups
            );

            //delete all patch name matches from the boundaryField
            forAll(patchIDs, i)
            {
                modified =
                (
                    boundaryFieldDict.remove(bmesh[patchIDs[i]].name())
                 || modified
                );
            }
            //delete entry matching current input entry
            //prevents existing ordering from modifying user intent
            modified = (boundaryFieldDict.remove(e.keyword()) || modified);
        }
    }

    //Add valid CBC entries to boundaryField
    forAllConstIter(dictionary, cbcDict, iter)
    {
        const entry& e = iter();
        word ename(e.keyword());

        if (activeEntry[ename])
        {
            //add CBC entry to boundaryField
            boundaryFieldDict.add(e.keyword(), e.dict(), true);

            //update CBC application reference
            // Remove subdict first, to prevent merging
            cbcRefDict.remove(e.keyword());
            cbcRefDict.add(e.keyword(), e.dict(), true);
        }
    }

    return modified;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields
(
    const dictionary& dict,
    const bool cbcSwitch
)
{
    Internal::readField(dict, "internalField");

    if (cbcSwitch)
    {
        if (applyCentralizedBoundaryConditions(const_cast<dictionary&>(dict)))
        {
            Info<< "CBC modified field " << this->name() << nl << endl;
        }
    }

    boundaryField_.readField(*this, dict.subDict("boundaryField"));

    if (dict.found("referenceLevel"))
    {
        Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));

        Field<Type>::operator+=(fieldAverage);

        forAll(boundaryField_, patchi)
        {
            boundaryField_[patchi].forceAssign(boundaryField_[patchi] + fieldAverage);
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields
(
     const bool cbcSwitch
)
{
    const localIOdictionary dict
    (
        IOobject
        (
            this->name(),
            this->instance(),
            this->local(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        typeName
    );

    this->close();

    readFields(dict, cbcSwitch);
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readIfPresent()
{
    if
    (
        this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        WarningInFunction
            << "read option IOobject::MUST_READ or MUST_READ_IF_MODIFIED"
            << " suggests that a read constructor for field " << this->name()
            << " would be more appropriate." << endl;
    }
    else if
    (
        this->readOpt() == IOobject::READ_IF_PRESENT
     && this->template typeHeaderOk<GeometricField<Type, PatchField, GeoMesh>>
        (
            true
        )
    )
    {
        readFields(); //cbcSwitch set to true by default

        // Check compatibility between field and mesh
        if (this->size() != GeoMesh::size(this->mesh()))
        {
            FatalIOErrorInFunction(this->readStream(typeName))
                << "   number of field elements = " << this->size()
                << " number of mesh elements = "
                << GeoMesh::size(this->mesh())
                << exit(FatalIOError);
        }

        readOldTimeIfPresent();

        return true;
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readOldTimeIfPresent()
{
    // Read the old time field if present
    IOobject field0
    (
        this->name()  + "_0",
        this->time().timeName(),
        this->db(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        this->registerObject()
    );

    if
    (
        field0.template typeHeaderOk<GeometricField<Type, PatchField, GeoMesh>>
        (
            true
        )
    )
    {
        if (debug)
        {
            InfoInFunction << "Reading old time level for field"
                << endl << this->info() << endl;
        }

        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            field0,
            this->mesh()
        );

        // Ensure the old time field oriented flag is set to the parent's state
        // Note: required for backwards compatibility in case of restarting from
        // an old run where the oriented state may not have been set
        field0Ptr_->oriented() = this->oriented();

        field0Ptr_->timeIndex_ = timeIndex_ - 1;

        if (!field0Ptr_->readOldTimeIfPresent())
        {
            field0Ptr_->oldTime();
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType, actualPatchTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& name,
    const objectRegistry& db,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType,
    const wordList& actualPatchTypes
)
:
    GeometricField
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            db
        ),
        mesh,
        ds,
        patchFieldType,
        actualPatchTypes
    )
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType, actualPatchTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    boundaryField_.forceAssign(dt.value());

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    if (debug)
    {
        InfoInFunction << "Creating temporary" << endl << this->info() << endl;
    }

    boundaryField_.forceAssign(dt.value());

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& name,
    const objectRegistry& db,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Type& value,
    const word& patchFieldType,
    const wordList& actualPatchTypes
)
:
    GeometricField
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            db
        ),
        mesh,
        dimensioned<Type>(name, ds, value),
        patchFieldType,
        actualPatchTypes
    )
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Internal& diField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    Internal(io, diField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, ptfl)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    Internal(io, mesh, ds, iField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, ptfl)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from components" << endl << this->info() << endl;
    }

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const bool readOldTime,
    const bool cbcSwitch
)
:
    Internal(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary())
{
    readFields(cbcSwitch);

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorInFunction(this->readStream(typeName))
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    if (readOldTime)
    {
        readOldTimeIfPresent();
    }

    if (debug)
    {
        InfoInFunction
            << "Finishing read-construction of" << endl << this->info() << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    Istream& is
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dimless),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, DimensionedField<Type, GeoMesh>::readField(is))
{
    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorInFunction(is)
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    readOldTimeIfPresent();

    if (debug)
    {
        Info<< "Finishing read-construct of "
               "GeometricField<Type, PatchField, GeoMesh>"
            << endl << this->info() << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& dict,
    const bool cbcSwitch
)
:
    Internal(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary())
{
    readFields(dict, cbcSwitch);

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalErrorInFunction
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    if (debug)
    {
        InfoInFunction
            << "Finishing dictionary-construct of "
            << endl << this->info() << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy" << endl << this->info() << endl;
    }

    if (gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            *gf.field0Ptr_
        );
    }

    this->writeOpt() = IOobject::NO_WRITE;
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal
    (
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp" << endl << this->info() << endl;
    }

    this->writeOpt() = IOobject::NO_WRITE;

    tgf.clear();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_ && !this->name().endsWith("PrevIter"))
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting IO params"
            << endl << this->info() << endl;
    }

    tgf.clear();

    readIfPresent();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(newName, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting name"
            << endl << this->info() << endl;
    }

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            newName + "_0",
            *gf.field0Ptr_
        );
    }
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal
    (
        newName,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting name"
            << endl << this->info() << endl;
    }

    tgf.clear();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const word& patchFieldType
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params"
            << endl << this->info() << endl;
    }

    boundaryField_.forceAssign(gf.boundaryField_);

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes

)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    if (debug)
    {
        InfoInFunction
            << "Constructing as copy resetting IO params and patch types"
            << endl << this->info() << endl;
    }

    boundaryField_.forceAssign(gf.boundaryField_);

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


#ifndef NoConstructFromTmp
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal
    (
        io,
        const_cast<GeometricField<Type, PatchField, GeoMesh>&>(tgf()),
        tgf.isTmp()
    ),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    if (debug)
    {
        InfoInFunction
            << "Constructing from tmp resetting IO params and patch types"
            << endl << this->info() << endl;
    }

    boundaryField_.forceAssign(tgf().boundaryField_);

    tgf.clear();
}
#endif


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::cloneUnSliced() const
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name(),
                this->mesh().thisDb().time().timeName(),
                this->mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            GeometricField<Type, PatchField, GeoMesh>::Patch::calculatedType()
        )
    );
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::~GeometricField()
{
    deleteDemandDrivenData(field0Ptr_);
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal&
Foam::GeometricField<Type, PatchField, GeoMesh>::ref()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal::FieldType&
Foam::GeometricField<Type, PatchField, GeoMesh>::primitiveFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return boundaryField_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::isOldTime() const
{
    return
        this->name().size() > 2
     && this->name()(this->name().size()-2, 2) == "_0";
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTimes() const
{
    if (field0Ptr_ && timeIndex_ != this->time().timeIndex() && !isOldTime())
    {
        storeOldTime();
    }

    // Correct time index
    timeIndex_ = this->time().timeIndex();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTime() const
{
    if (field0Ptr_)
    {
        field0Ptr_->storeOldTime();

        if (debug)
        {
            InfoInFunction
                << "Storing old time field for field" << endl
                << this->info() << endl;
        }

        field0Ptr_->forceAssign(*this);
        field0Ptr_->timeIndex_ = timeIndex_;

        if (field0Ptr_->field0Ptr_)
        {
            field0Ptr_->writeOpt() = this->writeOpt();
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::label Foam::GeometricField<Type, PatchField, GeoMesh>::nOldTimes() const
{
    if (field0Ptr_)
    {
        return field0Ptr_->nOldTimes() + 1;
    }
    else
    {
        return 0;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime() const
{
    if (!field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + "_0",
                this->time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                this->registerObject()
            ),
            *this
        );

        if (debug)
        {
            InfoInFunction
                << "created old time field " << field0Ptr_->info() << endl;

            if (debug&2)
            {
                error::printStack(Info);
            }
        }
    }
    else
    {
        storeOldTimes();
    }

    return *field0Ptr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime()
{
    static_cast<const GeometricField<Type, PatchField, GeoMesh>&>(*this)
        .oldTime();

    return *field0Ptr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime(const label n) const
{
    if (n == 0)
    {
        return *this;
    }
    else
    {
        return oldTime().oldTime(n - 1);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime(const label n)
{
    if (n == 0)
    {
        return *this;
    }
    else
    {
        return oldTime().oldTime(n - 1);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storePrevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        if (debug)
        {
            InfoInFunction
                << "Allocating previous iteration field" << endl
                << this->info() << endl;
        }

        fieldPrevIterPtr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + "PrevIter",
                this->time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                this->registerObject()
            ),
            *this
        );
    }
    else
    {
        fieldPrevIterPtr_->forceAssign(*this);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::prevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorInFunction
            << "previous iteration field" << endl << this->info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}

// Return non-const access to previous iteration field
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::prevIter()
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorInFunction
            << "previous iteration field" << endl << this->info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::
updateBoundaryCoeffs()
{
    label currentState = this->eventNo();
    boundaryFieldRef().updateCoeffs();
    this->eventNo() = currentState;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::
correctBoundaryConditions()
{
    this->setUpToDate();
    storeOldTimes();
    boundaryField_.evaluate();
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::needReference() const
{
    // Search all boundary conditions, if any are
    // fixed-value or mixed (Robin) do not set reference level for solution.

    bool needRef = true;

    forAll(boundaryField_, patchi)
    {
        if (boundaryField_[patchi].fixesValue())
        {
            needRef = false;
            break;
        }
    }

    reduce(needRef, andOp<bool>());

    return needRef;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readInternalFields()
{
    const localIOdictionary dict
    (
        IOobject
        (
            this->name(),
            this->instance(),
            this->local(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        typeName
    );

    this->close();

    Internal::readField(dict, "internalField");
}




template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax(const scalar alpha)
{
    if (debug)
    {
        InfoInFunction
           << "Relaxing" << endl << this->info() << " by " << alpha << endl;
    }

    if (alpha != 1.0)
    {
        forceAssign(prevIter() + alpha*(*this - prevIter()));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax()
{
    word name = this->name();

    if
    (
        this->mesh().data::template lookupOrDefault<bool>
        (
            "finalIteration",
            false
        )
    )
    {
        name += "Final";
    }

    if (this->mesh().solution().relaxField(name))
    {
        relax(this->mesh().solution().fieldRelaxationFactor(name));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::word Foam::GeometricField<Type, PatchField, GeoMesh>::select
(
    bool final
) const
{
    if (final)
    {
        return this->name() + "Final";
    }
    else
    {
        return this->name();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::writeMinMax
(
    Ostream& os
) const
{
    os  << "min/max(" << this->name() << ") = "
        << Foam::min(*this).value() << ", "
        << Foam::max(*this).value()
        << endl;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax(word name)
{
    scalar alpha = 0;

    if (this->mesh().solution().relaxField(name))
    {
        alpha = this->mesh().solution().fieldRelaxationFactor(name);
    }

    if (alpha > 0)
    {
        relax(alpha);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::
writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::T() const
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> result
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + ".T()",
                this->instance(),
                this->db()
            ),
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::T(result.ref().primitiveFieldRef(), primitiveField());
    Foam::T(result.ref().boundaryFieldRef(), boundaryField());

    return result;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
Foam::GeometricField<Type, PatchField, GeoMesh>::component
(
    const direction d
) const
{
    tmp<GeometricField<cmptType, PatchField, GeoMesh>> Component
    (
        new GeometricField<cmptType, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + ".component(" + Foam::name(d) + ')',
                this->instance(),
                this->db()
            ),
            this->mesh(),
            this->dimensions()
        )
    );

    Foam::component(Component.ref().primitiveFieldRef(), primitiveField(), d);
    Foam::component(Component.ref().boundaryFieldRef(), boundaryField(), d);

    return Component;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
     >& gcf
)
{
    primitiveFieldRef().replace(d, gcf.primitiveField());
    boundaryFieldRef().replace(d, gcf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const dimensioned<cmptType>& ds
)
{
    primitiveFieldRef().replace(d, ds.value());
    boundaryFieldRef().replace(d, ds.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::max
(
    const dimensioned<Type>& dt
)
{
    Foam::max(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::max(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::min
(
    const dimensioned<Type>& dt
)
{
    Foam::min(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::min(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::negate()
{
    primitiveFieldRef().negate();
    boundaryFieldRef().negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, gf, "=");

    // Only assign field contents not ID

    ref() = gf();
    boundaryFieldRef() = gf.boundaryField();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    if (this == &(tgf()))
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "=");

    // Only assign field contents not ID

    this->dimensions() = gf.dimensions();
    this->oriented() = gf.oriented();

    if (tgf.isTmp())
    {
        // Transfer the storage from the tmp
        primitiveFieldRef().transfer
        (
            const_cast<Field<Type>&>(gf.primitiveField())
        );
    }
    else
    {
        primitiveFieldRef() = gf.primitiveField();
    }

    boundaryFieldRef() = gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    ref() = dt;
    boundaryFieldRef() = dt.value();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::forceAssign
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    const GeometricField<Type, PatchField, GeoMesh>& gf = tgf();

    checkField(*this, gf, "==");

    // Only assign field contents not ID

    ref() = gf();
    boundaryFieldRef().forceAssign(gf.boundaryField());

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::forceAssign
(
    const dimensioned<Type>& dt
)
{
    ref() = dt;
    boundaryFieldRef().forceAssign(dt.value());
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const GeometricField<TYPE, PatchField, GeoMesh>& gf                        \
)                                                                              \
{                                                                              \
    checkField(*this, gf, #op);                                                \
                                                                               \
    ref() op gf();                                                             \
    boundaryFieldRef() op gf.boundaryField();                                  \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const tmp<GeometricField<TYPE, PatchField, GeoMesh>>& tgf                  \
)                                                                              \
{                                                                              \
    operator op(tgf());                                                        \
    tgf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    ref() op dt;                                                               \
    boundaryFieldRef() op dt.value();                                          \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    //write boundary reference data if available
    if (gf.db().template foundObject<IOdictionary>(gf.name() + ".reference"))
    {
        gf.db().template lookupObject<IOdictionary>
        (
            gf.name() + ".reference"
        ).regIOobject::write();
    }

    gf().writeData(os, "internalField");
    os  << nl;
    gf.boundaryField().writeEntry("boundaryField", os);

    os.check(FUNCTION_NAME);
    return os;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    os << tgf();
    tgf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/GeometricFields/GeometricField/GeometricFieldFunctions.C"

// ************************************************************************* //
