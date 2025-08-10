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
    (c) 2016 OpenCFD Ltd.
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/codedFixedValue/codedFixedValueFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCode.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCodeContext.H"
#include "primitives/strings/stringOps/stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeTemplateC
    = "fixedValueFvPatchFieldTemplate.C";

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeTemplateH
    = "fixedValueFvPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::setFieldTemplates
(
    dynamicCode& dynCode
)
{
    word fieldType(pTraits<Type>::typeName);

    // template type for fvPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for fvPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::IOdictionary& Foam::codedFixedValueFvPatchField<Type>::dict() const
{
    const objectRegistry& obr = this->db();

    if (obr.foundObject<IOdictionary>("codeDict"))
    {
        return obr.lookupObject<IOdictionary>("codeDict");
    }
    else
    {
        return obr.store
        (
            new IOdictionary
            (
                IOobject
                (
                    "codeDict",
                    this->db().time().system(),
                    this->db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


template<class Type>
Foam::dlLibraryTable& Foam::codedFixedValueFvPatchField<Type>::libs() const
{
    return const_cast<dlLibraryTable&>(this->db().time().libs());
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // take no chances - typeName must be identical to name_
    dynCode.setFilterVariable("typeName", name_);

    // set TemplateType and FieldType filter variables
    // (for fvPatchField)
    setFieldTemplates(dynCode);

    // compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // copy filtered H template
    dynCode.addCopyFile(codeTemplateH);


    // debugging: make BC verbose
    //  dynCode.setFilterVariable("verbose", "true");
    //  Info<<"compile " << name_ << " sha1: "
    //      << context.sha1() << endl;

    // define Make/options
    dynCode.setOptions
    (
        "-I$(FOAM_PROJECT_DIR)/src/finiteVolume "
        "-I$(FOAM_PROJECT_DIR)/src/meshTools " +
        context.options()
    );
    dynCode.setLinkLibs
    (
        "-lOpenFOAM -lfiniteVolume -lmeshTools " + context.libs()
    );
}


template<class Type>
const Foam::dictionary& Foam::codedFixedValueFvPatchField<Type>::codeDict()
const
{
    // use system/codeDict or in-line
    return
    (
        dict_.found("code")
      ? dict_
      : this->dict().subDict(name_)
    );
}


template<class Type>
Foam::string Foam::codedFixedValueFvPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->internalField().name();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::clearRedirect() const
{
    // remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    codedBase(),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    codedBase(),
    dict_(dict),
    name_
    (
        dict.found("redirectType")
      ? dict.lookup("redirectType")
      : dict.lookup("name")
    ),
    redirectPatchFieldPtr_()
{
    updateLibrary(name_);
}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fvPatchField<Type>&
Foam::codedFixedValueFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OStringStream os;
        os.writeEntry("type", name_);
        static_cast<const Field<Type>&>(*this).writeEntry("value", os);
        IStringStream is(os.str());
        dictionary dict(is);

        redirectPatchFieldPtr_.set
        (
            fvPatchField<Type>::New
            (
                this->patch(),
                this->internalField(),
                dict
            ).ptr()
        );
    }
    return redirectPatchFieldPtr_();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through value
    this->forceAssign(fvp);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).evaluate(commsType);

    fixedValueFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    os.writeEntry("name", name_);

    codedBase::writeCodeDict(os, dict_);
}


// ************************************************************************* //
