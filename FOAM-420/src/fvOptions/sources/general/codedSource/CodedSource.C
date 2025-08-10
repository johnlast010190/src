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
    (c) 2016 OpenCFD Ltd.
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "CodedSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCode.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCodeContext.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::CodedSource<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    word sourceType(pTraits<Type>::typeName);

    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);
    dynCode.setFilterVariable("TemplateType", sourceType);
    dynCode.setFilterVariable("SourceType", sourceType + "Source");

    //dynCode.removeFilterVariable("code");
    dynCode.setFilterVariable("codeCorrect", codeCorrect_);
    dynCode.setFilterVariable("codeAddSup", codeAddSup_);
    dynCode.setFilterVariable("codeSetValue", codeSetValue_);

    // compile filtered C template
    dynCode.addCompileFile("codedFvOptionTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("codedFvOptionTemplate.H");

    // debugging: make BC verbose
    //         dynCode.setFilterVariable("verbose", "true");
    //         Info<<"compile " << name_ << " sha1: "
    //             << context.sha1() << endl;

    // define Make/options
    dynCode.setOptions
    (
        "-I$(FOAM_PROJECT_DIR)/src/finiteVolume \
        -I$(FOAM_PROJECT_DIR)/src/meshTools \
        -I$(FOAM_PROJECT_DIR)/src/sampling \
        -I$(FOAM_PROJECT_DIR)/src/fvOptions " + context.options()
    );
    dynCode.setLinkLibs
    (
        "-lOpenFOAM -lfiniteVolume -lmeshTools -lfvOptions -lsampling" +
        context.libs()
    );
}


template<class Type>
Foam::dlLibraryTable& Foam::fv::CodedSource<Type>::libs() const
{
    return const_cast<Time&>(mesh_.time()).libs();
}


template<class Type>
Foam::string Foam::fv::CodedSource<Type>::description() const
{
    return "fvOption:: " + name_;
}


template<class Type>
void Foam::fv::CodedSource<Type>::clearRedirect() const
{
    redirectFvOptionPtr_.clear();
}


template<class Type>
const Foam::dictionary& Foam::fv::CodedSource<Type>::codeDict() const
{
    return coeffs_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::CodedSource<Type>::CodedSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fv::option& Foam::fv::CodedSource<Type>::redirectFvOption() const
{
    if (!redirectFvOptionPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);
        constructDict.changeKeyword(modelType_ & "Coeffs", name_ & "Coeffs");

        redirectFvOptionPtr_ = option::New
        (
            name_,
            constructDict,
            mesh_
        );
    }
    return redirectFvOptionPtr_();
}


template<class Type>
void Foam::fv::CodedSource<Type>::correct
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::correct for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().correct(field);
}


template<class Type>
void Foam::fv::CodedSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().addSup(eqn, fieldi);
}


template<class Type>
void Foam::fv::CodedSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().addSup(rho, eqn, fieldi);
}


template<class Type>
void Foam::fv::CodedSource<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::constrain for source " << name_ << endl;
    }

    updateLibrary(name_);
    redirectFvOption().constrain(eqn, fieldi);
}


// ************************************************************************* //
