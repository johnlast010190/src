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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fvMesh/fvPatches/derived/indirectWall/indirectWallFvPatch.H"
#include "fvMesh/fvPatches/derived/invisible/invisibleFvPatch.H"
#include "fvMesh/fvPatches/derived/inactive/inactiveFvPatch.H"
#include "fvMesh/fvPatches/constraint/cyclicAMI/cyclicAMIFvPatch.H"
#include "fvMesh/fvPatches/constraint/cyclicACMI/cyclicACMIFvPatch.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetryPlane/symmetryPlaneFvPatch.H"
#include "fvMesh/fvPatches/constraint/cyclic/cyclicFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/basic/generic/genericFvPatch.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "fvMesh/fvPatches/derived/nonConformalOrig/nonConformalOrigFvPatch.H"

#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "fields/fvPatchFields/constraint/cyclic/cyclicFvPatchField.H"
#include "fields/fvPatchFields/constraint/cyclicAMI/cyclicAMIFvPatchField.H"
#include "fields/fvPatchFields/constraint/cyclicACMI/cyclicACMIFvPatchField.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"
#include "fields/fvPatchFields/constraint/symmetry/symmetryFvPatchField.H"
#include "fields/fvPatchFields/constraint/symmetryPlane/symmetryPlaneFvPatchField.H"
#include "fields/fvPatchFields/constraint/wedge/wedgeFvPatchField.H"
#include "fields/fvPatchFields/constraint/nonConformalCyclic/nonConformalCyclicFvPatchField.H"
#include "fields/fvPatchFields/constraint/nonConformalError/nonConformalErrorFvPatchField.H"
#include "fields/fvPatchFields/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void boundarySetup<Type>::defaultBoundaries
(
    const fvMesh& mesh,
    const DimensionedField<Type, volMesh>& iField
)
{
    forAll(*this, patchI)
    {
        const fvPatch& cPatch(mesh.boundary()[patchI]);
        word patchType(cPatch.type());

        if
        (
            patchType == fvPatch::typeName
         || patchType == invisibleFvPatch::typeName
         || patchType == inactiveFvPatch::typeName
         || patchType == "indirectWall"
         || patchType == "inlet"
         || patchType == "outlet"
         || patchType == mappedFvPatch::typeName
         || patchType == nonConformalOrigFvPatch::typeName
        )
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    zeroGradientFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (isA<wallFvPatch>(cPatch))
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    zeroGradientFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if
        (
            isA<cyclicFvPatch>(cPatch)
        && !isA<nonConformalCyclicFvPatch>(cPatch)
        )
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    cyclicFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == cyclicAMIFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    cyclicAMIFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == cyclicACMIFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    cyclicACMIFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if
        (
            isA<processorFvPatch>(cPatch)
        && !isA<nonConformalProcessorCyclicFvPatch>(cPatch)
        )
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    processorFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == symmetryFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    symmetryFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == symmetryPlaneFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    symmetryPlaneFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == wedgeFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    wedgeFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == emptyFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    emptyPolyPatch::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == nonConformalCyclicFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    nonConformalCyclicFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == nonConformalErrorFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    nonConformalErrorFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == nonConformalProcessorCyclicFvPatch::typeName)
        {
            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    nonConformalProcessorCyclicFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else
        {
            WarningInFunction
                << cPatch.name() << ": patch type " << patchType
                << " not supported."
                << " Applying zeroGradient boundary conditions by default."
                << endl;

            this->set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    zeroGradientFvPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
    }

    // Default boundaries have been constructed without a dictionary,
    // leaving patch field uninitialised so far

    GeometricField<Type, fvPatchField, volMesh>::Boundary::evaluatePatchFields
    (
        mesh,
        *this
    );
}


template<class Type>
void boundarySetup<Type>::defaultBoundaries
(
    const fvMesh& mesh
)
{
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& cPatch(mesh.boundary()[patchI]);
        word patchType(cPatch.type());
        word patchName(cPatch.name());

        if
        (
            patchType == fvPatch::typeName
         || patchType == invisibleFvPatch::typeName
         || patchType == inactiveFvPatch::typeName
         || patchType == "indirectWall"
         || patchType == "inlet"
         || patchType == "outlet"
         || patchType == mappedFvPatch::typeName
         || patchType == nonConformalOrigFvPatch::typeName
        )
        {
            dictionary patchDict;

            patchDict.add("type", zeroGradientFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (isA<wallFvPatch>(cPatch))
        {
            dictionary patchDict;

            patchDict.add("type", zeroGradientFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if
        (
            isA<cyclicFvPatch>(cPatch)
        && !isA<nonConformalCyclicFvPatch>(cPatch)
        )
        {
            dictionary patchDict;

            patchDict.add("type", cyclicFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == cyclicAMIFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", cyclicAMIFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == cyclicACMIFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", cyclicACMIFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if
        (
            isA<processorFvPatch>(cPatch)
        && !isA<nonConformalProcessorCyclicFvPatch>(cPatch)
        )
        {
            dictionary patchDict;

            patchDict.add("type", processorFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == symmetryFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", symmetryFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == symmetryPlaneFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", symmetryPlaneFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == wedgeFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", wedgeFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == emptyFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", emptyPolyPatch::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == nonConformalCyclicFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add
            (
                "type",
                nonConformalCyclicFvPatchField<Type>::typeName
            );

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == nonConformalErrorFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add
            (
                "type",
                nonConformalErrorFvPatchField<Type>::typeName
            );

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == nonConformalProcessorCyclicFvPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add
            (
                "type",
                nonConformalProcessorCyclicFvPatchField<Type>::typeName
            );

            dictPtr_->add(patchName, patchDict);
        }
        else
        {
            WarningInFunction
                << cPatch.name() << ": patch type " << patchType
                << " not supported."
                << " Applying zeroGradient boundary conditions by default."
                << endl;

            dictionary patchDict;

            patchDict.add("type", zeroGradientFvPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
    }
}


template<class Type>
void boundarySetup<Type>::typeBoundaries
(
    const fvMesh& mesh,
    const DimensionedField<Type, volMesh>& iField,
    const dictionary& dict
)
{
    forAll(*this, patchI)
    {
        const fvPatch& cPatch(mesh.boundary()[patchI]);
        word patchType(cPatch.type());

        //physical type interpreter
        if
        (
            patchType == fvPatch::typeName
         || patchType == mappedFvPatch::typeName
        )
        {
            if (cPatch.patch().physicalType() == "inlet")
            {
                patchType = "inlet";
            }
            else if (cPatch.patch().physicalType() == "outlet")
            {
                patchType = "outlet";
            }
        }

        if (dict.found(patchType))
        {
            if (patchType != emptyPolyPatch::typeName)
            {
                this->set
                (
                    patchI,
                    fvPatchField<Type>::New
                    (
                        cPatch,
                        iField,
                        dict.subDict(patchType)
                    )
                );
            }
        }
    }
}

template<class Type>
void boundarySetup<Type>::typeBoundaries
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& cPatch(mesh.boundary()[patchI]);
        word patchType(cPatch.type());
        word patchName(cPatch.name());

        //physical type interpreter
        if
        (
            patchType == fvPatch::typeName
         || patchType == mappedFvPatch::typeName
        )
        {
            if (cPatch.patch().physicalType() == "inlet")
            {
                patchType = "inlet";
            }
            else if (cPatch.patch().physicalType() == "outlet")
            {
                patchType = "outlet";
            }
        }

        if (dict.found(patchType))
        {
            if (patchType != emptyPolyPatch::typeName)
            {
                dictPtr_->set(patchName, dict.subDict(patchType));
            }
        }
    }
}

template<class Type>
void boundarySetup<Type>::namedBoundaries
(
    const fvMesh& mesh,
    const DimensionedField<Type, volMesh>& iField,
    const dictionary& dict,
    bool strictPatchNameChecking
)
{

    boolList modifiedPatch(mesh.boundary().size(), false);


    //ref to boundary mesh
    const fvBoundaryMesh& bmesh(mesh.boundary());

    //Order of importance
    //1. Explicit
    //2. Group
    //3. RegExp
    //Where the same interface gets multiple hits, the last
    //entry will win.

    // 1. Handle explicit patch names.
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict() && !iter().keyword().isPattern())
        {
            word patchName = iter().keyword();
            label patchI = bmesh.findPatchID(patchName);

            if (patchI != -1)
            {
                const fvPatch& cPatch(bmesh[patchI]);
                word patchType(cPatch.type());

                if
                (
                    patchType != emptyPolyPatch::typeName
                    && patchType != wedgeFvPatch::typeName
                )
                {

                    this->set
                    (
                        patchI,
                        fvPatchField<Type>::New
                        (
                            cPatch,
                            iField,
                            iter().dict()
                        )
                    );

                    modifiedPatch[patchI] = true;
                }
            }
            else if (strictPatchNameChecking)
            {
                FatalErrorInFunction
                    << "No patch matching boundary "
                    << "name: " << iter().keyword()
                    << ". Disable strictPatchNameChecking to prevent"
                    << " this exception."
                    << exit(FatalError);
            }
        }
    }

    // 2. Patch-groups.
    if (dict.size())
    {
        for
        (
            IDLList<entry>::const_reverse_iterator iter = dict.rbegin();
            iter != dict.rend();
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

                if (!patchIDs.size() && strictPatchNameChecking)
                {
                    FatalErrorInFunction
                        << "No patch matching boundary "
                        << "name: " << iter().keyword()
                        << ". Disable strictPatchNameChecking to prevent"
                        << " this exception."
                        << exit(FatalError);
                }

                forAll(patchIDs, i)
                {
                    label patchI = patchIDs[i];

                    if (!modifiedPatch[patchI])
                    {

                        const fvPatch& cPatch(bmesh[patchI]);
                        word patchType(cPatch.type());

                        if
                        (
                            patchType != emptyPolyPatch::typeName
                            && patchType != processorFvPatch::typeName
                            && patchType != wedgeFvPatch::typeName
                        )
                        {
                            this->set
                            (
                                patchI,
                                fvPatchField<Type>::New
                                (
                                    cPatch,
                                    iField,
                                    iter().dict()
                                )
                            );
                            modifiedPatch[patchI] = true;
                        }

                    }
                }
            }
        }
    }

    // 3. Wildcard patch overrides
    if (dict.size())
    {
        forAllReverse(bmesh, patchI)
        {
            if (!modifiedPatch[patchI])
            {
                if (dict.found(bmesh[patchI].name()))
                {
                    const fvPatch& cPatch(bmesh[patchI]);
                    word patchType(cPatch.type());

                    if
                    (
                        patchType != emptyPolyPatch::typeName
                        && patchType != processorFvPatch::typeName
                        && patchType != wedgeFvPatch::typeName
                    )
                    {
                        this->set
                        (
                            patchI,
                            fvPatchField<Type>::New
                            (
                                cPatch,
                                iField,
                                dict.subDict(bmesh[patchI].name())
                            )
                        );
                        modifiedPatch[patchI] = true;
                    }
                }
            }
        }
    }
}


template<class Type>
void boundarySetup<Type>::namedBoundaries
(
    const fvMesh& mesh,
    const dictionary& dict,
    bool strictPatchNameChecking
)
{

    //ref to boundary mesh
    const fvBoundaryMesh& bmesh(mesh.boundary());

    //to exactly match instantiation order and items, this must match
    //the geometricField::Boundary constructor in terms of interpretation

    //explicit names && non-Group reg-exp
    forAllConstIter(dictionary, dict, iter)
    {
        const entry& e = iter();

        if (e.isDict())
        {
            const labelList patchIDs = bmesh.findIndices
            (
                e.keyword(),
                false                    // use patchGroups
            );

            //delete all patch name matches from the boundaryField
            forAll(patchIDs, i)
            {
                dictPtr_->remove(bmesh[patchIDs[i]].name());
            }
        }
    }

    //groups
    forAllConstIter(dictionary, dict, iter)
    {
        const entry& e = iter();

        if (e.isDict() && !e.keyword().isPattern())
        {
            const labelList patchIDs = bmesh.findIndices
            (
                e.keyword(),
                true                    // use patchGroups
            );

            //delete all patch name matches from the boundaryField
            forAll(patchIDs, i)
            {
                dictPtr_->remove(bmesh[patchIDs[i]].name());
            }
        }
    }

    //merge input dict
    dictPtr_->merge(dict);

}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct dictionary only
template<class Type>
boundarySetup<Type>::boundarySetup
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict,
    bool strictPatchNameChecking,
    bool printBoundaries
)
:
    PtrList<fvPatchField<Type>> (0)
{
    dictPtr_.reset(new dictionary);

    defaultBoundaries(mesh);

    if (dict.found("boundaryTypeDefaults"))
    {
        typeBoundaries
        (
            mesh,
            dict.subDict("boundaryTypeDefaults")
        );
    }

    if (dict.found("boundaryConditions"))
    {
        namedBoundaries
        (
            mesh,
            dict.subDict("boundaryConditions"),
            strictPatchNameChecking
        );
    }

    if (printBoundaries)
    {
        Info<< nl << name << endl;
        Info<< token::BEGIN_BLOCK << incrIndent << nl;
        Info<< indent << dictPtr_() << endl;
        Info  << decrIndent << token::END_BLOCK << nl << endl;
    }
}


// Construct from components
template<class Type>
boundarySetup<Type>::boundarySetup
(
    const fvMesh& mesh,
    const DimensionedField<Type, volMesh>& iField,
    const dictionary& dict,
    bool strictPatchNameChecking,
    bool printBoundaries
)
:
    PtrList<fvPatchField<Type>> (mesh.boundary().size())
{
    dictPtr_.reset(new dictionary);

    defaultBoundaries(mesh, iField);

    if (dict.found("boundaryTypeDefaults"))
    {
        typeBoundaries
        (
            mesh,
            iField,
            dict.subDict("boundaryTypeDefaults")
        );
    }
    if (dict.found("boundaryConditions"))
    {
        namedBoundaries
        (
            mesh,
            iField,
            dict.subDict("boundaryConditions"),
            strictPatchNameChecking
        );
    }

    if (printBoundaries)
    {
        Info<< nl << iField.name() << endl;

        Info<< token::BEGIN_BLOCK << incrIndent << nl;

        forAll(*this, patchi)
        {
            Info<< indent << this->operator[](patchi).patch().name() << nl
                 << indent << token::BEGIN_BLOCK << nl
                 << incrIndent << this->operator[](patchi) << decrIndent
                 << indent << token::END_BLOCK << endl;
        }

        Info  << decrIndent << token::END_BLOCK << nl << endl;
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy list elements
template<class Type>
boundarySetup<Type>::~boundarySetup()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const dictionary& boundarySetup<Type>::boundaryDict() const
{
    if (dictPtr_.valid())
    {
        return dictPtr_();
    }
    else
    {
        FatalErrorInFunction << "Boundary dictionary empty"
            << abort(FatalError);

        return dictPtr_();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

// ************************************************************************* //
