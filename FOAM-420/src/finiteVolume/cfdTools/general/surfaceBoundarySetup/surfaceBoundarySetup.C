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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/surfaceBoundarySetup/surfaceBoundarySetup.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
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

#include "fields/fvsPatchFields/basic/calculated/calculatedFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/cyclic/cyclicFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/cyclicAMI/cyclicAMIFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/cyclicACMI/cyclicACMIFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/processor/processorFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/symmetry/symmetryFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/symmetryPlane/symmetryPlaneFvsPatchField.H"
#include "fields/fvsPatchFields/constraint/wedge/wedgeFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void surfaceBoundarySetup<Type>::defaultBoundaries
(
    const surfaceMesh& mesh,
    const DimensionedField<Type, surfaceMesh>& iField
)
{
    forAll(*this, patchI)
    {
        const fvPatch& cPatch(mesh().boundary()[patchI]);
        word patchType(cPatch.type());

        if
        (
            patchType == fvPatch::typeName
         || patchType == "inlet"
         || patchType == "outlet"
         || patchType == mappedFvPatch::typeName
        )
        {
            this->set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    calculatedFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    calculatedFvsPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (isA<cyclicFvPatch>(cPatch))
        {
            this->set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    cyclicFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    cyclicAMIFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    cyclicACMIFvsPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (isA<processorFvPatch>(cPatch))
        {
            this->set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    processorFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    symmetryFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    symmetryPlaneFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    wedgeFvsPatchField<Type>::typeName,
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
                fvsPatchField<Type>::New
                (
                    emptyPolyPatch::typeName,
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
                << " Applying calculated boundary conditions by default."
                << endl;

            this->set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    calculatedFvsPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
    }
}

template<class Type>
void surfaceBoundarySetup<Type>::typeBoundaries
(
    const surfaceMesh& mesh,
    const DimensionedField<Type, surfaceMesh>& iField,
    const dictionary& dict
)
{
    forAll(*this, patchI)
    {
        const fvPatch& cPatch(mesh().boundary()[patchI]);
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
                    fvsPatchField<Type>::New
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
void surfaceBoundarySetup<Type>::namedBoundaries
(
    const surfaceMesh& mesh,
    const DimensionedField<Type, surfaceMesh>& iField,
    const dictionary& dict,
    bool strictPatchNameChecking
)
{

    boolList modifiedPatch(mesh().boundary().size(), false);


    //ref to boundary mesh
    const fvBoundaryMesh& bmesh(mesh().boundary());

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
                        fvsPatchField<Type>::New
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
                                fvsPatchField<Type>::New
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
                            fvsPatchField<Type>::New
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





// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //


// Construct from components
template<class Type>
surfaceBoundarySetup<Type>::surfaceBoundarySetup
(
    const surfaceMesh& mesh,
    const DimensionedField<Type, surfaceMesh>& iField,
    const dictionary& dict,
    bool strictPatchNameChecking
)
:
    PtrList<fvsPatchField<Type>> (mesh().boundary().size())
{

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

}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy list elements
template<class Type>
surfaceBoundarySetup<Type>::~surfaceBoundarySetup()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

// ************************************************************************* //
