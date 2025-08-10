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

#include "cfdTools/general/pointBoundarySetup/pointBoundarySetup.H"
#include "meshes/pointMesh/pointPatches/derived/wall/wallPointPatch.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatch/cyclicAMIPointPatch.H"
#include "AMIInterpolation/patches/cyclicACMI/cyclicACMIPointPatch/cyclicACMIPointPatch.H"
#include "mappedPatches/mappedPointPatch/mappedPointPatch.H"
#include "meshes/pointMesh/pointPatches/constraint/processor/processorPointPatch.H"
#include "meshes/pointMesh/pointPatches/constraint/symmetry/symmetryPointPatch.H"
#include "meshes/pointMesh/pointPatches/constraint/symmetryPlane/symmetryPlanePointPatch.H"
#include "meshes/pointMesh/pointPatches/constraint/cyclic/cyclicPointPatch.H"
#include "meshes/pointMesh/pointPatches/constraint/empty/emptyPointPatch.H"
#include "meshes/pointMesh/pointPatches/constraint/wedge/wedgePointPatch.H"
#include "meshes/pointMesh/pointPatches/basic/generic/genericPointPatch.H"
#include "mappedPatches/mappedPointPatch/mappedPointPatch.H"
#include "regionCoupled/patches/regionCoupledPointPatch/regionCoupledPointPatch.H"

#include "fields/pointPatchFields/basic/zeroGradient/zeroGradientPointPatchField.H"
#include "fields/pointPatchFields/constraint/cyclic/cyclicPointPatchField.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatchField/cyclicAMIPointPatchField.H"
#include "AMIInterpolation/patches/cyclicACMI/cyclicACMIPointPatchField/cyclicACMIPointPatchField.H"
#include "fields/pointPatchFields/constraint/processor/processorPointPatchField.H"
#include "fields/pointPatchFields/constraint/symmetry/symmetryPointPatchField.H"
#include "fields/pointPatchFields/constraint/symmetryPlane/symmetryPlanePointPatchField.H"
#include "fields/pointPatchFields/constraint/wedge/wedgePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void pointBoundarySetup<Type>::defaultBoundaries
(
    const pointMesh& mesh
)
{
    forAll(mesh.boundary(), patchI)
    {
        const pointPatch& cPatch(mesh.boundary()[patchI]);
        const polyPatch& pPatch(mesh().boundaryMesh()[patchI]);
        word physicalType(pPatch.physicalType());
        word patchType(cPatch.type());
        word patchName(cPatch.name());

        if
        (
            patchType == pointPatch::typeName
         || patchType == "invisible"
         || patchType == "inactive"
         || patchType == "indirectWall"
         || patchType == "inlet"
         || physicalType == "inlet"
         || patchType == "outlet"
         || physicalType == "outlet"
         || patchType == mappedPointPatch::typeName
        )
        {
            dictionary patchDict;

            patchDict.add("type", zeroGradientPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (isA<wallPointPatch>(cPatch))
        {
            dictionary patchDict;

            patchDict.add("type", zeroGradientPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (isA<cyclicPointPatch>(cPatch))
        {
            dictionary patchDict;

            patchDict.add("type", cyclicPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == cyclicAMIPointPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", cyclicAMIPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == cyclicACMIPointPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", cyclicACMIPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (isA<processorPointPatch>(cPatch))
        {
            dictionary patchDict;

            patchDict.add("type", processorPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == symmetryPointPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", symmetryPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == symmetryPlanePointPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", symmetryPlanePointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == wedgePointPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", wedgePointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
        else if (patchType == emptyPointPatch::typeName)
        {
            dictionary patchDict;

            patchDict.add("type", emptyPolyPatch::typeName);

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

            patchDict.add("type", zeroGradientPointPatchField<Type>::typeName);

            dictPtr_->add(patchName, patchDict);
        }
    }
}

template<class Type>
void pointBoundarySetup<Type>::defaultBoundaries
(
    const pointMesh& mesh,
    const DimensionedField<Type, pointMesh>& iField
)
{
    forAll(*this, patchI)
    {
        const pointPatch& cPatch(mesh.boundary()[patchI]);
        word patchType(cPatch.type());

        if
        (
            patchType == facePointPatch::typeName
         || patchType == "indirectWall"
         || patchType == "inlet"
         || patchType == "outlet"
         || patchType == regionCoupledPointPatch::typeName
         || patchType == mappedPointPatch::typeName
        )
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    zeroGradientPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (isA<wallPointPatch>(cPatch))
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    zeroGradientPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (isA<cyclicPointPatch>(cPatch))
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    cyclicPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == cyclicAMIPointPatch::typeName)
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    cyclicAMIPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == cyclicACMIPointPatch::typeName)
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    cyclicACMIPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (isA<processorPointPatch>(cPatch))
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    processorPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == symmetryPointPatch::typeName)
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    symmetryPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == symmetryPlanePointPatch::typeName)
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    symmetryPlanePointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == wedgePointPatch::typeName)
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    wedgePointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
        else if (patchType == emptyPointPatch::typeName)
        {
            this->set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    emptyPointPatch::typeName,
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
                pointPatchField<Type>::New
                (
                    zeroGradientPointPatchField<Type>::typeName,
                    cPatch,
                    iField
                )
            );
        }
    }
}


template<class Type>
void pointBoundarySetup<Type>::typeBoundaries
(
    const pointMesh& mesh,
    const dictionary& dict
)
{
    forAll(mesh.boundary(), patchI)
    {
        const pointPatch& cPatch(mesh.boundary()[patchI]);
        const polyPatch& pPatch(mesh().boundaryMesh()[patchI]);
        word patchType(cPatch.type());
        word patchName(cPatch.name());

        //physical type interpreter
        if
        (
            patchType == fvPatch::typeName
         || patchType == mappedFvPatch::typeName
        )
        {
            if (pPatch.physicalType() == "inlet")
            {
                patchType = "inlet";
            }
            else if (pPatch.physicalType() == "outlet")
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
void pointBoundarySetup<Type>::typeBoundaries
(
    const pointMesh& mesh,
    const DimensionedField<Type, pointMesh>& iField,
    const dictionary& dict
)
{
    forAll(*this, patchI)
    {
        const pointPatch& cPatch(mesh.boundary()[patchI]);
        const polyPatch& pPatch(mesh().boundaryMesh()[patchI]);
        word patchType(cPatch.type());

        //physical type interpreter
        if
        (
            patchType == fvPatch::typeName
         || patchType == mappedFvPatch::typeName
        )
        {
            if (pPatch.physicalType() == "inlet")
            {
                patchType = "inlet";
            }
            else if (pPatch.physicalType() == "outlet")
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
                    pointPatchField<Type>::New
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
void pointBoundarySetup<Type>::namedBoundaries
(
    const pointMesh& mesh,
    const dictionary& dict,
    bool strictPatchNameChecking
)
{

    //ref to boundary mesh
    const pointBoundaryMesh& bmesh(mesh.boundary());

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


template<class Type>
void pointBoundarySetup<Type>::namedBoundaries
(
    const pointMesh& mesh,
    const DimensionedField<Type, pointMesh>& iField,
    const dictionary& dict,
    bool strictPatchNameChecking
)
{

    boolList modifiedPatch(mesh.boundary().size(), false);


    //ref to boundary mesh
    const polyBoundaryMesh& bmesh(mesh().boundaryMesh());
    const pointBoundaryMesh& pbmesh(mesh.boundary());

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
                const pointPatch& cPatch(pbmesh[patchI]);
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
                        pointPatchField<Type>::New
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

                        const pointPatch& cPatch(pbmesh[patchI]);
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
                                pointPatchField<Type>::New
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
                    const pointPatch& cPatch(pbmesh[patchI]);
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
                            pointPatchField<Type>::New
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

template<class Type>
pointBoundarySetup<Type>::pointBoundarySetup
(
    const pointMesh& mesh,
    const word& name,
    const dictionary& dict,
    bool strictPatchNameChecking
)
:
    PtrList<pointPatchField<Type>> (0)
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
}




// Construct from components
template<class Type>
pointBoundarySetup<Type>::pointBoundarySetup
(
    const pointMesh& mesh,
    const DimensionedField<Type, pointMesh>& iField,
    const dictionary& dict,
    bool strictPatchNameChecking
)
:
    PtrList<pointPatchField<Type>> (mesh.boundary().size())
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

}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy list elements
template<class T>
pointBoundarySetup<T>::~pointBoundarySetup()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const dictionary& pointBoundarySetup<Type>::boundaryDict() const
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
