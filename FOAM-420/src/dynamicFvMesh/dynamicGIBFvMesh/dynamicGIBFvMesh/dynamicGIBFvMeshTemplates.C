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
    (c) 2015-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/syncTools/syncTools.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void StoreOldTimeFields
(
    const fvMesh& mesh
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());
        field.storeOldTimes();

    }
}

template
<
    class Type
>
void CorrectBCs
(
    const fvMesh& mesh
)
{
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, fvPatchField, volMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        if
        (
            (field.name() == "k") ||
            (field.name() == "epsilon")
        )
        {
            forAll(field, cI)
            {
                if (field[cI] < pTraits<Type>::zero)
                {
                    field[cI] = pTraits<Type>::zero;
                }
            }
            field.correctBoundaryConditions();
        }
        else
        {
            field.correctBoundaryConditions();
        }
/*
        if
        (
            (field.name() == "omega")
        )
        {
            forAll(field, cI)
            {
                if (field[cI] <= pTraits<Type>::zero)
                {
                    field[cI] = pTraits<Type>::zero;
                }
            }
            field.correctBoundaryConditions();
        }
*/
    }
}

template
<
    class Type
>
void CorrectProcessorBCs
(
    const fvMesh& mesh
)
{
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, fvPatchField, volMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        label nReq = Pstream::nRequests();

        forAll(field.boundaryField(), patchI)
        {
            fvPatchField<Type>& pf = field.boundaryFieldRef()[patchI];
            if (isA<processorFvPatchField<Type>>(pf))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
            }
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(field.boundaryField(), patchI)
        {
            fvPatchField<Type>& pf = field.boundaryFieldRef()[patchI];
            if (isA<processorFvPatchField<Type>>(pf))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void dynamicGIBFvMesh::PopFields
(
    const labelList& fPop,
    const surfaceScalarField& phiPop,
    const boolList& popIndi,
    const bool& popType
)
{
    const label dbgCell = -1;
    const label dbgProc = -1;
    const string dbgField = "";

    const boolList& popUpC = this->popUpCells();
    boolList popUpCNbr;
    syncTools::swapBoundaryCellList(*this, popUpC, popUpCNbr);

    const labelList& cReg = cRegion();
    labelList cRegNbr;
    syncTools::swapBoundaryCellList(*this, cReg, cRegNbr);

    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        this->thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());
        //label nT = field.nOldTimes();

        // Have not synced boundary values yet
        List<Type> fieldNbr;
        syncTools::swapBoundaryCellList(*this, field.primitiveField(), fieldNbr);

        const scalarField& cV0 = this->V0();
        scalarField V0 = cV0;

        if (Pstream::myProcNo() == dbgProc && (!dbgField.size() || field.name() == dbgField))
        {
            Pout<< "Before PopFields: " << field.name() << "[" << dbgCell << "] = " << field[dbgCell] << endl;
        }
        if (true)
        {
            for (int iter=0; iter < popSubSteps_; iter++)
            {
                scalarField V0Start = V0;
                GeometricField<Type, PatchField, GeoMesh> fieldStart = field;
                // Midpoint rule
                for (label m = 0; m < 2; m++)
                {
                    scalar dd = 1.0/popSubSteps_;
                    if (m == 0)
                    {
                        dd /= 2;
                    }
                    Field<Type> df = Field<Type>
                    (
                        this->cells().size(), pTraits<Type>::zero
                    );
                    Type ff = pTraits<Type>::zero;
                    Field<Type> minff = Field<Type>
                    (
                        this->cells().size(), pTraits<Type>::one*GREAT
                    );
                    Field<Type> maxff = Field<Type>
                    (
                        this->cells().size(), -pTraits<Type>::one*GREAT
                    );

                    scalarField sV = scalarField(this->cells().size(), 0.0);

                    forAll(fPop, fI)
                    {
                        const label& fPopI = fPop[fI];
                        if (fPopI<this->nInternalFaces())
                        {
                            const label& on = this->owner()[fPopI];
                            const label& nb = this->neighbour()[fPopI];

                            //ff = 0.5*(field[on]+field[nb]);
                            ff = pTraits<Type>::zero;;
                            scalar sVf = dd*phiPop[fPopI]*this->time().deltaTValue();

                            if (cRegion()[on] == cRegion()[nb])
                            {
                                if (popUpC[on]==1 && popUpC[nb]!=1)
                                {
                                    ff = field[nb];
                                }
                                else if (popUpC[nb]==1 && popUpC[on]!=1)
                                {
                                    ff = field[on];
                                }
                                else if (popUpC[nb]==1 && popUpC[on]==1)
                                {
                                    ff = 0.5*(field[on]+field[nb]);
                                }
                                else if (popUpC[nb]!=1 && popUpC[on]!=1)
                                {
                                    ff = 0.5*(field[on]+field[nb]);
                                }

                                df[on] += ff*sVf;
                                df[nb] -= ff*sVf;

                                if (boundPopValues_)
                                {
                                    minff[on] = min(minff[on], ff);
                                    maxff[on] = max(maxff[on], ff);
                                    minff[nb] = min(minff[nb], ff);
                                    maxff[nb] = max(maxff[nb], ff);
                                }
                            }
                            else
                            {
                                // Don't try to conserve across the regions
                                // We should perhaps do something to ensure global
                                // conservation in each region, but allowing it
                                // to flow between them seems to make no sense.
                                // Setting it to zero here would do this but
                                // locally would probably produce spikes.
                                df[on] += field[on]*sVf;
                                df[nb] -= field[nb]*sVf;

                                if (boundPopValues_)
                                {
                                    minff[on] = min(minff[on], field[on]);
                                    maxff[on] = max(maxff[on], field[on]);
                                    minff[nb] = min(minff[nb], field[nb]);
                                    maxff[nb] = max(maxff[nb], field[nb]);
                                }
                            }

                            if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && (on == dbgCell || nb == dbgCell))
                            {
                                Pout<< on << " " << nb << " face: " << popUpC[on] << " " << popUpC[nb] << " " << cRegion()[on] << " " << cRegion()[nb] << " " << field[on] << " " << field[nb] << " " << sVf << endl;
                            }

                            sV[on] += sVf;
                            sV[nb] -= sVf;
                        }
                        else
                        {
                            const label& on = this->faceOwner()[fPopI];
                            const label patchI = this->boundaryMesh().whichPatch(fPopI);

                            if (!isA<indirectPolyPatch>(this->boundary()[patchI].patch()))
                            {
                                const label lpfI = fPopI - this->boundaryMesh()[patchI].start();
                                const label bfI = fPopI - this->nInternalFaces();
                                const scalar sVf = dd*
                                    (
                                        phiPop.boundaryField()[patchI][lpfI]
                                    )*
                                    this->time().deltaTValue();

                                if (this->boundary()[patchI].coupled())
                                {
                                    ff = pTraits<Type>::zero;

                                    if (cRegion()[on] == cRegNbr[bfI])
                                    {
                                        if (popUpC[on]==1 && popUpCNbr[bfI]!=1)
                                        {
                                            ff = fieldNbr[bfI];
                                        }
                                        else if (popUpCNbr[bfI]==1 && popUpC[on]!=1)
                                        {
                                            ff = field[on];
                                        }
                                        else if (popUpCNbr[bfI]!=1 && popUpC[on]!=1)
                                        {
                                            ff = 0.5*(field[on]+fieldNbr[bfI]);
                                        }
                                        else if (popUpCNbr[bfI]==1 && popUpC[on]==1)
                                        {
                                            ff = 0.5*(field[on]+fieldNbr[bfI]);
                                        }

                                        if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && on == dbgCell)
                                        {
                                            Pout<< on << " " << " pface: " << fPopI << " " << popUpC[on] << " " << popUpCNbr[bfI] << " " << cRegion()[on] << " " << cRegNbr[bfI] << " " << field[on] << " " << fieldNbr[bfI] << " " << ff << " " << ff*sVf << endl;
                                        }

                                        df[on] += ff*sVf;

                                        if (boundPopValues_)
                                        {
                                            minff[on] = min(minff[on], ff);
                                            maxff[on] = max(maxff[on], ff);
                                        }
                                    }
                                }
                                else
                                {
                                    if (popUpC[on]==1)
                                    {
                                        if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && on == dbgCell)
                                        {
                                            Pout<< on << " bface: " << popUpC[on] << endl;
                                        }

                                        ff = field[on];

                                        df[on] += ff*sVf;

                                        if (boundPopValues_)
                                        {
                                            minff[on] = min(minff[on], ff);
                                            maxff[on] = max(maxff[on], ff);
                                        }

                                        sV[on] += sVf;
                                    }
                                }
                            }
                        }
                    }
                    if (m == 1)
                    {
                        V0 = V0Start;
                        field.ref() = fieldStart.ref();
                    }
                    forAll(popIndi, cI)
                    {
                        if (popIndi[cI])
                        {
                            scalar newVol = V0[cI]+sV[cI];

                            if (mag(newVol) > SMALL)
                            {
                                field[cI] = fieldStart[cI]*V0[cI] + df[cI];

                                if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && cI == dbgCell)
                                {
                                    Pout<< cI << " " << field[cI]/newVol << " " << fieldStart[cI] << " " << V0[cI] << " " << newVol << endl;
                                }

                                V0[cI] = newVol;
                                field[cI] /= V0[cI];
                            }
                            else
                            {
                              field[cI] = pTraits<Type>::zero;
                              V0[cI] = newVol;
                            }

                            if (boundPopValues_)
                            {
                                // Clamp the value so that it can't exceed the min/max
                                // of face values. Can occur if overall change in volume small
                                // compared to sum of abs values of face volume sweeps - then we
                                // end up with a large imbalance of incoming/outgoing fluxes due to slight
                                // variations in face values.
                                // TODO: redistribute conservatively rather
                                if (minff[cI] <= maxff[cI]) // Failsafe in case no faces had any influence
                                {
                                    field[cI] = max(minff[cI], field[cI]);
                                    field[cI] = min(maxff[cI], field[cI]);
                                    if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && cI == dbgCell)
                                    {
                                        Pout<< cI << " " << minff[cI] << " " << maxff[cI] << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if
        (
            debugMode_ &&
            (
                field.name() == "T" ||
                field.name() == "rho"
            )
        )
        {
            forAll(field, i)
            {
                if (refCast<scalarField>(field)[i] < 0)
                {
                    Pout<< "Neg " << field.name() << " after pop: " << i << " " << field[i] << endl;
                }
            }
        }
    }
}



template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void dynamicGIBFvMesh::PopFields2
(
    const labelList& fPop,
    const surfaceScalarField& phiPop,
    const boolList& popIndi,
    const scalarField& oldVolgrow
)
{
    const label dbgCell = -1;
    const label dbgProc = -1;
    const string dbgField = "";

    const boolList& popUpC = this->popUpCells();
    boolList popUpCNbr;
    syncTools::swapBoundaryCellList(*this, popUpC, popUpCNbr);

    const labelList& cReg = cRegion();
    labelList cRegNbr;
    syncTools::swapBoundaryCellList(*this, cReg, cRegNbr);

    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        this->thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());
        //label nT = field.nOldTimes();

        if (Pstream::myProcNo() == dbgProc && (!dbgField.size() || field.name() == dbgField))
        {
            Pout<< "Before PopFields2: " << field.name() << "[" << dbgCell << "] = " << field[dbgCell] << endl;
        }

//        forAll(popUpC, cI)
//        {
//            if (popUpC[cI] == 1)
//            {
//                field[cI] = pTraits<Type>::zero;
//            }
//        }

        // Initialise the faces of the pop cells by sweeping values
        // from non-pop cells in the same region (should maybe use FaceCellWave?)
        boolList popUpToDo(popUpC);
        while (true)
        {
            bool somethingChanged = false;
            bool allDone = true;
            boolList popUpToDoNext(popUpToDo);
            boolList popUpToDoNbr;
            syncTools::swapBoundaryCellList(*this, popUpToDo, popUpToDoNbr);
            // Have not synced boundary values yet
            List<Type> fieldNbr;
            syncTools::swapBoundaryCellList(*this, field.primitiveField(), fieldNbr);
            forAll(popUpToDo, cI)
            {
                if (popUpToDo[cI])
                {
                    //Pout<< cI << endl;
                    label count = 0;
                    Type val = pTraits<Type>::zero;
                    forAll(this->cells()[cI], cfI)
                    {
                        const label fI = this->cells()[cI][cfI];
                        if (fI < this->nInternalFaces())
                        {
                            const label on = this->owner()[fI];
                            const label nb = this->neighbour()[fI];
                            if (field.name() == dbgField && cI == dbgCell && Pstream::myProcNo() == dbgProc)
                            {
                                Pout<< fI << " " << cReg[on] << " " << cReg[nb] << " " << popUpToDo[on] << " " << popUpToDo[nb] << endl;
                            }
                            if (on == cI && !popUpToDo[nb] && cReg[nb] == cReg[cI])
                            {
                                val += field[nb];
                                count++;
                            }
                            else if (nb == cI && !popUpToDo[on] && cReg[on] == cReg[cI])
                            {
                                val += field[on];
                                count++;
                            }
                        }
                        else
                        {
                            const label patchI = this->boundaryMesh().whichPatch(fI);
                            if (this->boundary()[patchI].coupled())
                            {
                                const label bfI = fI - this->nInternalFaces();
                                if (!popUpToDoNbr[bfI] && cRegNbr[bfI] == cReg[cI])
                                {
                                    val += fieldNbr[bfI];
                                    count++;
                                }
                            }
                        }
                    }
                    if (field.name() == dbgField && cI == dbgCell && Pstream::myProcNo() == dbgProc)
                    {
                        Pout<< count << endl;
                    }
                    if (count)
                    {
                        val /= count;
                        field[cI] = val;
                        popUpToDoNext[cI] = false;
                        somethingChanged = true;
                    }
                    else
                    {
                        if (debugMode_)
                        {
                            Pout<< cI << tab
                                 << this->C()[cI] << tab
                                 << cReg[cI] << tab
                                 << this->cells()[cI] << tab
                                 << endl;
                         }

                        allDone = false;

                        // Patch fix for sudden pop cells (boundary or internal)
                        // without any neighbours to take information from
                        // This happens mostly in the solid region
                        // Now we do not care about the solution inside the
                        // solid region. This happens rarely
                        // Most probably wrong value will be assign.
                        // However mean value will be a "safe" value to avoid
                        // assigning a zero density in compressible cases
                        // Otherwise another information has to be imported from
                        // the user or derived with another method

                        field[cI] = val;
                        field[cI] = average(field.primitiveField());
                    }
                }
            }
            popUpToDo = popUpToDoNext;
            reduce(allDone, andOp<bool>());
            if (!allDone)
            {
                reduce(somethingChanged, orOp<bool>());
                if (!somethingChanged)
                {
                    WarningInFunction << "Unable to propagate initial values to all pop cells" << nl << endl;
                    break;
                }
                else if (debugMode_)
                {
                    Pout<< "Need another sweep for pop cell initialisation" << endl;
                }
            }
            else
            {
                break;
            }
        }

        // Have not synced boundary values yet
        List<Type> fieldNbr;
        syncTools::swapBoundaryCellList(*this, field.primitiveField(), fieldNbr);

        if (field.name() == dbgField && Pstream::myProcNo() == dbgProc)
        {
            Pout<< field[dbgCell] << endl;
        }

        scalarField oldVolg = oldVolgrow;

        if (true)
        {
            for (int iter=0; iter < popSubSteps_; iter++)
            {
                scalarField volStart = oldVolg;
                GeometricField<Type, PatchField, GeoMesh> fieldStart = field;
                // Midpoint rule
                for (label m = 0; m < 2; m++)
                {
                    scalar dd = 1.0/popSubSteps_;
                    if (m == 0)
                    {
                        dd /= 2;
                    }
                    Field<Type> df = Field<Type>
                    (
                        this->cells().size(), pTraits<Type>::zero
                    );
                    Type ff = pTraits<Type>::zero;
                    Field<Type> minff = Field<Type>
                    (
                        this->cells().size(), pTraits<Type>::one*GREAT
                    );
                    Field<Type> maxff = Field<Type>
                    (
                        this->cells().size(), -pTraits<Type>::one*GREAT
                    );

                    scalarField sV = scalarField(this->cells().size(), 0.0);

                    forAll(fPop, fI)
                    {
                        const label& fPopI = fPop[fI];

                        if (fPopI < this->nInternalFaces())
                        {
                            const label& on = this->owner()[fPopI];
                            const label& nb = this->neighbour()[fPopI];

                            ff = pTraits<Type>::zero;
                            scalar sVf = dd*phiPop[fPopI]*this->time().deltaTValue();

                            if (cRegion()[on] == cRegion()[nb])
                            {
                                if (popUpC[on]==1 && popUpC[nb]!=1)
                                {
                                    ff = field[nb];
                                }
                                else if (popUpC[nb]==1 && popUpC[on]!=1)
                                {
                                    ff = field[on];
                                }
                                else if (popUpC[nb]!=1 && popUpC[on]!=1)
                                {
                                    ff = 0.5*(field[on]+field[nb]);
                                }
                                else if (popUpC[nb]==1 && popUpC[on]==1)
                                {
                                    ff = 0.5*(field[on]+field[nb]);
                                }

                                if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && (on == dbgCell || nb == dbgCell))
                                {
                                    Pout<< on << " " << nb << " face: " << fPopI << " " << popUpC[on] << " " << popUpC[nb] << " " << cRegion()[on] << " " << cRegion()[nb] << " " << field[on] << " " << field[nb] << " " << ff << " " << ff*sVf << endl;
                                }

                                df[on] += ff*sVf;
                                df[nb] -= ff*sVf;

                                if (boundPopValues_)
                                {
                                    minff[on] = min(minff[on], ff);
                                    maxff[on] = max(maxff[on], ff);
                                    minff[nb] = min(minff[nb], ff);
                                    maxff[nb] = max(maxff[nb], ff);
                                }
                            }
                            else
                            {
                                // Don't try to conserve across the regions
                                // We should perhaps do something to ensure global
                                // conservation in each region, but allowing it
                                // to flow between them seems to make no sense.
                                // Setting it to zero here would do this but
                                // locally would probably produce spikes.
                                if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && (on == dbgCell || nb == dbgCell))
                                {
                                    Pout<< on << " " << nb << " face: " << fPopI << " " << popUpC[on] << " " << popUpC[nb] << " " << cRegion()[on] << " " << cRegion()[nb] << " " << field[on] << " " << field[nb] << " " << field[on]*sVf << " " << field[nb]*sVf << endl;
                                }

                                df[on] += field[on]*sVf;
                                df[nb] -= field[nb]*sVf;

                                if (boundPopValues_)
                                {
                                    minff[on] = min(minff[on], field[on]);
                                    maxff[on] = max(maxff[on], field[on]);
                                    minff[nb] = min(minff[nb], field[nb]);
                                    maxff[nb] = max(maxff[nb], field[nb]);
                                }
                            }

                            sV[on] += sVf;
                            sV[nb] -= sVf;
                        }
                        else
                        {
                            const label& on = this->faceOwner()[fPopI];
                            label patchI = this->boundaryMesh().whichPatch(fPopI);

                            if (!isA<indirectPolyPatch>(this->boundary()[patchI].patch()))
                            {
                                const label lpfI = fPopI - this->boundaryMesh()[patchI].start();
                                const label bfI = fPopI - this->nInternalFaces();
                                const scalar sVf = dd*
                                    (
                                        phiPop.boundaryField()[patchI][lpfI]
                                    )*
                                    this->time().deltaTValue();

                                if (this->boundary()[patchI].coupled())
                                {
                                    ff = pTraits<Type>::zero;

                                    if (cRegion()[on] == cRegNbr[bfI])
                                    {
                                        if (popUpC[on]==1 && popUpCNbr[bfI]!=1)
                                        {
                                            ff = fieldNbr[bfI];
                                        }
                                        else if (popUpCNbr[bfI]==1 && popUpC[on]!=1)
                                        {
                                            ff = field[on];
                                        }
                                        else if (popUpCNbr[bfI]!=1 && popUpC[on]!=1)
                                        {
                                            ff = 0.5*(field[on]+fieldNbr[bfI]);
                                        }
                                        else if (popUpCNbr[bfI]==1 && popUpC[on]==1)
                                        {
                                            ff = 0.5*(field[on]+fieldNbr[bfI]);
                                        }

                                        if (field.name() == dbgField && Pstream::myProcNo() == dbgProc && on == dbgCell)
                                        {
                                            Pout<< on << " " << " pface: " << fPopI << " " << popUpC[on] << " " << popUpCNbr[bfI] << " " << cRegion()[on] << " " << cRegNbr[bfI] << " " << field[on] << " " << fieldNbr[bfI] << " " << ff << " " << ff*sVf << endl;
                                        }

                                        df[on] += ff*sVf;
                                        sV[on] += sVf;

                                        if (boundPopValues_)
                                        {
                                            minff[on] = min(minff[on], ff);
                                            maxff[on] = max(maxff[on], ff);
                                        }
                                    }
                                }
                                else
                                {
                                    if (popUpC[on]==1)
                                    {
                                        ff = field[on];

                                        df[on] += ff*sVf;

                                        if (boundPopValues_)
                                        {
                                            minff[on] = min(minff[on], ff);
                                            maxff[on] = max(maxff[on], ff);
                                        }

                                        sV[on] += sVf;
                                    }
                                }
                            }
                        }
                    }
                    if (m == 1)
                    {
                        oldVolg = volStart;
                        field.ref() = fieldStart.ref();
                    }
                    forAll(popIndi, cI)
                    {
                        if (popIndi[cI])
                        {
                            scalar newVol = oldVolg[cI] + sV[cI];
                            scalar oldVol = oldVolg[cI];

                            if (mag(newVol) > SMALL)
                            {
                                field[cI] = (field[cI]*oldVol + df[cI])/newVol;
                            }
                            else
                            {
                                //- bad fix for now at proc
                                const labelList& cellFaces(this->cells()[cI]);
                                forAll(cellFaces, fI)
                                {
                                    label faceI = cellFaces[fI];
                                    if (faceI < this->nInternalFaces())
                                    {
                                        const label& on = this->owner()[faceI];
                                        const label& nb = this->neighbour()[faceI];
                                        if (popUpC[on]!=1)
                                        {
                                            field[cI] = field[on];
                                        }
                                        if (popUpC[nb]!=1)
                                        {
                                            field[cI] = field[nb];
                                        }
                                    }
                                }
                                //field[cI] = pTraits<Type>::zero;
                            }

                            if (boundPopValues_)
                            {
                                // Clamp the value so that it can't exceed the min/max
                                // of face values. Can occur if overall change in volume small
                                // compared to sum of abs values of face volume sweeps - then we
                                // end up with a large imbalance of incoming/outgoing fluxes due to slight
                                // variations in face values.
                                // TODO: redistribute conservatively rather
                                if (minff[cI] <= maxff[cI]) // Failsafe in case no faces had any influence
                                {
                                    field[cI] = max(minff[cI], field[cI]);
                                    field[cI] = min(maxff[cI], field[cI]);
                                }
                            }

                            if
                            (
                                field.name() == dbgField &&
                                Pstream::myProcNo() == dbgProc &&
                                cI == dbgCell
                            )
                            {
                                Pout<< oldVol << tab
                                     << "newVol: "
                                     << newVol << tab
                                     << field[cI] << tab
                                     << df[cI] << endl;
                            }

                            oldVolg[cI] = newVol;
                        }
                    }
                }
            }

            // Deal with ordinary boundaries attached to pop
            // cells: Their current value will not be meaningful
            // since they have changed regions; therefore
            // we initialise with extrapolated values
            forAll(fPop, fI)
            {
                const label& fPopI = fPop[fI];
                if (fPopI >= this->nInternalFaces())
                {
                    const label& on = this->faceOwner()[fPopI];
                    if (popIndi[on])
                    {
                        const label patchI = this->boundaryMesh().whichPatch(fPopI);
                        if (!isA<indirectPolyPatch>(this->boundary()[patchI].patch()))
                        {
                            const label lpfI = fPopI - this->boundaryMesh()[patchI].start();

                            // Copy and reassign the whole patch in order to
                            // respect any overriding of operator=
                            Field<Type> bf(field.boundaryField()[patchI]);
                            bf[lpfI] = field[on];
                            field.boundaryFieldRef()[patchI] = bf;
                        }
                    }
                }
            }
        }
        const labelList& fsc = fullSnapCells();
        // Have not synced final boundary values yet
        syncTools::swapBoundaryCellList
        (
            *this,
            field.primitiveField(),
            fieldNbr
        );
        forAll(fsc, cI)
        {
            label gcI = fsc[cI];
            if (popUpC[gcI] == 1)
            {
                field[gcI] = pTraits<Type>::zero;
                scalar counter = 0;
                forAll(this->cells()[gcI], fI)
                {
                    label faceI = this->cells()[gcI][fI];
                    if (faceI < this->nInternalFaces())
                    {
                        const label& on = this->owner()[faceI];
                        const label& nb = this->neighbour()[faceI];
                        if (cReg[nb] == cReg[on])
                        {
                            if (on == gcI)
                            {

                                if
                                (
                                    field.name() == dbgField &&
                                    Pstream::myProcNo() == dbgProc &&
                                    gcI == dbgCell
                                )
                                {
                                    Pout<< nb << " " << field[nb] << endl;
                                }

                                field[gcI] += field[nb];
                                counter += 1;
                            }
                            else
                            {

                                if
                                (
                                    field.name() == dbgField &&
                                    Pstream::myProcNo() == dbgProc &&
                                    gcI == dbgCell
                                )
                                {
                                    Pout<< on << " " << field[on] << endl;
                                }

                                field[gcI] += field[on];
                                counter += 1;
                            }
                        }
                    }
                    else
                    {
                        const label patchI =
                            this->boundaryMesh().whichPatch(faceI);
                        if (this->boundaryMesh()[patchI].coupled())
                        {
                            const label bfaceI = faceI-this->nInternalFaces();
                            if (cReg[gcI] == cRegNbr[bfaceI])
                            {
                                field[gcI] += fieldNbr[bfaceI];
                                counter++;
                            }
                        }
                    }
                }
                if (!counter)
                {
                    if (cRegion()[gcI] == 0)
                    {
                        FatalErrorInFunction
                            << "Found no usable neighbours for full snap pop cell: "
                            << endl
                            << "label: " << gcI << endl
                            << "C(): " << this->C()[gcI] << endl
                            << "CRefion: " << cRegion()[gcI] << endl
                            << "CRefion0: " << cRegion0()[gcI]
                            << exit(FatalError);
                    }
                    else
                    {
                        if (debugMode_)
                        {
                            Pout<< "Found no usable neighbours for full snap pop cell: "
                            << endl
                            << "However the cell is inside inactive zone" << endl
                            << "Initialization with average" << endl
                            << "label: " << gcI << endl
                            << "C(): " << this->C()[gcI] << endl
                            << "CRefion: " << cRegion()[gcI] << endl
                            << "CRefion0: " << cRegion0()[gcI]
                            << endl;
                        }
                        field[gcI] = average(field.primitiveField());
                        counter = 1;
                    }
                }

                field[gcI] /= counter;

                if
                (
                    (field.name() == dbgField) &&
                    (Pstream::myProcNo() == dbgProc) &&
                    (gcI == dbgCell)
                )
                {
                    Pout<< "fsc: " << field[gcI] << endl;
                }
            }
        }
        if
        (
            debugMode_ &&
            (
                (field.name() == "T") ||
                (field.name() == "rho")
            )
        )
        {
            forAll(field, i)
            {
                if (refCast<scalarField>(field)[i] < 0)
                {
                    Pout<< "Neg" << tab
                         << field.name() << tab
                         << "after pop2: "
                         << i << tab
                         << field[i]
                         << endl;
                }
            }
        }

        if (Pstream::myProcNo() == dbgProc && (!dbgField.size() || field.name() == dbgField))
        {
            Pout<< "After PopFields2: " << field.name() << "[" << dbgCell << "] = " << field[dbgCell] << endl;
        }
    }

}


template
<
    class Type
>
void Foam::dynamicGIBFvMesh::writeScalarField
(
    word name,
    const List<Type>& iF,
    bool writeNow
) const
{
    volScalarField outField
    (
        IOobject
        (
            name,
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("0", dimless, 0),
        "calculated"
    );

    //forAll loop because some = operations are
    //not general for instance boolList->scalarField
    //outField.internalField() = iF;

    forAll(outField, cI)
    {
        outField[cI] = iF[cI];
    }

    if (time().outputTime()||writeNow)
    {
        outField.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// end namespace FOAM
}
// ************************************************************************* //
