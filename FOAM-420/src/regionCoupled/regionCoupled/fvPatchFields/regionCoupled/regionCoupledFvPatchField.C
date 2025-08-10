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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "fields/Fields/symmTransformField/symmTransformField.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::regionCoupledFvPatchField<Type>::lookupInactiveRegions()
{
    const objectRegistry& obr = this->db().time();
    if (obr.found("regionProperties"))
    {
        return
            !obr.lookupObject<IOdictionary>
            (
                "regionProperties"
            ).lookupOrDefault
            (
                "inactiveRegions",
                wordList()
            ).found(sampleRegion());
    }

    return true;
}


template<class Type>
Foam::fvPatchField<Type>*
Foam::regionCoupledFvPatchField<Type>::setHeFromT
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!regionCoupled())
    {
        NotImplemented;
        return nullptr;
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::regionCoupledFvPatchField<Type>::regionCoupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    regionCoupledPatch_(refCast<const mappedFvPatch>(p)),
    isRegionCoupled_(lookupInactiveRegions()),
    neighbourFieldName_(iF.name()),
    neighbourPhaseName_(iF.group()),
    neighbourRegion_(word::null),
    nonCoupledRegionBoundary_(setHeFromT(p, iF))
{}


template<class Type>
Foam::regionCoupledFvPatchField<Type>::regionCoupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    regionCoupledPatch_(refCast<const mappedFvPatch>(p)),
    isRegionCoupled_(lookupInactiveRegions()),
    neighbourFieldName_(iF.name()),
    neighbourPhaseName_(iF.group()),
    neighbourRegion_(word::null),
    nonCoupledRegionBoundary_(setHeFromT(p, iF))
{}


template<class Type>
Foam::regionCoupledFvPatchField<Type>::regionCoupledFvPatchField
(
    const regionCoupledFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    regionCoupledPatch_(refCast<const mappedFvPatch>(p)),
    isRegionCoupled_(lookupInactiveRegions()),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourPhaseName_(ptf.neighbourPhaseName_),
    neighbourRegion_(ptf.neighbourRegion_),
    nonCoupledRegionBoundary_
    (
        (!isRegionCoupled_ && ptf.nonCoupledRegionBoundary().rawPtr())
      ? fvPatchField<Type>::New
        (
            ptf.nonCoupledRegionBoundary().clone()(),
            p,
            iF,
            mapper
        ).ptr()
      : nullptr
    )
{}


template<class Type>
Foam::regionCoupledFvPatchField<Type>::regionCoupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    regionCoupledPatch_(refCast<const mappedFvPatch>(p)),
    isRegionCoupled_(lookupInactiveRegions()),
    neighbourFieldName_(dict.lookupOrDefault("neighbourFieldName", iF.name())),
    neighbourPhaseName_(dict.lookupOrDefault("neighbourPhaseName", iF.group())),
    neighbourRegion_(dict.lookupOrDefault("neighbourRegion", word::null))
{
    if (!regionCoupled())
    {
        if (dict.found("uncoupledBoundary"))
        {
            nonCoupledRegionBoundary_ =
                fvPatchField<Type>::New(p, iF, dict.subDict("uncoupledBoundary")).ptr();
        }
        else
        {
            nonCoupledRegionBoundary_ =
                fvPatchField<Type>::New("zeroGradient", p, iF).ptr();
        }
    }
    if (neighbourPhaseName_ == "none")
    {
        neighbourPhaseName_ = word::null;
    }
}


template<class Type>
Foam::regionCoupledFvPatchField<Type>::regionCoupledFvPatchField
(
    const regionCoupledFvPatchField<Type>& ptf
)
:
    coupledFvPatchField<Type>(ptf),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    isRegionCoupled_(lookupInactiveRegions()),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourPhaseName_(ptf.neighbourPhaseName_),
    nonCoupledRegionBoundary_
    (
        (!isRegionCoupled_ && ptf.nonCoupledRegionBoundary().rawPtr())
      ? ptf.nonCoupledRegionBoundary()().clone().ptr()
      : nullptr
    )
{}


template<class Type>
Foam::regionCoupledFvPatchField<Type>::regionCoupledFvPatchField
(
    const regionCoupledFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    isRegionCoupled_(lookupInactiveRegions()),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourPhaseName_(ptf.neighbourPhaseName_),
    nonCoupledRegionBoundary_
    (
        (!isRegionCoupled_ && ptf.nonCoupledRegionBoundary().rawPtr())
      ? ptf.nonCoupledRegionBoundary()().clone().ptr()
      : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::regionCoupledFvPatchField<Type>::coupledField() const
{
    if (neighbourFieldName_ != this->internalField().name())
    {
        return neighbourFieldName_;
    }
    else if (neighbourPhaseName_ != this->internalField().group())
    {
        return
            IOobject::groupName
            (
                this->internalField().member(),
                neighbourPhaseName_
            );
    }
    else
    {
        return neighbourFieldName_;
    }

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    if (!regionCoupled())
    {
        return this->patchInternalField();
    }

    return
        interpolateFromNeighbour
        (
            neighbourFvPatchField().patchInternalField(),
            this->patchInternalField()
        );
}


template<class Type>
const Foam::objectRegistry&
Foam::regionCoupledFvPatchField<Type>::neighbourDb() const
{
    if (!regionCoupled())
    {
        FatalErrorInFunction
            << "Neighbour can't be called for deactivated regions."
            << exit(FatalError);
    }

    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    const fvMesh& nbrMesh = nbrPatch.boundaryMesh().mesh();

    // Select the one identified by neighbourRegion, if supplied
    if (neighbourRegion_ != word::null)
    {
        const fvSolutionRegistry* nbrReg =
            nbrMesh.lookupObjectPtr<fvSolutionRegistry>(neighbourRegion_);
        return nbrReg ? nbrReg->registry() : nbrMesh;
    }

    // Otherwise, if there is more than one solution registry in neighbour mesh,
    // check for one that has a region coupled patch, if neighbourRegion not
    // specified
    HashTable<const fvSolutionRegistry*> nbrRegistries =
        nbrMesh.lookupClass<fvSolutionRegistry>();
    if (!nbrRegistries.size())
    {
        return nbrMesh;
    }
    if (nbrRegistries.size() == 1)
    {
        return nbrRegistries.begin()()->registry();
    }
    const objectRegistry* foundRegistry(nullptr);
    forAllConstIters(nbrRegistries, nbrRegIter)
    {
        const GeometricField<Type, fvPatchField, volMesh>* nbrFieldPtr =
            nbrRegIter()->registry().lookupObjectPtr
            <
                GeometricField<Type, fvPatchField, volMesh>
            >(this->coupledField());
        if (nbrFieldPtr)
        {
            const fvPatchField<Type>& pf =
                nbrPatch.patchField
                <
                    GeometricField<Type, fvPatchField, volMesh>, Type
                >(*nbrFieldPtr);
            if (isA<regionCoupledFvPatchField>(pf))
            {
                if (foundRegistry)
                {
                    FatalErrorInFunction
                        << "Found multiple solution regions coupled to patch "
                        << this->patch().name() << " in field "
                        << this->internalField().name() << ", region "
                        << this->internalField().db().parent().name() << "/"
                        << this->internalField().db().name() << "." << nl;
                    FatalErrorInFunction
                        << "Please use the 'neighbourRegion' keyword to "
                        << "specify the correct neighbour region."
                    << exit(FatalError);
                }
                foundRegistry = &nbrRegIter()->registry();
            }
        }
    }
    if (!foundRegistry)
    {
        // Don't give an error, because this could be called before the neigh-
        // bouring region has been constructed
        return nbrMesh;
    }
    else
    {
        return *foundRegistry;
    }
}


template<class Type>
const Foam::regionCoupledFvPatchField<Type>&
Foam::regionCoupledFvPatchField<Type>::neighbourFvPatchField() const
{
    if (!regionCoupled())
    {
        FatalErrorInFunction
            << "Neighbour can't be called for deactivated regions."
            << exit(FatalError);
    }

    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    return
        refCast<const regionCoupledFvPatchField<Type>>
        (
            nbrPatch.template lookupPatchFieldInDb
            <
                GeometricField<Type, fvPatchField, volMesh>,
                Type
            >
            (
                neighbourDb(),
                this->coupledField()
            )
        );
}


template<class Type>
void Foam::regionCoupledFvPatchField<Type>::evaluate(const Pstream::commsTypes commsType)
{
    if (!regionCoupled())
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }

        nonCoupledRegionBoundary_().evaluate(commsType);
        fvPatchField<Type>::evaluate(commsType);
        return;
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().snGrad(deltaCoeffs);
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::snGrad() const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().snGrad();
    }

    // Coupled version: should not be called
    NotImplemented;
}



template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().gradientInternalCoeffs(deltaCoeffs);
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::gradientInternalCoeffs() const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().gradientInternalCoeffs();
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().gradientBoundaryCoeffs();
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().gradientBoundaryCoeffs(deltaCoeffs);
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& test
) const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().valueInternalCoeffs(test);
    }

    // Coupled version: should not be called
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& test
) const
{
    if (!regionCoupled())
    {
        return nonCoupledRegionBoundary_().valueBoundaryCoeffs(test);
    }

    // Coupled version: should not be called
    NotImplemented;
}



template<class Type>
void Foam::regionCoupledFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    if (!regionCoupled())
    {
        return;
    }

    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    const labelUList& nbrFaceCells = nbrPatch.faceCells();
    // Default to symmetry in case of holes in AMI mapping
    // (low weight)
    scalarField pnf
    (
        interpolateFromNeighbour
        (
            scalarField(psiInternal, nbrFaceCells),
            symmetryValue
            (
                this->patchInternalField()().component(cmpt), this->patch().Sf()
            ),
            cmpt
        )
    );

    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::regionCoupledFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    if (!regionCoupled())
    {
        return;
    }

    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    const labelUList& nbrFaceCells = nbrPatch.faceCells();
    // Default to symmetry in case of holes in AMI mapping
    // (low weight)
    Field<Type> pnf
    (
        interpolateFromNeighbour
        (
            Field<Type>(psiInternal, nbrFaceCells),
            symmetryValue
            (
                this->patchInternalField(), this->patch().Sf()
            )
        )
    );

    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchField<Type>::interpolateFromNeighbour
(
    const scalarField& pf,
    const UList<scalar>& defaultValues,
    const direction cmpt
) const
{
    if (!regionCoupled())
    {
        FatalErrorInFunction
            << "Neighbour can't be called for deactivated regions."
            << exit(FatalError);
    }

    // Do the AMI interp
    return regionCoupledPatch_.interpolate(pf, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::interpolateFromNeighbour
(
    const Field<Type>& pf, const UList<Type>& defaultValues
) const
{
    if (!regionCoupled())
    {
        FatalErrorInFunction
            << "Neighbour can't be called for deactivated regions."
            << exit(FatalError);
    }

    // Do the AMI interp
    return regionCoupledPatch_.interpolate(pf, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::regionCoupledFvPatchField<Type>::interpolateToNeighbour
(
    const Field<Type>& pf, const UList<Type>& defaultValues
) const
{
    if (!regionCoupled())
    {
        FatalErrorInFunction
            << "Neighbour can't be called for deactivated regions."
            << exit(FatalError);
    }

    return neighbourFvPatchField().interpolateFromNeighbour(pf, defaultValues);
}


template<class Type>
void Foam::regionCoupledFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (neighbourFieldName_ != this->internalField().name())
    {
        os.writeEntry("neighbourFieldName", neighbourFieldName_);
    }
    if (neighbourPhaseName_ != this->internalField().group())
    {
        os.writeEntry
        (
            "neighbourPhaseName",
            (neighbourPhaseName_ == word::null ? "none" : neighbourPhaseName_)
        );
    }
    if (nonCoupledRegionBoundary_)
    {
        os.beginBlock("uncoupledBoundary");
        nonCoupledRegionBoundary_().write(os);
        os.endBlock();
    }
}


// ************************************************************************* //
