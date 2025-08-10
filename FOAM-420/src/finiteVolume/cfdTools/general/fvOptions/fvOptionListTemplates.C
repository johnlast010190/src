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
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif
#include "algorithms/serialThreads/serialThreads.H"
#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const bool block
)
{
    return this->operator()(field, field.name(), block);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    serialThreads::pauseSwitching();

    checkApplied();

    checkCorrect("operator", field.name());

    const dimensionSet ds = field.dimensions()/dimTime*dimVolume;

    tmp<fvMatrix<Type>> tmtx(new fvMatrix<Type>(field, ds));
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(regionOptions_, regioni)
    {
        forAll(regionOptions_[regioni], i)
        {
            option& source = regionOptions_[regioni].operator[](i);

            label fieldi =
                source.applyToField(fieldName, thisRegionName());

            if (fieldi != -1)
            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fvopt, "fvOption()." + source.name());
                #endif

                source.setApplied(fieldi);

                if (source.isActive())
                {
                    if (debug)
                    {
                        Info<< "Applying source " << source.name()
                            << " to field " << fieldName << endl;
                    }

                    source.addSup(mtx, fieldi);
                }
            }
        }
    }

    serialThreads::resumeSwitching();

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName,
    const bool block
)
{
    checkApplied();

    checkCorrect("operator", field.name());

    const dimensionSet ds = field.dimensions()/dimTime*dimVolume;

    tmp<fvBlockMatrix<Type>> tmtx(new fvBlockMatrix<Type>(field, ds));

    fvBlockMatrix<Type>& mtx = tmtx.ref();

    forAll(regionOptions_, regioni)
    {
        forAll(regionOptions_[regioni], i)
        {
            option& source = regionOptions_[regioni].operator[](i);

            label fieldi =
                source.applyToField(fieldName, thisRegionName());

            if (fieldi != -1)
            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fvopt, "fvOption()." + source.name());
                #endif

                source.setApplied(fieldi);

                if (source.isActive())
                {
                    if (debug)
                    {
                        Info<< "Applying source " << source.name() << " to field "
                            << fieldName << endl;
                    }

                    source.addSup(mtx, fieldi);
                }
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const bool block
)
{
    return this->operator()(rho, field, field.name(), block);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    serialThreads::pauseSwitching();

    checkApplied();

    checkCorrect("operator", field.name());

    const dimensionSet ds
    (
        rho.dimensions()*field.dimensions()/dimTime*dimVolume
    );

    tmp<fvMatrix<Type>> tmtx(new fvMatrix<Type>(field, ds));
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(regionOptions_, regioni)
    {
        forAll(regionOptions_[regioni], i)
        {
            option& source = regionOptions_[regioni].operator[](i);

            label fieldi =
                source.applyToField(fieldName, thisRegionName());

            if (fieldi != -1)
            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fvopt, "fvOption()." + source.name());
                #endif

                source.setApplied(fieldi);

                if (source.isActive())
                {
                    if (debug)
                    {
                        Info<< "Applying source " << source.name() << " to field "
                            << fieldName << endl;
                    }

                    source.addSup(rho, mtx, fieldi);
                }
            }
        }
    }

    serialThreads::resumeSwitching();

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvBlockMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName,
    const bool block
)
{
    checkApplied();

    checkCorrect("operator", field.name());

    const dimensionSet ds
    (
        rho.dimensions()*field.dimensions()/dimTime*dimVolume
    );

    tmp<fvBlockMatrix<Type>> tmtx(new fvBlockMatrix<Type>(field, ds));

    fvBlockMatrix<Type>& mtx = tmtx.ref();

    forAll(regionOptions_, regioni)
    {
        forAll(regionOptions_[regioni], i)
        {
            option& source = regionOptions_[regioni].operator[](i);

            label fieldi =
                source.applyToField(fieldName, thisRegionName());

            if (fieldi != -1)
            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fvopt, "fvOption()." + source.name());
                #endif

                source.setApplied(fieldi);

                if (source.isActive())
                {
                    if (debug)
                    {
                        Info<< "Applying source " << source.name() << " to field "
                            << fieldName << endl;
                    }

                    source.addSup(rho, mtx, fieldi);
                }
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(alpha, rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    serialThreads::pauseSwitching();

    checkApplied();

    checkCorrect("operator", field.name());

    const dimensionSet ds
    (
        alpha.dimensions()*rho.dimensions()*field.dimensions()
       /dimTime*dimVolume
    );

    tmp<fvMatrix<Type>> tmtx(new fvMatrix<Type>(field, ds));
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(regionOptions_, regioni)
    {
        forAll(regionOptions_[regioni], i)
        {
            option& source = regionOptions_[regioni].operator[](i);

            label fieldi =
                source.applyToField(fieldName, thisRegionName());

            if (fieldi != -1)
            {
                #if !defined( WIN32 ) && !defined( WIN64 )
                addProfiling(fvopt, "fvOption()." + source.name());
                #endif

                source.setApplied(fieldi);

                if (source.isActive())
                {
                    if (debug)
                    {
                        Info<< "Applying source " << source.name() << " to field "
                            << fieldName << endl;
                    }
                    source.addSup(alpha, rho, mtx, fieldi);
                }
            }
        }
    }

    serialThreads::resumeSwitching();

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    volScalarField one
    (
        IOobject
        (
            "one",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar("one", dimless, 1.0)
    );

    return this->operator()(alpha, one, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fv::optionList::boundarySources
(
    const word& fieldName,
    const label patchID,
    const Field<Type>& pf,
    Field<Type>& bsourceDerivs
)
{
    serialThreads::pauseSwitching();

    checkBoundaryApplied();

    checkCorrect("boundary", fieldName);

    const word patchName = this->mesh_.boundary()[patchID].name();
    tmp<Field<Type>> bsources
    (
        new Field<Type>
        (
            this->mesh_.boundary()[patchID].size(),
            pTraits<Type>::zero
        )
    );
    bsourceDerivs =
        Field<Type>
        (
            this->mesh_.boundary()[patchID].size(),
            pTraits<Type>::zero
        );

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldPatchi =
            source.applyToBoundaryFieldAndPatch(fieldName, patchID);

        if (fieldPatchi != -1)
        {
            source.setBoundaryApplied(fieldName, fieldPatchi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying boundary source " << source.name()
                        << " to field " << fieldName
                        << ", patch " << patchName << endl;
                }
                source.addBoundarySource
                (
                    fieldName,
                    patchID,
                    pf,
                    bsources.ref(),
                    bsourceDerivs
                );
            }
        }
    }

    serialThreads::resumeSwitching();

    return bsources;
}


template<class Type>
void Foam::fv::optionList::constrain(fvMatrix<Type>& eqn)
{
    serialThreads::pauseSwitching();

    checkApplied();

    checkCorrect("constrain", eqn.psi().name());

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(eqn.psi().name(), thisRegionName());

        if (fieldi != -1)
        {
            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(fvopt, "fvOption::constrain." + eqn.psi().name());
            #endif

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying constraint " << source.name()
                        << " to field " << eqn.psi().name() << endl;
                }

                source.constrain(eqn, fieldi);
            }
        }
    }

    serialThreads::resumeSwitching();
}


template<class Type>
void Foam::fv::optionList::constrain(fvBlockMatrix<Type>& eqn)
{
    serialThreads::pauseSwitching();

    checkApplied();

    checkCorrect("constrain", eqn.psi().name());

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(eqn.psi().name(), thisRegionName());

        if (fieldi != -1)
        {
            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(fvopt, "fvOption::constrain." + eqn.psi().name());
            #endif

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying constraint " << source.name()
                        << " to field " << eqn.psi().name() << endl;
                }

                source.constrain(eqn, fieldi);
            }
        }
    }

    serialThreads::resumeSwitching();
}


template<class Type>
void Foam::fv::optionList::correct
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    serialThreads::pauseSwitching();

    const word& fieldName = field.name();

    checkCorrect("correct", fieldName);

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(fieldName, thisRegionName());

        if (fieldi != -1)
        {
            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(fvopt, "fvOption::correct." + source.name());
            #endif

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Correcting source " << source.name()
                        << " for field " << fieldName << endl;
                }

                source.correct(field);
            }
        }
    }

    serialThreads::resumeSwitching();
}



// ************************************************************************* //
