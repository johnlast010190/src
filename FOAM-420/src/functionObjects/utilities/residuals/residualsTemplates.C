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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "residuals/residuals.H"
#include "fields/volFields/volFields.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::residuals::writeFileBlockHeader
(
    Ostream& os,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (foundObject<fieldType>(fieldName))
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            const word resultName = fieldName + std::to_string(cmpt);
            writeDelimited
            (
                os,
                resultName
            );
        }
    }
}

template<class Type>
void Foam::functionObjects::residuals::writeFileHeader
(
    Ostream& os,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (foundObject<fieldType>(fieldName))
    {
        typename pTraits<Type>::labelType validComponents
        (
            mesh_.validComponents<Type>()
        );

        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            if (component(validComponents, cmpt) != -1)
            {
                const word resultName =
                    fieldName + word(pTraits<Type>::componentNames[cmpt]);
                writeDelimited
                (
                    os,
                    resultName
                );
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::residuals::initialiseField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    if (foundObject<volFieldType>(fieldName))
    {
        const Foam::dictionary& solverDict = mesh_.solverPerformanceDict();

        if (solverDict.found(fieldName))
        {
            typename pTraits<Type>::labelType validComponents
            (
                mesh_.validComponents<Type>()
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                if (component(validComponents, cmpt) != -1)
                {
                    const word resultName =
                        fieldName + word(pTraits<Type>::componentNames[cmpt]);

                    createField(resultName);
                }
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::residuals::initialiseBlockField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    if (foundObject<volFieldType>(fieldName))
    {
        const Foam::dictionary& solverDict = mesh_.blockSolverPerformanceDict();

        if (solverDict.found(fieldName))
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                const word resultName = fieldName + std::to_string(cmpt);

                createField(resultName);
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::residuals::writeResidual(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (foundObject<fieldType>(fieldName))
    {

        const Foam::dictionary& solverDict = mesh_.solverPerformanceDict();

        if (solverDict.found(fieldName))
        {
            const List<SolverPerformance<Type>> sp
            (
                solverDict.lookup(fieldName)
            );

            const Type& residual = sp.first().initialResidual();

            typename pTraits<Type>::labelType validComponents
            (
                mesh_.validComponents<Type>()
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                if (component(validComponents, cmpt) != -1)
                {
                    const word resultName =
                        fieldName + word(pTraits<Type>::componentNames[cmpt]);

                    const scalar r = component(residual, cmpt);

                    if (Pstream::master())
                    {
                        file() << token::TAB << r;
                    }
                    setResult(resultName, r);

                    writeField(resultName);
                }
            }
        }
    }
}

template<class Type>
void Foam::functionObjects::residuals::writeBlockResidual(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (foundObject<fieldType>(fieldName))
    {
        const Foam::dictionary& solverDict = mesh_.blockSolverPerformanceDict();

        if (solverDict.found(fieldName))
        {
            const List<BlockSolverPerformance<Type>> sp
            (
                solverDict.lookup(fieldName)
            );

            const Type& residual = sp.first().initialResidual();

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                const word resultName = fieldName +  std::to_string(cmpt);

                const scalar r = component(residual, cmpt);

                if (Pstream::master())
                {
                    file() << token::TAB << r;
                }
                setResult(resultName, r);

                writeField(resultName);
            }
        }
    }
}


// ************************************************************************* //
