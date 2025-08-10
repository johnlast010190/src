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
    (c) 2016 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "EquationFlux.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::fv::EquationFlux<Type>::writeFileHeader
(
    Ostream& os,
    const fvMatrix<Type>& eqn
)
const
{
    // Add headers to output data
    writeFile::writeHeader(os, "'" + equation_  + "' equation flux density");
    ios_base::fmtflags oldflags =
        os.setf(ios_base::fmtflags(0), ios_base::floatfield);
    writeHeaderValue(os, "Dimensions - min/max", eqn.dimensions()/dimArea);
    writeHeaderValue(os, "Dimensions - integrated", eqn.dimensions());
    writeCommented(os, "Time");
    for (auto& pf : eqn.psi().boundaryField())
    {
        if (!pf.coupled() && !isA<emptyPolyPatch>(pf.patch().patch()))
        {
            word patchName = pf.patch().name();
            writeDelimited(os, "min("+patchName+")");
            writeDelimited(os, "max("+patchName+")");
            writeDelimited(os, "integrated("+patchName+")");
        }
    }
    os  << endl;
    os.setf(oldflags, ios_base::floatfield);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::EquationFlux<Type>::EquationFlux
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(name, modelType, dict, obr),
    writeFile(obr, name, typeName, dict),
    writeField_(false),
    storeField_(false),
    headerWritten(false)
{
    this->read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::EquationFlux<Type>::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        writeFile::read(dict);

        equation_ = word(dict.lookup("equation"));

        outputFieldName_ =
            dict.lookupOrDefault
            (
                "outputFieldName",
                word(equation_+"EquationFlux")
            );
        writeField_ =
            dict.lookupOrDefault("writeField", dict.found("outputFieldName"));
        storeField_ =
            dict.lookupOrDefault("storeField", false);

        fieldNames_.setSize(1, equation_);
        applied_.setSize(1, false);

        mesh().schemes().setFluxRequired(equation_);

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::fv::EquationFlux<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInformation
        << "EquationFlux<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    tmp<GeometricField<Type, fvPatchField, volMesh>> tflux;
    if
    (
        !storeField_
     || !obr_.foundObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            outputFieldName_
        )
    )
    {
        tflux =
            tmp<GeometricField<Type, fvPatchField, volMesh>>
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    outputFieldName_,
                    eqn.psi().db(),
                    mesh_,
                    eqn.dimensions()/dimArea,
                    pTraits<Type>::zero
                )
            );
    }
    if (storeField_)
    {
        tflux.ptr()->store();
    }
    GeometricField<Type, fvPatchField, volMesh>& flux =
    (
        tflux.valid()
      ? tflux.ref()
      : obr_.lookupObjectRef<GeometricField<Type, fvPatchField, volMesh>>
        (
            outputFieldName_
        )
    );

    if (Pstream::master())
    {
        if (!headerWritten)
        {
            writeFileHeader(file(), eqn);
        }
        writeTime(file());
    }
    headerWritten = true;

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> teqnFlux = eqn.flux();

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();
    forAll(psi.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pf = psi.boundaryField()[patchi];
        if (!pf.coupled() && !isA<emptyPolyPatch>(pf.patch().patch()))
        {
            Field<Type>& pFlux = flux.boundaryFieldRef()[patchi];
            pFlux = teqnFlux->boundaryField()[patchi]/pf.patch().magSf();

            Type minfp = gMin(pFlux);
            Type maxfp = gMax(pFlux);
            Type integralfp = gSum(pFlux*pf.patch().magSf());

            if (Pstream::master())
            {
                file()
                    << token::TAB << minfp
                    << token::TAB << maxfp
                    << token::TAB << integralfp;
            }
        }
    }
    if (Pstream::master())
    {
        file() << endl;
    }

    if (writeField_ && mesh_.time().outputTime())
    {
        flux.write();
    }
}


// ************************************************************************* //
