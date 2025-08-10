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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSurface/sampledSurface/sampledSurface.H"
#include "db/IOstreams/Fstreams/OFstream.H"

template<class Type>
bool Foam::sampledSurface::checkFieldSize(const Field<Type>& field) const
{
    if (faces().empty() || field.empty())
    {
        return false;
    }

    if (field.size() != faces().size())
    {
        FatalErrorInFunction
            << "size mismatch: "
            << "field (" << field.size()
            << ") != surface (" << faces().size() << ")"
            << exit(FatalError);
    }

    return true;
}


template<class Type>
Type Foam::sampledSurface::integrate(const Field<Type>& field) const
{
    Type value = Zero;

    if (checkFieldSize(field))
    {
        value = sum(field*magSf());
    }

    reduce(value, sumOp<Type>());
    return value;
}


template<class Type>
Type Foam::sampledSurface::integrate(const tmp<Field<Type>>& field) const
{
    Type value = integrate(field());
    field.clear();
    return value;
}


template<class Type>
Type Foam::sampledSurface::average(const Field<Type>& field) const
{
    Type value = Zero;

    if (checkFieldSize(field))
    {
        value = sum(field*magSf());
    }

    reduce(value, sumOp<Type>());

    // avoid divide-by-zero
    if (area())
    {
        return value/area();
    }
    else
    {
        return Zero;
    }
}


template<class Type>
Type Foam::sampledSurface::average(const tmp<Field<Type>>& field) const
{
    Type value = average(field());
    field.clear();
    return value;
}


template<class Type>
void Foam::sampledSurface::statistics
(
    const fileName& outputDir,
    const word& fieldName,
    const Field<Type>& field
) const
{
    checkFieldSize(field);

    scalar maxVal = -GREAT;
    scalar minVal = GREAT;
    scalar aveVal = 0;
    scalar aveNormal = 0;
    scalar stdDev = 0;
    scalar homogeneity = 0;
    scalar homogeneityN = 0;

    if (string(pTraits<Type>::typeName) == string("scalar"))
    {
        aveVal = sum(cmptAv(field)*magSf());
        maxVal = max(cmptAv(field));
        minVal = min(cmptAv(field));
        aveNormal = aveVal;
    }
    else
    {
        aveVal = sum(mag(field)*magSf());
        maxVal = max(mag(field));
        minVal = min(mag(field));
        if (string(pTraits<Type>::typeName) == string("vector"))
        {
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                aveNormal += sum(field.component(cmpt) * Sf().component(cmpt));
            }
        }
    }

    reduce(
        std::tie(aveVal, aveNormal, maxVal, minVal),
        ParallelOp<sumOp<scalar>, sumOp<scalar>, maxOp<scalar>, minOp<scalar>>{}
    );

    if (area())
    {
        aveVal /= area();
        aveNormal /= area();
    }
    else
    {
        aveVal = 0.;
    }

    if (string(pTraits<Type>::typeName) == string("scalar"))
    {
        stdDev = sum(sqr(cmptAv(field)-aveVal)*magSf());
        homogeneity = sum(mag(cmptAv(field)-aveVal)*magSf());
    }
    else
    {
        stdDev = sum(sqr(mag(field) - aveVal)*magSf());
        homogeneity = sum(mag(mag(field) - aveVal)*magSf());

        if (string(pTraits<Type>::typeName) == string("vector"))
        {
            scalarField vn(field.size(), 0.);
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                vn += field.component(cmpt) * Sf().component(cmpt);
            }
            homogeneityN = sum(mag((vn/(magSf()+SMALL))-aveNormal)*magSf());
        }
    }

    reduce(
        std::tie(stdDev, homogeneity, homogeneityN),
        UniformParallelOp<sumOp<scalar>, 3>{}
    );

    if (area())
    {
        stdDev /= area();
        homogeneity = ((mag(aveVal) > SMALL) ?
                       (1. - (homogeneity/(2.*mag(aveVal)*area()))) : 0.0);
        if (string(pTraits<Type>::typeName) == string("vector"))
        {
            homogeneityN = ((mag(aveNormal) > SMALL) ?
                           (1. - (homogeneityN/(2.*mag(aveNormal)*area()))) : 0.0);
        }
    }
    else
    {
        stdDev = 0.;
        homogeneity = 0.;
        homogeneityN = 0.;
    }

    stdDev = sqrt(stdDev);

    fileName logFileName = outputDir/name()+"_"+fieldName+"_surfaceStatistics.dat";
    OFstream logFile(logFileName);

    if (string(pTraits<Type>::typeName) == string("vector"))
    {
        logFile << "# Field "<<fieldName
                <<": area, min, max, ave, stdDev, homogeneity, homogeneity (surface normal)"<<endl;
        logFile << area() <<" "<<minVal <<" "<<maxVal<<" "<<aveVal<<" "<<stdDev
                <<" "<<homogeneity<<" "<<homogeneityN<<endl;

        Info<<"Surface : "  << name() <<" : Field : "<< fieldName
             <<" : area "<< area() <<" : min , "<< minVal <<" : max , "<< maxVal
             <<" : ave , "<< aveVal <<" : stdDev , "<< stdDev
             <<" : homogeneity , "<< homogeneity <<" : homogeneity (surface normal) , "
             << homogeneityN <<endl;
    }
    else
    {
        logFile << "# Field "<<fieldName
                <<": area, min, max, ave, stdDev, homogeneity"<<endl;
        logFile << area() <<" "<<minVal <<" "<<maxVal<<" "<<aveVal<<" "<<stdDev
                <<" "<<homogeneity<<endl;

        Info<<"Surface : "  << name() <<" : Field : "<< fieldName
             <<" : area "<< area() <<" : min , "<< minVal <<" : max , "<< maxVal
             <<" : ave , "<< aveVal <<" : stdDev , "<< stdDev
             <<" : homogeneity , "<< homogeneity
             <<endl;
    }

    return;
}


template<class Type>
void Foam::sampledSurface::statistics(const tmp<Field<Type>>& field) const
{
    statistics(field());
    field.clear();
    return;
}


template<class ReturnType, class Type>
void Foam::sampledSurface::project
(
    Field<ReturnType>& res,
    const Field<Type>& field
) const
{
    if (checkFieldSize(field))
    {
        const vectorField& norm = Sf();

        forAll(norm, facei)
        {
            res[facei] = field[facei] & (norm[facei]/mag(norm[facei]));
        }
    }
    else
    {
        res.clear();
    }
}


template<class ReturnType, class Type>
void Foam::sampledSurface::project
(
    Field<ReturnType>& res,
    const tmp<Field<Type>>& field
) const
{
    project(res, field());
    field.clear();
}


template<class ReturnType, class Type>
Foam::tmp<Foam::Field<ReturnType>>
Foam::sampledSurface::project
(
    const tmp<Field<Type>>& field
) const
{
    tmp<Field<ReturnType>> tRes(new Field<ReturnType>(faces().size()));
    project(tRes(), field);
    return tRes;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::sampledSurface::pointAverage
(
    const GeometricField<Type, pointPatchField, pointMesh>& pfld
) const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(pfld.mesh()());

    tmp<GeometricField<Type, fvPatchField, volMesh>> tcellAvg
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "cellAvg",
                mesh.time().timeName(),
                pfld.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>("zero", dimless, Zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& cellAvg = tcellAvg.ref();

    labelField nPointCells(mesh.nCells(), 0);
    {
        for (label pointi = 0; pointi < mesh.nPoints(); pointi++)
        {
            const labelList& pCells = mesh.pointCells(pointi);

            forAll(pCells, i)
            {
                label celli = pCells[i];

                cellAvg[celli] += pfld[pointi];
                nPointCells[celli]++;
            }
        }
    }
    forAll(cellAvg, celli)
    {
        cellAvg[celli] /= nPointCells[celli];
    }
    // Give value to calculatedFvPatchFields
    cellAvg.correctBoundaryConditions();

    return tcellAvg;
}


// ************************************************************************* //
