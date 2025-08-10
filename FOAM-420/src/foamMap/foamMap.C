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
    (c) 2019 Esi Ltd.


Class
    Foam::foamMap

Group
    grpfoamMap

Description
    A set of classes designed to map an arbitrary number of flow fields from an unstructured
    mesh region into a uniform grid, or another unstructured mesh region.

SourceFile
    foamMap.C

\*---------------------------------------------------------------------------*/

#include "foamMap.H"
#include "db/IOobjectList/IOobjectList.H"

using namespace Foam;

//template<class Type, template<class> class PatchField, class GeoMesh>
//const word GeometricField<Type, PatchField, GeoMesh>::typeName;

foamMap::foamMap()
{
    wdistMap_   = false;
    norm_       = false;
    parmodel_   = 2; //memory-efficient mode
    nwdist_     = 100;
    extcoef_    = 0.02;
    rhoRefSource_ = 1.205;
    rhoRefTarget_ = 1.205;
    myid_       = 0;
    nprocs_     = 1;
    masterid_   = 0;
    deltaBox_   = 0.1;
    interp_     = false;
    alphaMax_   = 0.5;
    errorBound_ = 1.0e-7;
    mapTimeName_= "latestTime";
    tgtTimeName_= "latestTime";
    sourceCase_ = word::null;
    caseType_   = "vehicle";
    mapBoundary_= true;
    mapFixedBC_ = false;
    function_   = "mapping"; //or validate
    UrotDegreeFromSource_ = 0; //rotating degree frpm source to target, ani-clockwise positive
    bboxType_   = "smBoundBox"; //or boundBox
    boxCellNum_ = 40000;
    sourceRegions_ = wordList(1, fvMesh::defaultRegion);
    targetRegions_ = wordList(1, fvMesh::defaultRegion);
    targetRegionDirs_ = wordList(1, word::null);
    allRegions_ = false;
    reportMem_  = false;
    scaleSource_= false;
}



scalar foamMap::map_error
(
    const word& name
)
{
    const word unknown("unknown");
    const word& type = fieldTypes_.lookup(name, unknown);

    if (type == "unknown")
    {
        Info<< "error validate field mapping for " << name << " unknown field." << endl;
        return 0;
    }

    scalar perror = 0;
    scalar small = 1.0e-25;
    const fvMesh* mesh = target().mesh;
    scalar coef = target().scale(name);
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);
        scalar alphai = target().wdists_[i];
        std::vector<label> nb;
        std::vector<scalar>dist;
        source().search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            perror += 1.0; //totally missed
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        scalar value = source().scalarfields_[name][ic];
        scalar pi = value*coef;
        scalar prefi = target().getScalarField(name, i);

        scalar perrori = fabs(pi - prefi)/(fabs(prefi) + small);

        target().getScalarField(name).primitiveFieldRef()[i] = pi;
        perror += perrori;

    }
    Info<< "Total relative mapping error for scalar field: " << name << "  is "
        << perror << endl;
    return perror;
}


scalar foamMap::Umap_error
(
    const word& name
)
{
    const word unknown("unknown");
    const word& type = fieldTypes_.lookup(name, unknown);

    if (type == "unknown")
    {
        Info<< "error validate field mapping for " << name
            << " unknown field." << endl;
        return 0;
    }

    scalar perror = 0;
    scalar small = 1.0e-25;
    const fvMesh* mesh = target().mesh;
    scalar coef = target().scale(name);
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);
        scalar alphai = target().wdists_[i];
        std::vector<label> nb;
        std::vector<scalar>dist;
        source().search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            perror += 1.0; //totally missed
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        point value = source().vectorfields_[name][ic];
        point pi = value*coef;
        point prefi = target().getVectorField(name, i);

        scalar perrori = mag(pi-prefi)/(mag(prefi)+small);
        target().getVectorField(name).primitiveFieldRef()[i] = pi;
        perror += perrori;

    }
    Info<< "Total relative mapping error for vector field: " << name
        << " is " << perror << endl;
    return perror;
}

label foamMap::search_error()
{

    //seach error from taget field
    label cnt = 0;
    const fvMesh* mesh = target().mesh;
    label  ncells = mesh->C().size();
    List<label> ref(ncells);

    for (label i = 0; i < ncells; i++)
    {
        ref[i] = i;
    }

    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();

        x = xbar(x);
        y = ybar(y);
        z = zbar(z);
        scalar alphai = target().wdists_[i];

        std::vector<label> nb;
        std::vector<scalar> dist;

        source().search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            cnt++;
            continue;
        }
        label k = nb[0];
        label ic = source().cellIndex(alphai, k);


        if (ic != ref[i])
        {
            cnt++;
        }
    }//forAll

    Info<< "number of mis-matched cells:" << cnt << endl;
    return cnt;
}

void foamMap::stringtok
(
    std::vector<word>& container,
    const word& in0,
    const word& del
)
{
    word in = in0;
    if (del == ",")
    {
        replace_str(in);
    }

    label len = in.length();
    label i = 0;

    if (container.size() >= 1)
    {
        container.resize(0);
    }
    while (i < len)
    {
        i = in.find_first_not_of(del, i);
        if (i < 0)
        {
            return;
        }
        label  j = in.find_first_of(del, i);

        if (j < 0)
        {
            word ss = in.substr(i);
            trim(ss, " ");
            container.push_back (ss);
            break;
        }
        else
        {
            word ss1 = in.substr(i, j-i);
            trim(ss1, " ");
            container.push_back(ss1);
        }

        i = j + 1;
    }

    for (unsigned int k = 0; k < container.size(); k++)
    {
        word sk = container[k];
        trim(sk, " ");
        container[k] = sk;
    }
    return;
}

label foamMap::trimleft
(
    word& str,
    const word& del
)
{
    label i0 = str.find_first_not_of(del);
    if (i0 > 0)
    {
        label sz = str.size() - i0;
        str = str.substr(i0, sz);
    }

    return i0;
}

label foamMap::trimright
(
    word& str,
    const word& del
)
{
    label i0 = str.find_last_not_of(del);
    label sz = str.size() - 1;
    if (i0 > 0 && i0 < sz)
    {
        str = str.substr(0, i0 + 1);
    }

    return i0;
}

void foamMap::trim(word& str, const word& del)
{
    trimleft(str, del);
    trimright(str, del);
}

void foamMap::replace_str
(
    word& s1
)
{
    word nstr;
    word s2;
    char c0 = ' ';
    //remove \t,\r,\n
    bool tostart = false;
    label i2 = s1.find_last_not_of(' ');

    for (label i = 0; i <= i2; i++)
    {
        char c = s1[i];
        if (!tostart)
        {
            if (c != ' ' && c != '\r' && c != '\t' && c != '\n')
            {
                tostart = true;
                if (c0 == ',' && c == ',')
                {
                    s2.push_back('-');
                }
                c0 = c;
                s2.push_back(c);
            }
        }
        else
        {
            if (c == '\t' || c == '\r' || c == '\n')
            {
                continue;
            }

            if (c == ' ' && c0 == ' ')
            {
                continue;
            }

            if (c0 == ',' && c == ',')
            {
                s2.push_back('-');
            }
            s2.push_back(c);
            c0 = c;
        }
    }
    s1 = s2;
    return;
}

void foamMap::mapFields()
{
    //scalar fields
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();

        if (parRun())
        {
            Info<< "mapping scalar par... " << scname << endl;
            mapScalarFieldPar(scname);
        }
        else
        {
            mapScalarField(scname);
        }
    }

    wait();

    //vector field
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& vecname = iter();

        if (parRun())
        {
            Info<< "mapping vector par... " << vecname << endl;
            mapVectorFieldPar(vecname);
        }
        else
        {
            mapVectorField(vecname);
        }
    }
}

void foamMap::mapScalarFieldPar
(
    const word& scname
)
{
    const fvMesh* mesh = target().mesh;

    scalar coef = target().scale(scname);
    //internal field map
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();

        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai =
            wdistMap_
          ? target().wdists_[i]
          : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;

        search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            continue;
        }

        label k = nb[0];

        if (k < 0 || k >= sourceScalarFields_[scname].size())
        {
            continue;
        }

        label ic = cellIndex(alphai, k); //properly considered alpha>alphaMax situation

        if (ic < 0 || ic >= sourceScalarFields_[scname].size())
        {
            continue;
        }
        scalar value = sourceScalarFields_[scname][ic];

        scalar pi = value*coef;
        target().getScalarField(scname).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        //boundary field
        Info<< "mapping boundary field..." << endl;
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if (type == "empty" || type == "symmetryPlane" || type == "symmetry")
            {
                continue;
            }

            if
            (
                !mapFixedBC_
             && isA<fixedValueFvPatchScalarField>
                (
                    target().getScalarField(scname).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << scname
                    << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getScalarField(scname).boundaryField()[r].type()
                    << endl;
                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                search_nearest
                (
                    nb,
                    dist,
                    x,
                    y,
                    z
                );

                if (nb.empty())
                {
                    continue;
                }
                label k = nb[0];
                label ic = cellMap_[k];

                scalar pi = sourceScalarFields_[scname][ic];
                pi *= coef;
                target().getScalarField(scname).boundaryFieldRef()[r][j] = pi;
            }//forAll j
        }//forAll r
    }

    Info<< "Scalar field " << scname << " mapped." << endl;
}

void foamMap::getFieldTypes
(
    const fvMesh& mesh,
    const word& timeName
)
{
    IOobjectList objects(mesh, timeName);

    IOobjectList fields  = objects.lookupClass(volScalarField::typeName);
    IOobjectList vfields = objects.lookupClass(volVectorField::typeName);

    mapScalarFields_ = HashSet<word>();
    mapVectorFields_ = HashSet<word>();

    if (mapFieldNames_.size() == 0)
    {
        forAllIter(IOobjectList, fields, fieldIter)
        {
            mapScalarFields_.insert(fieldIter()->name());
        }
        forAllIter(IOobjectList, vfields, fieldIter)
        {
            mapVectorFields_.insert(fieldIter()->name());
        }
    }
    else
    {
        forAll(mapFieldNames_, i)
        {
            const word& name = mapFieldNames_[i];
            if (fields.lookup(name) != nullptr)
            {
                mapScalarFields_.insert(name);
            }
            else if (vfields.lookup(name) != nullptr)
            {
                mapVectorFields_.insert(name);
            }
        }
    }
}

void foamMap::scalePressure
(
    const dimensionSet& sourceDim,
    const dimensionSet& targetDim
)
{
    // Unset parallel flag to avoid problems with rho field missing on some procesors
    const bool oldParRun = Pstream::parRun();
    Pstream::parRun() = false;

    if
    (
        sourceDim == dimensionSet(1, -1, -2, 0, 0, 0, 0)
     && targetDim == dimensionSet(0,  2, -2, 0, 0, 0, 0)
    )
    {
        Info<< "Dividing source pressure by rho field..." << endl;

        volScalarField rhoField
        (
            IOobject
            (
                "rho",
                source().runTime_->timeName(),
                *(source().mesh),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            *(source().mesh),
            dimensionedScalar
            (
                "one",
                dimMass/dimVolume,
                1
            )
        );

        Info<< "rho dimensions: " << rhoField.dimensions() << endl;

        for (label i = 0; i < rhoField.primitiveField().size(); i++)
        {
            source().scalarfields_["p"][i] /= rhoField.primitiveField()[i];
        }

        //boundaries
        if (mapBoundary_)
        {
            label cnt = rhoField.primitiveField().size();

            forAll(source().mesh->boundaryMesh(), r)
            {
                word type = source().mesh->boundaryMesh()[r].type();
                if (type == "empty" || type == "symmetryPlane")
                {
                    continue;
                }

                forAll(source().mesh->boundaryMesh()[r], j)
                {
                    scalar sc = rhoField.boundaryField()[r][j];
                    source().scalarfields_["p"][cnt] /= sc;
                    cnt++;
                }
            }
        }
    }
    else if
    (
        targetDim == dimensionSet(1, -1, -2, 0, 0, 0, 0)
     && sourceDim == dimensionSet(0,  2, -2, 0, 0, 0, 0)
    )
    {
        Info<< "Multipling source pressure by rho (rhoRefSource = "
            << rhoRefSource_ << ")" << endl;

        for (label i = 0; i < source().scalarfields_["p"].size(); i++)
        {
            source().scalarfields_["p"][i] *= rhoRefSource_;
        }
    }

    Pstream::parRun() = oldParRun;
    Info<< endl;
}

void foamMap::checkFieldDimensions()
{
    Info<< "Checking field dimensions..." << endl;

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const volScalarField& sourceField = source().getScalarField(name);
        const dimensionSet& sourceDim = sourceField.dimensions();

        const volScalarField& targetField = target().getScalarField(name);
        const dimensionSet& targetDim = targetField.dimensions();

        if (sourceDim != targetDim)
        {
            if (name == "p")
            {
                Info<< "source field " << name << " dimensions: " << sourceDim << endl;
                Info<< "target field " << name << " dimensions: " << targetDim << endl;
                scalePressure(sourceDim, targetDim);
            }
            else
            {
                FatalErrorInFunction
                    << "Field " << name << " dimensions do not match." << nl
                    << "source dimensions: " << sourceDim << nl
                    << "target dimensions: " << targetDim << nl
                    << exit(FatalError);
            }
        }
    }

    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const volVectorField& sourceField = source().getVectorField(name);
        const dimensionSet& sourceDim = sourceField.dimensions();

        const volVectorField& targetField = target().getVectorField(name);
        const dimensionSet& targetDim = targetField.dimensions();

        if (sourceDim != targetDim)
        {
            FatalErrorInFunction
                << "Field " << name << " dimensions do not match." << nl
                << "source dimensions: " << sourceDim << nl
                << "target dimensions: " << targetDim << nl
                << exit(FatalError);
        }
    }
}


void foamMap::mapScalarField
(
    const word& scname
)
{
    const fvMesh* mesh = target().mesh;
    scalar coef = target().scale(scname);
    //internal field map
    label cnt = 0;
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai =
            wdistMap_
          ? target().wdists_[i]
          : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;
        nb.resize(2);
        dist.resize(2);
        source().search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            cnt += 1;
            continue;
        }

        label k = nb[0];

        label ic = source().cellIndex(alphai, k);

        if (ic != i)
        {
            cnt++;
        }

        scalar value = source().scalarfields_[scname][ic];
        scalar pi = value*coef;
        target().getScalarField(scname).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        //boundary field
        Info<< "mapping boundary field..." << endl;
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if (type == "empty" || type == "symmetryPlane" || type == "symmetry")
            {
                continue;
            }

            if
            (
                !mapFixedBC_
             && isA<fixedValueFvPatchScalarField>
                (
                    target().getScalarField(scname).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << scname
                    << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getScalarField(scname).boundaryField()[r].type()
                    << endl;
                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                source().search_nearest_bface
                (
                    nb,
                    dist,
                    x,
                    y,
                    z
                );

                if (nb.empty())
                {
                    continue;
                }
                label k = nb[0];
                label ic = source().bfaceIndex(k);

                scalar pi = source().scalarfields_[scname][ic];
                pi *= coef;
                target().getScalarField(scname).boundaryFieldRef()[r][j] = pi;
            }//forAll j
        }//forAll r
    }

    Info<< "Scalar field " << scname << " mapped." << " " << cnt << endl;
}

void foamMap::mapVectorFieldPar
(
    const word& vecname
)
{
    const fvMesh* mesh = target().mesh;
    scalar coef = target().scale(vecname);
    //internal field map
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai =
            wdistMap_
          ? target().wdists_[i]
          : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;
        search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            continue;
        }

        label k = nb[0];

        if (k < 0 || k >= sourceVectorFields_[vecname].size())
        {
            continue;
        }

        label ic = cellIndex(alphai, k);
        point value = sourceVectorFields_[vecname][ic];

        point pi = value*coef;
        target().getVectorField(vecname).primitiveFieldRef()[i] = pi;
    }


    if (mapBoundary_)
    {
        //boundary field
        Info<< "mapping boundary field..." << endl;
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if (type == "empty" || type == "symmetryPlane" || type == "symmetry")
            {
                continue;
            }

            if
            (
                !mapFixedBC_
             && isA<fixedValueFvPatchVectorField>
                (
                    target().getVectorField(vecname).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << vecname
                    << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getVectorField(vecname).boundaryField()[r].type()
                    << endl;
                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                search_nearest
                (
                    nb,
                    dist,
                    x,
                    y,
                    z
                );

                if (nb.empty())
                {
                    continue;
                }

                label k = nb[0];
                label ic = cellMap_[k];
                point pi = sourceVectorFields_[vecname][ic];
                pi *= coef;
                target().getVectorField(vecname).boundaryFieldRef()[r][j] = pi;
            }//forAll j
        }//forAll r
    }

    wait();
    Info<< "vector field " << vecname << " mapped." << endl;
    return;
}

void foamMap::mapVectorField
(
    const word& vecname
)
{
    const fvMesh* mesh = target().mesh;
    scalar coef = target().scale(vecname);
    scalar beta = UrotDegreeFromSource_*3.1415926/180.0;

    //scalar uycoef=0;
    //internal field map
    forAll(mesh->C(), i)
    {
        scalar x = mesh->C()[i].x();
        scalar y = mesh->C()[i].y();
        scalar z = mesh->C()[i].z();
        x = xbar(x);
        y = ybar(y);
        z = zbar(z);

        scalar alphai =
            wdistMap_
          ? target().wdists_[i]
          : GREAT;

        std::vector<label> nb;
        std::vector<scalar> dist;

        source().search_nearest
        (
            nb,
            dist,
            alphai,
            x,
            y,
            z
        );

        if (nb.empty())
        {
            continue;
        }

        label k = nb[0];
        label ic = source().cellIndex(alphai, k);

        point value = source().vectorfields_[vecname][ic];
        scalar pmag = mag(value);
        if (fabs(beta) > 1.0e-4 && pmag > 1.0e-3)
        {
            scalar alpha0 = asin(value.y()/pmag);
            scalar vx = pmag*cos(alpha0 - beta);
            scalar vy = pmag*sin(alpha0 - beta);
            scalar vz = value.z();
            value = point(vx, vy, vz);
        }
        point pi = value*coef;

        target().getVectorField(vecname).primitiveFieldRef()[i] = pi;
    }

    if (mapBoundary_)
    {
        //boundary field
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if (type == "empty" || type == "symmetryPlane" || type == "symmetry")
            {
                continue;
            }

            if
            (
                !mapFixedBC_
             && isA<fixedValueFvPatchVectorField>
                (
                    target().getVectorField(vecname).boundaryField()[r]
                )
            )
            {
                Info<< "Not mapping " << vecname
                    << " on patch "
                    << mesh->boundaryMesh()[r].name()
                    << " with boundary condition "
                    << target().getVectorField(vecname).boundaryField()[r].type()
                    << endl;
                continue;
            }

            forAll(mesh->boundaryMesh()[r], j)
            {
                point ptf = mesh->Cf().boundaryField()[r][j];
                scalar x = ptf.x();
                scalar y = ptf.y();
                scalar z = ptf.z();
                x = xbar(x);
                y = ybar(y);
                z = zbar(z);

                std::vector<label> nb;
                std::vector<scalar> dist;

                source().search_nearest_bface
                (
                    nb,
                    dist,
                    x,
                    y,
                    z
                );


                if (nb.empty())
                {
                    continue;
                }
                label k = nb[0];
                label ic = source().bfaceIndex(k);

                point pi = source().vectorfields_[vecname][ic];
                scalar pmag = mag(pi);
                if (fabs(beta) > 1.0e-4 && pmag > 1.0e-3)
                {
                    scalar alpha0 = asin(pi.y()/pmag);
                    scalar vx = pmag*cos(alpha0 - beta);
                    scalar vy = pmag*sin(alpha0 - beta);
                    scalar vz = pi.z();
                    pi = point(vx, vy, vz);
                }
                pi *= coef;
                target().getVectorField(vecname).boundaryFieldRef()[r][j] = pi;
            }//forAll j
        }//forAll r
    }
    Info<< "vector field " << vecname << " mapped." << endl;
}

void foamMap::setInput
(
    const dictionary& dict
)
{
    if
    (
        sourceCase_ == word::null // foamSample initializes sourceCase to "./"
    )
    {
        sourceCase_ = fileName(dict.lookup("sourceCase"));
    }

    if (dict.found("fields"))
    {
        mapFieldNames_ = wordList(dict.lookup("fields"));
    }

    if (dict.found("mapScalarFields"))
    {
        mapScalarFields_ = HashSet<word>(dict.lookup("mapScalarFields"));
    }

    if (dict.found("mapVectorFields"))
    {
        mapVectorFields_ = HashSet<word>(dict.lookup("mapVectorFields"));
    }

    nwdist_ = dict.lookupOrDefault<label>("nwdist", 80);

    rhoRefSource_ = dict.lookupOrDefault<scalar>("rhoRefSource", 1.205);

    rhoRefTarget_ = dict.lookupOrDefault<scalar>("rhoRefTarget", 1.205);

    vector uref(30, 0, 0);
    UrefSource_ = dict.lookupOrDefault<vector>("UrefSource", uref);
    UrefTarget_ = dict.lookupOrDefault<vector>("UrefTarget", uref);

    interp_ = dict.lookupOrDefault<bool>("interpolation", false);

    UrotDegreeFromSource_ = dict.lookupOrDefault<scalar>("UrotDegreeFromSource", 0);

    alphaMax_ = dict.lookupOrDefault<scalar>("alphaMax", 0.5);

    errorBound_ = dict.lookupOrDefault<scalar>("errorBound", 1.0e-7);

    if (dict.found("sourceTime"))
    {
        ITstream sourceTimeEntry = dict.lookup("sourceTime");

        if (sourceTimeEntry[0] == word("latestTime"))
        {
            mapTimeName_ = word("latestTime");
        }
        else
        {
            scalar sourceTimeScalar = readScalar(dict.lookup("sourceTime"));
            mapTimeName_ = name(sourceTimeScalar);
        }
    }

    if (dict.found("targetTime"))
    {
        ITstream targetTimeEntry = dict.lookup("targetTime");

        if (targetTimeEntry[0] == word("latestTime"))
        {
            tgtTimeName_ = word("latestTime");
        }
        else
        {
            scalar targetTimeScalar = readScalar(dict.lookup("targetTime"));
            tgtTimeName_ = name(targetTimeScalar);
        }
    }

    if (dict.found("sourceRegion"))
    {
        sourceRegions_ = wordList(1, word(dict.lookup("sourceRegion")));
    }

    if (dict.found("targetRegion"))
    {
        targetRegions_ = wordList(1, word(dict.lookup("targetRegion")));
    }

    if (dict.found("allRegions"))
    {
        allRegions_ = readBool(dict.lookup("allRegions"));
    }

    if (dict.found("regionMaps"))
    {
        List<Pair<word>> regionMap;
        regionMap = List<Pair<word>>(dict.lookup("regionMaps"));

        sourceRegions_.resize(0);
        targetRegions_.resize(0);

        forAll(regionMap, regioni)
        {
            sourceRegions_.append(regionMap[regioni].first());
            targetRegions_.append(regionMap[regioni].second());
        }

        targetRegionDirs_ = targetRegions_;
    }

    mapBoundary_ = dict.lookupOrDefault<bool>("mapBoundary", true);
    mapFixedBC_ = dict.lookupOrDefault<bool>("mapFixedBC", false);

    caseType_ = dict.lookupOrDefault<word>("caseType", "vehicle");
    wdistMap_ = dict.lookupOrDefault<bool>("wdist_map", false);
    extcoef_ = dict.lookupOrDefault<scalar>("extensionCoef", 0.1);
    parmodel_ = dict.lookupOrDefault<label>("parallelModel", 2);
    function_ = dict.lookupOrDefault<word>("function", "mapping");
    norm_ = dict.lookupOrDefault<bool>("coordsNormalise", false);

    //memory usage info
    reportMem_ = dict.lookupOrDefault<bool>("reportMemoryUsage", false);
    //debug
    bool debug = dict.lookupOrDefault<bool>("debug", false);
    if (debug)
    {
        reportMem_ = true;
    }

    if (dict.found("fieldTypes"))
    {
        fieldTypes_ = HashTable<word>(dict.lookup("fieldTypes"));
    }

    //wall patch ids for calculating wall distances
    if (dict.found("patches"))
    {
        wallPatchNames_ = wordReList(dict.lookup("patches"));
    }

    // smart boundBox for HPC run
    bboxType_ = dict.lookupOrDefault<word>("boundBoxType", "smBoundBox"); //boundBox
    boxCellNum_ = dict.lookupOrDefault<label>("gridCells", 1000000);
    //Info<<"bound box type:"<<bboxType_<<endl;

    scaleSource_ = dict.lookupOrDefault<bool>("scaleSource", false);

    return;
}

void foamMap::setInput
(
    const dictionary& dict,
    const argList& args
)
{
    if
    (
        sourceCase_ == word::null // foamSample initializes sourceCase to "./"
     && !args.optionFound("sourceCase")
    )
    {
        sourceCase_ = fileName(dict.lookup("sourceCase"));
    }

    if (dict.found("fields"))
    {
        mapFieldNames_ = wordList(dict.lookup("fields"));
    }

    if (dict.found("mapScalarFields"))
    {
        mapScalarFields_ = HashSet<word>(dict.lookup("mapScalarFields"));
    }

    if (dict.found("mapVectorFields"))
    {
        mapVectorFields_ = HashSet<word>(dict.lookup("mapVectorFields"));
    }

    nwdist_ = dict.lookupOrDefault<label>("nwdist", 80);

    rhoRefSource_ = dict.lookupOrDefault<scalar>("rhoRefSource", 1.205);

    rhoRefTarget_ = dict.lookupOrDefault<scalar>("rhoRefTarget", 1.205);

    vector uref(30, 0, 0);
    UrefSource_ = dict.lookupOrDefault<vector>("UrefSource", uref);
    UrefTarget_ = dict.lookupOrDefault<vector>("UrefTarget", uref);

    interp_ = dict.lookupOrDefault<bool>("interpolation", false);

    UrotDegreeFromSource_ = dict.lookupOrDefault<scalar>("UrotDegreeFromSource", 0);

    alphaMax_ = dict.lookupOrDefault<scalar>("alphaMax", 0.5);

    errorBound_ = dict.lookupOrDefault<scalar>("errorBound", 1.0e-7);

    if (dict.found("sourceTime"))
    {
        ITstream sourceTimeEntry = dict.lookup("sourceTime");

        if (sourceTimeEntry[0] == word("latestTime"))
        {
            mapTimeName_ = word("latestTime");
        }
        else
        {
            scalar sourceTimeScalar = readScalar(dict.lookup("sourceTime"));
            mapTimeName_ = name(sourceTimeScalar);
        }
    }

    if (dict.found("targetTime"))
    {
        ITstream targetTimeEntry = dict.lookup("targetTime");

        if (targetTimeEntry[0] == word("latestTime"))
        {
            tgtTimeName_ = word("latestTime");
        }
        else
        {
            scalar targetTimeScalar = readScalar(dict.lookup("targetTime"));
            tgtTimeName_ = name(targetTimeScalar);
        }
    }

    if (dict.found("sourceRegion"))
    {
        sourceRegions_ = wordList(1, word(dict.lookup("sourceRegion")));
    }

    if (dict.found("targetRegion"))
    {
        targetRegions_ = wordList(1, word(dict.lookup("targetRegion")));
    }

    if (dict.found("allRegions"))
    {
        allRegions_ = readBool(dict.lookup("allRegions"));
    }

    if (dict.found("regionMaps"))
    {
        List<Pair<word>> regionMap;
        regionMap = List<Pair<word>>(dict.lookup("regionMaps"));

        sourceRegions_.resize(0);
        targetRegions_.resize(0);

        forAll(regionMap, regioni)
        {
            sourceRegions_.append(regionMap[regioni].first());
            targetRegions_.append(regionMap[regioni].second());
        }

        targetRegionDirs_ = targetRegions_;
    }

    mapBoundary_ = dict.lookupOrDefault<bool>("mapBoundary", true);
    mapFixedBC_ = dict.lookupOrDefault<bool>("mapFixedBC", false);

    caseType_ = dict.lookupOrDefault<word>("caseType", "vehicle");
    wdistMap_ = dict.lookupOrDefault<bool>("wdist_map", false);
    extcoef_ = dict.lookupOrDefault<scalar>("extensionCoef", 0.1);
    parmodel_ = dict.lookupOrDefault<label>("parallelModel", 2);
    function_ = dict.lookupOrDefault<word>("function", "mapping");
    norm_ = dict.lookupOrDefault<bool>("coordsNormalise", false);

    //memory usage info
    reportMem_ = dict.lookupOrDefault<bool>("reportMemoryUsage", false);
    //debug
    bool debug = dict.lookupOrDefault<bool>("debug", false);
    if (debug)
    {
        reportMem_ = true;
    }

    if (dict.found("fieldTypes"))
    {
        fieldTypes_ = HashTable<word>(dict.lookup("fieldTypes"));
    }

    //wall patch ids for calculating wall distances
    if (dict.found("patches"))
    {
        wallPatchNames_ = wordReList(dict.lookup("patches"));
    }

    // smart boundBox for HPC run
    bboxType_ = dict.lookupOrDefault<word>("boundBoxType", "smBoundBox"); //boundBox
    boxCellNum_ = dict.lookupOrDefault<label>("gridCells", 1000000);
    //Info<<"bound box type:"<<bboxType_<<endl;

    scaleSource_ = dict.lookupOrDefault<bool> ("scaleSource", false);

    return;
}

void foamMap::getWallDistPatchs
(
    const fvMesh& mesh
)
{
    Info<< "\nCalculating wall distance for ";

    //wall patch ids for calculating wall distances
    if (wallPatchNames_.size()==0)
    {
        Info<< "all patches." <<endl;
        return;
    }

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    wallDistPatchs_ = bMesh.patchSet(wallPatchNames_);

    const wordList allPatchNames(bMesh.names());

    Info<< wallDistPatchs_.size() << " patch"
        << (wallDistPatchs_.size() > 1 ? "es:" : ":")
        << endl;

    forAllConstIter(labelHashSet, wallDistPatchs_, it)
    {
        Info<< token::TAB << allPatchNames[*it] <<endl;
    }

    Info<< nl;

    return;
}


void foamMap::setOptions
(
    const argList& args
)
{
    args.optionReadIfPresent("sourceCase", sourceCase_);

    if (args.optionFound("mapScalarFields"))
    {
        args.optionLookup("mapScalarFields")() >> mapScalarFields_;
    }

    if (args.optionFound("mapVectorFields"))
    {
        args.optionLookup("mapVectorFields")() >> mapVectorFields_;
    }

    if (args.optionFound("mapSurfaceScalarFields"))
    {
        args.optionLookup("mapSurfaceScalarFields")() >> mapSurfaceScalarFields_;
    }

    args.optionReadIfPresent("alphaMax", alphaMax_);

    args.optionReadIfPresent("errorBound", errorBound_);

    if (args.optionFound("mapBoundary"))
    {
        mapBoundary_ = true;
    }

    if (args.optionFound("mapFixedBC"))
    {
        mapFixedBC_ = true;
    }

    args.optionReadIfPresent("sourceTime", mapTimeName_);

    args.optionReadIfPresent("targetTime", tgtTimeName_);

    if (args.optionFound("sourceRegion"))
    {
        sourceRegions_ = wordList(1, args["sourceRegion"]);
    }

    if (args.optionFound("targetRegion"))
    {
        targetRegions_ = wordList(1, args["targetRegion"]);
        targetRegionDirs_ = targetRegions_;
    }

    if (args.optionFound("allRegions"))
    {
        allRegions_ = true;
    }

    if (args.optionFound("regionMaps"))
    {
        sourceRegions_.resize(0);
        targetRegions_.resize(0);

        List<Pair<word>> regionMap;
        args.optionLookup("regionMaps")() >> regionMap;

        forAll(regionMap, regioni)
        {
            sourceRegions_.append(regionMap[regioni].first());
            targetRegions_.append(regionMap[regioni].second());
        }

        targetRegionDirs_ = targetRegions_;
    }

    args.optionReadIfPresent("function", function_);

    label nndist = 0;
    args.optionReadIfPresent("nwdist", nndist);
    if (nndist > 2 && nndist <= 1001)
    {
        nwdist_ = nndist;
    }

    args.optionReadIfPresent("rhoRefSource", rhoRefSource_);

    args.optionReadIfPresent("rhoRefTarget", rhoRefTarget_);

    args.optionReadIfPresent("UrefSource", UrefSource_);

    args.optionReadIfPresent("UrefTarget", UrefTarget_);

    args.optionReadIfPresent("UrotDegreeFromSource", UrotDegreeFromSource_);

    if (args.optionFound("interpolation"))
    {
        interp_ = true;
    }

    args.optionReadIfPresent("fieldTypes", fieldTypes_);

    return;
}

void foamMap::setMapTime
(
    Time& runTime
)
{
    scalar itime = atof(runTime.timeName().c_str());
    scalar maptime = atof(mapTimeName_.c_str());
    runTime += (maptime - itime);
}

void foamMap::buildpKdTrees()
{
    Info<< "building parallel search trees..." << endl;

    sourceAlphaMax_ = alphaMax_;

    label treeNum = nwdist_ + 1;
    deltaY_ = 1.0/float(nwdist_);

    kdTrees_.resize(treeNum);
    alphaMap_.resize(treeNum);
    const point& pmin = target().gbox().min();
    const point& pmax = target().gbox().max();
    xmin_ = pmin.x();
    xmax_ = pmax.x();
    ymin_ = pmin.y();
    ymax_ = pmax.y();
    zmin_ = pmin.z();
    zmax_ = pmax.z();

    Array2d<point> pxyz; //put the cell into groups
    pxyz.resize(treeNum);

    label n2 = sourceWdists_.size();

    forAll(sourceXyz_, i)
    {
        const point& pt = sourceXyz_[i];
        if (i > n2 - 1 || sourceWdists_[i] > sourceAlphaMax_)
        {
            continue;
        }

        label ic = inGroup(i);
        if (ic < 0) ic = 0;
        if (ic > treeNum - 1) ic = treeNum - 1;

        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());
        point pt1 (xb, yb, zb);
        pxyz[ic].push_back(pt1);
        alphaMap_[ic].push_back(i);
    }

    label inactiveKdTrees = 0;
    if (wdistMap_)
    {
        //use relative coordinates to build kdTrees
        for (unsigned ic = 0; ic < pxyz.size(); ic++)
        {
            label npt = pxyz[ic].size();
            if (npt < 1)
            {
                kdTrees_[ic].deactivate();
                inactiveKdTrees++;
                continue;
            }
            farray2d trainData(npt);
            for (label i = 0; i < npt; i++)
            {
                trainData[i].resize(3);
                point pt = pxyz[ic][i];
                trainData[i][0] = pt.x();
                trainData[i][1] = pt.y();
                trainData[i][2] = pt.z();
            }

            if (!interp_)
            {
                kdTrees_[ic].k = 1;
            }

            kdTrees_[ic].eps = errorBound_;

            kdTrees_[ic].build(trainData);
        }
    }

    if
    (
        wdistMap_
     && inactiveKdTrees < label(kdTrees_.size())
    )
    {
        scalar coef = 0.8;
        construct_internal_knn(coef*alphaMax_);
    }
    else
    {
        wdistMap_ = false;

        if (!interp_)
        {
            internalMap_.k = 1;
        }
        construct_internal_knn(-0.1);
    }

    wait();
    Info<< "build parallel search trees completed." << endl;
    return;
}


void foamMap::createSource
(
    const fvMesh* mesh,
    const Time* runTime
)
{
    source_ =
        new readFields
        (
            mesh,
            runTime,
            this,
            "source"
        );

    source().setNorm(norm_);
    source().constructBoundBox(0.0, bboxType_, boxCellNum_);
}


void foamMap::createTarget
(
    const fvMesh* mesh,
    const Time* runTime
)
{
    target_ =
        new readFields
        (
            mesh,
            runTime,
            this,
            "target"
        );

    target().setNorm(norm_);
    target().constructBoundBox(extcoef_, bboxType_, boxCellNum_);
}

void foamMap::clearSourceTarget()
{
    internalMap_.deAllocate();

    if (wdistMap_)
    {
        for (size_t i = 0; i < kdTrees_.size(); i++)
        {
            if (kdTrees_[i].active)
            {
                kdTrees_[i].deAllocate();
            }
        }
    }

    sourceXyz_.clear();
    sourceScalarFields_.clear();
    sourceVectorFields_.clear();
    sourceWdists_.clear();

    source_.clear();
    target_.clear();
}

void foamMap::construct_internal_knn(scalar cutoff)
{

    label cnt = 0;

    if (!wdistMap_) cutoff = -0.1;
    farray2d trainData;
    //for gridMap (cutoff<0), no need to do alpha-mapping, internal map contains all cells
    label ncells = sourceXyz_.size();
    if (wdistMap_) ncells = min(sourceXyz_.size(), sourceWdists_.size());

    for (label i = 0; i < ncells; i++)
    {
        if (wdistMap_)
        {
            scalar alphai = sourceWdists_[i];
            if (alphai < cutoff)
            {
                continue;
            }
        }

        point pt = sourceXyz_[i];
        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());
        std::vector<scalar> pts(3);
        pts[0] = xb;
        pts[1] = yb;
        pts[2] = zb;
        trainData.push_back(pts);
        cellMap_[cnt] = i;
        cnt++;
    }

    internalMap_.eps = errorBound_;

    internalMap_.build(trainData);

    return;
}


void foamMap::parInit(label np)
{
    //initialization for parallel run
    if (!parRun())
    {
        return;
    }

    procs_.resize(np);
    for (label i = 0; i < np; i++)
    {
        procs_[i] = i;
    }

}

bool foamMap::inlist
(
    label k,
    const DynamicList<label>& list
)
{
    if (list.size() == 0)
    {
        return false;
    }
    for (label i = 0; i < list.size(); i++)
    {
        if (list[i] == k)
        {
            return true;
        }
    }
    return false;
}

bool foamMap::inlist
(
    const word& key,
    const std::vector<word>& list
)
{
    if (list.size() == 0)
    {
        return false;
    }
    for (unsigned i = 0; i < list.size(); i++)
    {
        if (list[i] == key)
        {
            return true;
        }
    }
    return false;
}


template<class T>
void Foam::interProcTrans
(
    DynamicList<T>& field,
    const DynamicList<point>& verts,
    const List<boundBox>& tboxs,
    const List<boundBox>& sboxs,
    const DynamicList<T>& source,
    const word& name
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    label myid = Pstream::myProcNo();
    for (label ib = 0; ib < Pstream::nProcs(); ib++)
    {
        const boundBox& tgtbox = tboxs[ib];
        const boundBox& mySourceBox = sboxs[myid];
        //skip if mySourceBox does not overlap target box
        DynamicList<T> cellField;
        bool boxInDomain = mySourceBox.overlaps(tgtbox);
        if (boxInDomain)
        {
            for (label i = 0; i < source.size(); i++)
            {
                const point& pt = verts[i];

                if (tgtbox.contains(pt))
                {
                    cellField.append(source[i]);
                }
            }
        }
        UOPstream toProc(ib, pBufs);
        toProc << cellField;
    }
    pBufs.finishedSends();

    //Info<<"get data lists..."<<endl;

    field.resize(0);
    for (label ib = 0; ib < Pstream::nProcs(); ib++)
    {
        DynamicList<T> posGlob;
        UIPstream fromProc(ib, pBufs);
        fromProc >> posGlob ;
        for (label i = 0; i < posGlob.size(); i++)
        {
            field.append(posGlob[i]);
        }
    }

    pBufs.finishedSends();


    Info<< "field transferred for: " << name
        << ", size:" << field.size() << endl;
}


void foamMap::reportMemory
(
    memInfo& mem,
    const word& desc,
    scalar& maxUsage
)
{
    scalar memory = scalar(mem.update().size())*1.0e-3;
    if (memory > maxUsage)
    {
        scalar maxUsage = memory;

        Pout<< "Maximum memory usage now: " << maxUsage << " MB," << desc << endl;
    }
}


void foamMap::getDomainFields()
{
    if (parModel() == 1)
    {
        return getDomainFields0();
    }

    const smBoundBox& tbox = target().bbox();

    List<boundBox> tboxs(Pstream::nProcs());
    label myid = Pstream::myProcNo();
    tboxs[myid] = tbox;

    Pstream::allGatherList(tboxs); //make it available for other processors

    const smBoundBox& sbox = source().bbox();

    List<boundBox> sboxs(Pstream::nProcs());
    sboxs[myid] = sbox;
    Pstream::allGatherList(sboxs); //make it available for other processors
    // cell centre coordinates

    scalar maxUsage = 0;
    word desc = "getDomainFields: before get sourceXYZ";
    memInfo mem;

    if (reportMemUsage())
    {
        reportMemory(mem, desc, maxUsage);
    }

    DynamicList<point> sourceXyz;

    interProcTrans
    (
        sourceXyz,
        source().xyz(),
        tboxs,
        sboxs,
        source().xyz(),
        "cellVerts"
    );

    List<bool> status(sourceXyz.size(), false);
    sourceXyz_.resize(0);
    label cellCount = 0;
    forAll(sourceXyz, i)
    {
        const point& pt = sourceXyz[i];
        if (tbox.pointInBox(pt))
        {
            status[i] = true;
            sourceXyz_.append(pt);
            cellCount++;
        }
    }

    if (reportMemUsage())
    {
        Pout<< "Number of cells in this processor: " << source().xyz().size() << endl;
        Pout<< "Points from all processors: " << sourceXyz.size() << endl;
    }
    sourceXyz.clear();


    if (reportMemUsage())
    {
        desc = "getDomainFields: after get sourceXYZ";
        reportMemory(mem, desc, maxUsage);
        Pout<< "sourceXyz_:" << sourceXyz_.size() << endl;
        Pout<< "Number of active sourceXyz_ count:" << cellCount << endl;
    }

    wait();

    if (wdistMap_)
    {
        DynamicList<scalar> sourceWdists;
        //wall distance
        interProcTrans
        (
            sourceWdists,
            source().xyz(),
            tboxs,
            sboxs,
            source().wdists_,
            "wallDistance"
        );

        forAll(sourceWdists,i)
        {
            if (i < status.size() && status[i])
            {
                sourceWdists_.append(sourceWdists[i]);
            }
        }

        sourceWdists.clear();

        if (reportMemUsage())
        {
            desc = "getDomainFields: after get wallDists";
            reportMemory(mem, desc, maxUsage);
        }

        wait();
    }

    //scalar fields
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const DynamicList<scalar>& scfld=source().scalarfields_[name];
        DynamicList<scalar> sourceScalarFields;
        interProcTrans
        (
            sourceScalarFields,
            source().xyz(),
            tboxs,
            sboxs,
            scfld,
            name
        );

        sourceScalarFields_[name].resize(0);
        forAll(sourceScalarFields, i)
        {
            if (i < status.size() && status[i])
            {
                sourceScalarFields_[name].append(sourceScalarFields[i]);
            }
        }

        sourceScalarFields.clear();


        if (reportMemUsage())
        {
            desc = "getDomainFields: after get scalar field: ";
            desc += name;
            reportMemory(mem, desc, maxUsage);
        }
    }

    wait();
    //vector fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();
        const DynamicList<vector>& vecfld = source().vectorfields_[name];
        DynamicList<vector> sourceVectorFields;
        interProcTrans
        (
            sourceVectorFields,
            source().xyz(),
            tboxs,
            sboxs,
            vecfld,
            name
        );

        sourceVectorFields_[name].resize(0);
        forAll(sourceVectorFields, i)
        {
            if (i < status.size() && status[i])
            {
                sourceVectorFields_[name].append(sourceVectorFields[i]);
            }
        }
        sourceVectorFields.clear();
    }

    if (reportMemUsage())
    {
        desc = "getDomainFields: after get vector fields";
        reportMemory(mem, desc, maxUsage);
    }

    wait();

    return;
}


void foamMap::getDomainFields0()
{
    if (!Pstream::parRun())
    {
        return;
    }

    Info<< "Single-block parallel model..." << endl;

    List<vectorField> posLoc(Pstream::nProcs());
    label myid = Pstream::myProcNo();
    posLoc[myid] = source().xyz();

    Pstream::allGatherList(posLoc); //make it available for other processors

    vectorField posGlob = ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());
    cellInDomain_.resize(posGlob.size());
    cellInDomain_ = false;
    sourceXyz_.resize(0);
    forAll(posGlob, i)
    {
        point pt = posGlob[i];
        bool indm = false;
        if (target().indomain(pt))
        {
            sourceXyz_.append(pt);
            indm = true;
        }
        cellInDomain_[i] = indm;
     }

    posGlob.resize(0);
    posLoc.clear();

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const DynamicList<scalar>& fld=source().scalarfields_[name];
        label sz=fld.size();

        scalarField  scfldLoc(sz);
        for (label i = 0; i < sz; i++)
        {
            scfldLoc[i] = fld[i];
        }

        List<scalarField> scalar_field(Pstream::nProcs());
        scalar_field[myid] = scfldLoc;
        Pstream::allGatherList(scalar_field);

        scalarField fldGlobal = ListListOps::combine<scalarField>(scalar_field,accessOp<scalarField>());
        label sz1 = fldGlobal.size();
        if (sz1 > cellInDomain_.size())
        {
            sz1 = cellInDomain_.size();
        }


        sourceScalarFields_[name].resize(0);

        for (label i = 0; i < sz1; i++)
        {
            if (cellInDomain_[i])
            {
                sourceScalarFields_[name].append(fldGlobal[i]);
            }
        }

        wait();
        fldGlobal.clear();
        scalar_field.clear();
    }
    //wall distance

    label sz = source().wdists_.size();

    scalarField  scfldLoc(sz);
    for (label i = 0; i < sz; i++)
    {
        scfldLoc[i] = source().wdists_[i]*source().maxWallDist_;//recover absolute value
    }

    List<scalarField> wdists(Pstream::nProcs());
    wdists[myid] = scfldLoc;
    Pstream::allGatherList(wdists);

    scalarField fldGlobal = ListListOps::combine<scalarField>(wdists, accessOp<scalarField>());
    scalar maxwdist = max(fldGlobal);

    fldGlobal /= maxwdist;

    label sz2 = fldGlobal.size();
    if (sz2 > cellInDomain_.size())
    {
        sz2 = cellInDomain_.size();
    }

    sourceWdists_.resize(0);
    for (label i = 0; i < sz2; i++)
    {
        if (cellInDomain_[i])
        {
            sourceWdists_.append(fldGlobal[i]);
        }
    }

    wait();

    fldGlobal.clear();
    wdists.clear();
    //velocity fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const DynamicList<point> &fld=source().vectorfields_[name];
        label sz = fld.size();

        vectorField vectLoc(sz);
        for (label i = 0; i < sz; i++)
        {
            vectLoc[i] = fld[i];
        }

        List<vectorField> vector_field(Pstream::nProcs());
        vector_field[myid] = vectLoc;
        Pstream::allGatherList(vector_field);

        vectorField fldGlobal = ListListOps::combine<vectorField>(vector_field, accessOp<vectorField>());
        label sz1 = fldGlobal.size();
        if (sz1 > cellInDomain_.size())
        {
            sz1 = cellInDomain_.size();
        }

        sourceVectorFields_[name].resize(0);

        for (label i = 0; i < sz1; i++)
        {
            if (cellInDomain_[i])
            {
                sourceVectorFields_[name].append(fldGlobal[i]);
            }
        }
        wait();
        vector_field.clear();
        fldGlobal.clear();
    }
}

void foamMap::generateNodes
(
    const fvMesh *mesh
)
{
    return;
}
