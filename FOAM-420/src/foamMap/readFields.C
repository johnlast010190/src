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
    Foam::readFields

Group
    grpfoamMap

Description
    Base class for mapping two unstructured flow fields or sampling data
    from unstructured fields into uniform grids

SourceFile
    readFields.C

\*---------------------------------------------------------------------------*/

#include "readFields.H"
#include "foamMap.H"
#include "db/IOobjectList/IOobjectList.H"

using namespace Foam;

readFields::readFields()
{
    type_ = "source";
    norm_ = false;
}

readFields::readFields
(
    const fvMesh* _mesh,
    const Time* runTime,
    foamMap* map,
    const word& type
)
{
    norm_       = false;
    mesh        = _mesh;
    runTime_    = runTime;
    type_       = type;
    fieldmap_   = map;
    scfields_   = fieldmap_->mapScalarFields_;
    vecfields_  = fieldmap_->mapVectorFields_;
    UrotDegreeFromSource_ = 0;
    rhoRef_     = 1;
    alphaMax_   = fieldmap_->alphaMax_;
    interp_     = fieldmap_->interp_;
    Uref_       = point(1,0,0);
    ncells_     = mesh->C().size();
    mapTime_    = runTime_->timeName();

    if (fieldmap_->mapTimeName_.length() >= 1)
    {
        mapTime_ = fieldmap_->mapTimeName_;
    }

    boundBox meshBb(mesh->points(), true);
    getBoundBox(meshBb);
}

void readFields::createFields(const word& tname)
{
    word timeName(tname);
    if (timeName.length() < 1)
    {
        timeName = runTime_->timeName();
    }

    IOobjectList objects(*mesh, timeName);

    HashSet<word> removeScalars;
    HashSet<word> removeVectors;

    forAllConstIter(HashSet<word>, fieldmap_->mapScalarFields_, iter)
    {
        const word& scname = iter();

        Info<< "Reading scalar field: "<<scname;

        if (objects.lookup(scname) == nullptr)
        {
            Info<< " ...not found" << endl;

            removeScalars.insert(scname);

            continue;
        }

        Info<< endl;

        autoPtr<volScalarField> scfieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    scname,
                    timeName,
                    *mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh
            )
        );

        scfieldPtr.ptr()->store();
    }

    forAllConstIter(HashSet<word>, fieldmap_->mapVectorFields_, iter)
    {
        const word& vecname = iter();

        Info<< "Reading vector field: " << vecname;

        if (objects.lookup(vecname) == nullptr)
        {
            Info<< " ...not found" << endl;

            removeVectors.insert(vecname);

            continue;
        }

        Info<< endl;

        autoPtr<volVectorField> vecfieldPtr
        (
            new volVectorField
            (
                IOobject
                (
                    vecname,
                    timeName,
                    *mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh
            )
        );

        vecfieldPtr.ptr()->store();
    }


    fieldmap_->mapScalarFields_ -= removeScalars;
    fieldmap_->mapVectorFields_ -= removeVectors;

    //flux map
    forAllConstIter(HashSet<word>, fieldmap_->mapSurfaceScalarFields_, iter)
    {
        const word& scname = iter();

        autoPtr<surfaceScalarField> scfieldPtr
        (
            new surfaceScalarField
            (
                IOobject
                (
                    scname,
                    timeName,
                    *mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh
            )
        );

        scfieldPtr.ptr()->store();
    }

    return;
}

void readFields::storeFields()
{
    if (type_ != "source") return;
    label nbface = numBoundFaces();

    label ncell = mesh->C().size();
    label sz = ncell;
    if (fieldmap_->mapBoundary_)
    {
        sz += nbface;
    }
    xyz_.setSize(sz);
    forAll(mesh->C(), i)
    {
        xyz_[i] = mesh->C()[i];
    }

    if (fieldmap_->mapBoundary_)
    {
        label cnt = ncell;
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if (type == "empty" || type == "symmetryPlane")
            {
                continue;
            }
            const vectorField& bndPoints = mesh->Cf().boundaryField()[r];
            forAll(bndPoints, j)
            {
                const point& pt = bndPoints[j];
                xyz_[cnt] = pt;
                cnt++;
            }
        }
    }

    forAllConstIter(HashSet<word>, scfields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);

        const volScalarField& sfield = getScalarField(name);
        label sz = sfield.primitiveField().size();
        if (fieldmap_->mapBoundary_)
        {
            sz += nbface;
        }

        scalarfields_[name].resize(sz);

        for (label i = 0; i < ncell; i++)
        {
            scalarfields_[name][i] = sfield.primitiveField()[i]/coef;
        }

        if (fieldmap_->mapBoundary_)
        {
            //boundaries
            label cnt = ncell;
            forAll(mesh->boundaryMesh(), r)
            {
                word type = mesh->boundaryMesh()[r].type();
                if (type == "empty" || type == "symmetryPlane")
                {
                    continue;
                }

                forAll(mesh->boundaryMesh()[r], j)
                {
                    scalar sc = sfield.boundaryField()[r][j];
                    scalarfields_[name][cnt] = sc/coef;
                    cnt++;
                }
            }
        }
    }

    //vector fields
    forAllConstIter(HashSet<word>, vecfields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);
        const volVectorField& vsfield = getVectorField(name);
        label sz = vsfield.primitiveField().size();
        if (fieldmap_->mapBoundary_)
        {
            sz += nbface;
        }

        vectorfields_[name].resize(sz);
        for (label i = 0; i < ncell; i++)
        {
            vectorfields_[name][i] = vsfield.primitiveField()[i]/coef;
        }

        if (fieldmap_->mapBoundary_)
        {
            //boundaries
            label cnt = ncell;
            forAll(mesh->boundaryMesh(), r)
            {
                word type = mesh->boundaryMesh()[r].type();
                if (type == "empty" || type == "symmetryPlane")
                {
                    continue;
                }
                forAll(mesh->boundaryMesh()[r], j)
                {
                    point vsc = vsfield.boundaryField()[r][j];
                    vectorfields_[name][cnt] = vsc/coef;
                    cnt++;
                }
            }
        }
    }
    return;
}

void readFields::storeFields
(
    const boundBox& sampleBox
)
{
    if (type_ != "source") return;
    xyz_.resize(0);
    std::vector<label> incells;
    typedef Tuple2<label, label> regface;
    std::vector<regface> bfaces;

    forAll(mesh->C(), i)
    {
        const point& pt = mesh->C()[i];
        if (sampleBox.contains(pt))
        {
            xyz_.append(pt);
            incells.push_back(i);
        }
    }

    if (fieldmap_->mapBoundary_)
    {
        forAll(mesh->boundaryMesh(), r)
        {
            word type = mesh->boundaryMesh()[r].type();
            if (type == "empty" || type == "symmetryPlane")
            {
                continue;
            }
            const vectorField& bndPoints = mesh->Cf().boundaryField()[r];
            forAll(bndPoints, j)
            {
                const point& pt = bndPoints[j];
                if (sampleBox.contains(pt))
                {
                    xyz_.append(pt);
                    regface face(r, j);
                    bfaces.push_back(face);
                }
            }
        }
    }

    //scalar fields
    forAllConstIter(HashSet<word>, scfields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);
        //internal field
        const volScalarField& sfield = getScalarField(name);
        label ncell = incells.size();
        label nbface = bfaces.size();
        scalarfields_[name].resize(ncell + nbface);
        for (label i = 0; i < ncell; i++)
        {
            label ic = incells[i];
            scalarfields_[name][i] = sfield.primitiveField()[ic]/coef;
        }

        //boundaries
        if (fieldmap_->mapBoundary_)
        {
            for (label jf = 0; jf < nbface; jf++)
            {
                const regface& face = bfaces[jf];
                label r = face.first();
                label j = face.second();
                scalar sc = sfield.boundaryField()[r][j];
                scalarfields_[name][ncell + jf] = sc/coef;
            }
        }
    }//isc

    //vector fields
    forAllConstIter(HashSet<word>, vecfields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);
        //internal field
        const volVectorField& vsfield = getVectorField(name);
        label ncell = incells.size();
        label nbface = bfaces.size();
        vectorfields_[name].resize(ncell + nbface);

        for (label i = 0; i < ncell; i++)
        {
            label ic = incells[i];
            vectorfields_[name][i] = vsfield.primitiveField()[ic]/coef;
        }

        //boundaries
        if (fieldmap_->mapBoundary_)
        {
            for (label jf = 0; jf < nbface; jf++)
            {
                const regface& face = bfaces[jf];
                label r = face.first();
                label j = face.second();
                point vsc = vsfield.boundaryField()[r][j];
                vectorfields_[name][ncell+jf] = vsc/coef;
            }
        }
    }

    return;
}

void readFields::setWallDist
(
    const volScalarField& y,
    scalar maxwdist
)
{
    maxWallDist_ = maxwdist;
    label ncells = y.primitiveField().size();
    label nbface = numBoundFaces();
    label sz = ncells;
    if (fieldmap_->mapBoundary_)
    {
        sz += nbface;
    }

    wdists_.setSize(sz);
    wdists_ = 1;
    for (label i = 0; i < ncells; i++)
    {
        wdists_[i] = y.primitiveField()[i]/maxWallDist_;
    }
    deltaY_ = 1.0/float(fieldmap_->nwdist_);
}

void readFields::buildKdTrees
(
    const volScalarField& y,
    scalar cutoff, /*= 0.8*/
    const word& excludeBndType /*= none*/
)
{
    maxWallDist_ = gMax(y.primitiveField());

    label ncells = y.primitiveField().size();
    wdists_.setSize(ncells);
    for (label i = 0; i < ncells; i++)
    {
        wdists_[i] = y.primitiveField()[i]/maxWallDist_;
    }
    deltaY_ = 1.0/float(fieldmap_->nwdist_);

    Array2d<point> xyz;
    xyz.resize(fieldmap_->nwdist_ + 1);

    alphaMap_.resize(fieldmap_->nwdist_ + 1);
    forAll(mesh->C(), i)
    {
        point pt = mesh->C()[i];
        if (wdists_[i] > fieldmap_->alphaMax_)
        {
            continue;
        }

        label ic = inGroup(i);
        if (ic < 0) ic = 0;
        if (ic > fieldmap_->nwdist_) ic = fieldmap_->nwdist_;
        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());
        point pt1(xb, yb, zb);
        xyz[ic].push_back(pt1);
        alphaMap_[ic].push_back(i);
    }

    if (type_ == "source")
    {
        kdTrees_.resize(fieldmap_->nwdist_ + 1);

        for (unsigned ic = 0; ic < xyz.size(); ic++)
        {
            label npt = xyz[ic].size();
            if (npt < 1) continue;
            farray2d trainData(npt);
            for (label i = 0; i < npt; i++)
            {
                trainData[i].resize(3);
                point pt = xyz[ic][i];
                trainData[i][0] = pt.x();
                trainData[i][1] = pt.y();
                trainData[i][2] = pt.z();
            }
            kdTrees_[ic].build(trainData);
        }
        //boundary search tree
        construct_internal_knn(cutoff);

        if (fieldmap_->mapBoundary_)
        {
            construct_boundary_knn();
        }
    }

    return;
}


void readFields::getBoundBox
(
    const boundBox& box
)
{
    const point& pmin = box.min();
    const point& pmax = box.max();
    xmin_ = pmin.x();
    ymin_ = pmin.y();
    zmin_ = pmin.z();
    xmax_ = pmax.x();
    ymax_ = pmax.y();
    zmax_ = pmax.z();
}

readFields::~readFields(){;}

void readFields::setInputs()
{
    if (type_ == "source")
    {
        Uref_ = fieldmap_->UrefSource_;
        rhoRef_ = fieldmap_->rhoRefSource_;
    }
    else
    {
        Uref_ = fieldmap_->UrefTarget_;
        rhoRef_ = fieldmap_->rhoRefTarget_;
    }
    interp_ = fieldmap_->interp_;
    return;
}

scalar readFields::scale(const word& fldname)
{
    const word unknown("unknown");
    const word& ftype = fieldmap_->fieldTypes_.lookup(fldname, unknown);

    if (ftype == "velocity")
    {
        return U0();
    }
    else if (ftype == "pressure")
    {
        return pref();
    }
    else if (ftype == "turbEnergy")
    {
        return U0()*U0();
    }

    return 1;
}


label readFields::numBoundFaces()
{
    label nbface = 0;
    forAll(mesh->boundaryMesh(), r)
    {
        label sz = mesh->boundaryMesh()[r].size();
        word type = mesh->boundaryMesh()[r].type();
        if (type == "empty" || type == "symmetryPlane")
        {
            continue;
        }
        nbface += sz;
    }
    return nbface;
}

void readFields::construct_internal_knn(scalar cutoff)
{
    Info<< "constructing internal KNN map..." << endl;
    label cnt = 0;
    cutoff = -0.1;

    farray2d trainData;
    //for gridMap (cutoff<0), no need to do alpha-mapping, internal map contains all cells
    forAll(mesh->C(), i)
    {
        if (cutoff >= 0)
        {
            scalar alphai = wdists_[i];
            if (alphai > cutoff)
            {
                continue;
            }
        }
        point pt = mesh->C()[i];
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
    Info<< "traindata size:" << trainData.size() << " " << cnt << endl;
    internalMap_.build(trainData);
    Info<< "internal map build." << endl;
    return;
}


void readFields::construct_boundary_knn()
{
    farray2d trainData;

    forAll(mesh->boundaryMesh(), r)
    {
        word type = mesh->boundaryMesh()[r].type();
        if (type == "empty" || type == "symmetryPlane")
        {
            continue;
        }

        word name = mesh->boundaryMesh()[r].name();

        forAll(mesh->boundaryMesh()[r], j)
        {
            point pt = mesh->Cf().boundaryField()[r][j];
            scalar xb = xbar(pt.x());
            scalar yb = ybar(pt.y());
            scalar zb = zbar(pt.z());
            std::vector<scalar> pts(3);
            pts[0] = xb;
            pts[1] = yb;
            pts[2] = zb;
            trainData.push_back(pts);
        }
    }

    if (trainData.size() >= 1)
    {
        boundaryMap_.build(trainData);
    }
    return;
}



void readFields::transformPoints
(
    const boundBox& tbox,
    DynamicList<point>& xyz
)
{
    //called after construct boundBox, note that tbox must be the original (not extended)
    if (myType() != "source")
    {
        return;
    }

    scalar xmint = tbox.min().x();
    scalar xmaxt = tbox.max().x();
    scalar xmins = gbox().min().x();
    scalar xmaxs = gbox().max().x();
    scalar xcoef = (xmaxt - xmint)/(xmaxs - xmins);

    scalar ymint = tbox.min().y();
    scalar ymaxt = tbox.max().y();
    scalar ymaxs = gbox().max().y();
    scalar ymins = gbox().min().y();
    scalar ycoef = (ymaxt - ymint)/(ymaxs - ymins);

    scalar zmint = tbox.min().z();
    scalar zmaxt = tbox.max().z();
    scalar zmaxs = gbox().max().z();
    scalar zmins = gbox().min().z();
    scalar zcoef = (zmaxt - zmint)/(zmaxs - zmins);

    forAll(xyz, i)
    {
        scalar x = xyz[i].x();
        x = xmint + (x - xmins)*xcoef;
        xyz[i].x() = x;

        scalar y = xyz[i].y();
        y = ymint + (y - ymins)*ycoef;
        xyz[i].y() = y;

        scalar z=xyz[i].z();
        z = zmint + (z - zmins)*zcoef;
        xyz[i].z() = z;
    }
}

void readFields::constructBoundBox
(
    scalar extcoef,
    const word type,
    const label numCells
)
{
    //global box
    gbox_ = new boundBox(mesh->C(), true);

    // local box
    boundBox lbox(mesh->C(), false);

    if (extcoef > SMALL)
    {
        lbox.inflate(extcoef);
    }

    bbox_ = new smBoundBox(lbox.min(), lbox.max(), type);

    if (type == "smBoundBox")
    {
        bbox_().createBox(mesh, numCells);
    }
}

