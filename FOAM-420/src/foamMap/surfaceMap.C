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
    Foam::surfaceMap

Group
    grpfoamMap

Description
    A class designed to create a symmetric scalar field from an unsymmetric field

SourceFile
    surfaceMap.C

\*---------------------------------------------------------------------------*/


#include "surfaceMap.H"
using namespace Foam;

void surfaceMap::defaults()
{
    imap_=true;
    bmap_=true;
    byType_=false;
    byNameList_=true;
    byPhysicalType_=false;
    bodySurfaceType_="wall";
}

surfaceMap::surfaceMap()
:
    foamMap()
{
    defaults();
}

surfaceMap::surfaceMap
(
    const fvMesh* mesh,
    const Time* runTime
)
:
    foamMap(),
    mesh_(mesh),
    runTime_(runTime)
{
    defaults();
}

surfaceMap::surfaceMap
(
    const fvMesh& mesh
)
:
    foamMap(),
    mesh_(&mesh),
    runTime_(&mesh.time())
{
    defaults();
}
void surfaceMap::setInput
(
    const dictionary& dict
)
{
    foamMap::setInput(dict);

    word identifyObject=dict.lookupOrDefault<word> ("identifyObject","byNameList");

    //Info<<"idenfyObject:"<<identifyObject<<endl;
    if (identifyObject=="byType")
    {
        byType_=true;
    }

    else if (identifyObject=="byPhysicalType")
    {
        byPhysicalType_=true;
    }
    else
    {
        byNameList_=true;
    }

    //Info<<"byType:"<<byType_<<endl;

    bodySurfaceType_=dict.lookupOrDefault<word> ("bodySurfaceType","wall");
    //Info<<"bodySurfaceType:"<<bodySurfaceType_<<endl;

    imap_=dict.lookupOrDefault<bool> ("internalSym",true);

    bmap_=dict.lookupOrDefault<bool> ("surfaceSym",false);

    point p0(0,0,0);
    vector norm(0,1,0);
    plane plane0(p0,norm);

    plane newplane=dict.lookupOrDefault<plane> ("mirror",plane0);

    mirror_=new plane(newplane);

    if (byType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(),r)
        {
            word type=mesh_->boundaryMesh()[r].type();
            word name=mesh_->boundaryMesh()[r].name();
            if (type==bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }

    }

    else if (byPhysicalType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(),r)
        {
            word type=mesh_->boundaryMesh()[r].physicalType();
            word name=mesh_->boundaryMesh()[r].name();
            if (type==bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }

    }
    else
    {
        const polyBoundaryMesh& bMesh = mesh_->boundaryMesh();

        labelHashSet includePatches;

        if (dict.found("patches"))
        {
            includePatches = bMesh.patchSet
                (
                    wordReList(dict.lookup("patches"))
                );
        }

        const wordList allPatchNames(bMesh.names());

        forAllConstIter(labelHashSet, includePatches, it)
        {
            sourcePatchNames_.push_back
                (
                    allPatchNames[*it]
                );
        }
    }

    Info<<"number of selected patch regions:"<<sourcePatchNames_.size()<<endl;
}
void surfaceMap::setInput
(
    const dictionary& dict,
    const argList& args
)
{
    foamMap::setInput(dict, args);

    word identifyObject=dict.lookupOrDefault<word> ("identifyObject","byNameList");

    //Info<<"idenfyObject:"<<identifyObject<<endl;
    if (identifyObject=="byType")
    {
        byType_=true;
    }

    else if (identifyObject=="byPhysicalType")
    {
        byPhysicalType_=true;
    }
    else
    {
        byNameList_=true;
    }

    //Info<<"byType:"<<byType_<<endl;

    bodySurfaceType_=dict.lookupOrDefault<word> ("bodySurfaceType","wall");
    //Info<<"bodySurfaceType:"<<bodySurfaceType_<<endl;

    imap_=dict.lookupOrDefault<bool> ("internalSym",true);

    bmap_=dict.lookupOrDefault<bool> ("surfaceSym",false);

    point p0(0,0,0);
    vector norm(0,1,0);
    plane plane0(p0,norm);

    plane newplane=dict.lookupOrDefault<plane> ("mirror",plane0);

    mirror_=new plane(newplane);

    if (byType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(),r)
        {
            word type=mesh_->boundaryMesh()[r].type();
            word name=mesh_->boundaryMesh()[r].name();
            if (type==bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }

    }

    else if (byPhysicalType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(),r)
        {
            word type=mesh_->boundaryMesh()[r].physicalType();
            word name=mesh_->boundaryMesh()[r].name();
            if (type==bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }

    }
    else
    {
        const polyBoundaryMesh& bMesh = mesh_->boundaryMesh();

        labelHashSet includePatches;

        if (dict.found("patches"))
        {
            includePatches = bMesh.patchSet
            (
                wordReList(dict.lookup("patches"))
            );
        }

        const wordList allPatchNames(bMesh.names());

        forAllConstIter(labelHashSet, includePatches, it)
        {
            sourcePatchNames_.push_back
            (
                allPatchNames[*it]
            );
        }
    }

    Info<<"number of selected patch regions:"<<sourcePatchNames_.size()<<endl;

}

void surfaceMap::getParallelFields()
{
    if (!Pstream::parRun())
    {
        return;
    }

    List<vectorField> posLoc(Pstream::nProcs());
    label myid=Pstream::myProcNo();

    if (internalSym())
    {
        DynamicList<point> xyz;
        xyz.resize(mesh_->C().size());
        forAll(mesh_->C(),i)
        {
            xyz[i]=mesh_->C()[i];
        }

        posLoc[myid] =xyz;

        Pstream::allGatherList(posLoc); //make it available for other processors
        vectorField posGlob=ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());

        sourceXyz_=posGlob;
        posGlob.resize(0);
    }
    Info<<"source_xyz:"<<sourceXyz_.size()<<endl;


    if (bndSym())
    {
        //boundary
        List<vectorField> posLoc(Pstream::nProcs());

        for (unsigned ir=0; ir<sourcePatchNames_.size();ir++)
        {
            const word &bname=sourcePatchNames_[ir];
            label r=patchId(bname);
            if (r<0)
            {
                 WarningInFunction<<"Warning: patch id "<<bname<<" not found.\n";
                 continue;
            }

            DynamicList<point> bxyz;
            forAll(mesh_->Cf().boundaryField()[r],j)
            {
                const point &pt=mesh_->Cf().boundaryField()[r][j];
                bxyz.append(pt);
            }

            posLoc[myid]=bxyz;
            Pstream::allGatherList(posLoc); //make it available for other processors
            vectorField posGlob=ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());
            bsource_xyz_[bname]=posGlob;

        }
    }

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();

        const volScalarField &fld=source().getScalarField(scname);
        label sz=fld.primitiveField().size();
        if (internalSym())
        {
            scalarField  scfldLoc(sz);
            forAll(fld.primitiveField(),i)
            {
                scfldLoc[i]=fld[i];
            }
            List<scalarField> scalar_field(Pstream::nProcs());
            scalar_field[myid]=scfldLoc;
            Pstream::allGatherList(scalar_field);
            scalarField fldGlobal=ListListOps::combine<scalarField>(scalar_field,accessOp<scalarField>());
            sourceScalarFields_[scname]=fldGlobal;

            wait();
        }


        if (bndSym())
        {
            for (unsigned ir=0; ir<sourcePatchNames_.size();ir++)
            {
                word pname=sourcePatchNames_[ir];
                word key=pname+','+scname;
                List<scalarField> bscalar_field(Pstream::nProcs());//local field on boundary
                label pid=patchId(pname);
                DynamicList<scalar> sclist;
                forAll(fld.boundaryField()[pid],j)
                {
                    scalar scb=fld.boundaryField()[pid][j];
                    sclist.append(scb);
                }
                bscalar_field[myid]=sclist;
                Pstream::allGatherList(bscalar_field);
                scalarField bfldGlobal=ListListOps::combine<scalarField>(bscalar_field,accessOp<scalarField>());
                bsourceScalarFields_[key]=bfldGlobal;
                bfldGlobal.resize(0);
            }
        }

    }//for isc
}

void surfaceMap::createMirrorFields
(
    const word &timeName
)
{
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();
        word mScalar=mirrorName(scname);
        volScalarField &fld=source().getScalarField(scname);

        autoPtr<volScalarField> scfieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    mScalar,
                    timeName,
                    *mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                //fld
                *mesh_,
                fld.dimensions()
            )
        );

        scfieldPtr.ptr()->store();
        volScalarField &mfld=source().getScalarField(mScalar);
        mfld=fld;
    }
}

void surfaceMap::buildParSearchTree()
{
    if (internalSym())
    {
        farray2d trainData;
        forAll(sourceXyz_,i)
        {
            const point &pt=sourceXyz_[i];
            scalar xb=source().xbar(pt.x());
            scalar yb=source().ybar(pt.y());
            scalar zb=source().zbar(pt.z());
            std::vector<scalar> pts(3);
            pts[0]=xb;
            pts[1]=yb;
            pts[2]=zb;
            trainData.push_back(pts);
        }

        internalMap_.build(trainData);
        wait();
    }

    if (bndSym())
    {
        std::map<word,DynamicList<point>>::iterator it;

        for (it=bsource_xyz_.begin();it!=bsource_xyz_.end();it++)
        {
            word bname=it->first;
            farray2d tData;
            const DynamicList<point> &bflist=it->second;
            for (label j=0;j<bflist.size();j++)
            {
                const point &pt=bflist[j];
                std::vector<scalar> pts(3);
                pts[0]=pt.x();
                pts[1]=pt.y();
                pts[2]=pt.z();
                tData.push_back(pts);
                faceIdMaps_[bname].push_back(j);
            }

            bndMaps_[bname].build(tData);
        }
        wait();
    }

}

void surfaceMap::buildSearchTree()
{
    if (Pstream::parRun())
    {
        buildParSearchTree();
        return;
    }

    if (internalSym())
    {
        farray2d trainData;
        forAll(mesh_->C(),i)
        {
            const point &pt=mesh_->C()[i];
            scalar xb=source().xbar(pt.x());
            scalar yb=source().ybar(pt.y());
            scalar zb=source().zbar(pt.z());
            std::vector<scalar> pts(3);
            pts[0]=xb;
            pts[1]=yb;
            pts[2]=zb;
            trainData.push_back(pts);
        }

        Info<<"Number of data samples:"<<trainData.size()<<endl;
        internalMap_.build(trainData);
    }

    if (bndSym())
    {
        for (unsigned ir=0; ir<sourcePatchNames_.size();ir++)
        {
            const word &bname=sourcePatchNames_[ir];
            label r=patchId(bname);
            if (r<0)
            {
                 WarningInFunction<<"Warning, patch id: "<<bname<<" not found.\n";
                 continue;
            }

            farray2d tData;
            std::vector<label> faceIds;

            forAll(mesh_->Cf().boundaryField()[r], j)
            {
                const point &pt=mesh_->Cf().boundaryField()[r][j];
                std::vector<scalar> pts(3);
                pts[0]=pt.x();
                pts[1]=pt.y();
                pts[2]=pt.z();
                tData.push_back(pts);
                faceIdMaps_[bname].push_back(j);
            }

            bndMaps_[bname].build(tData);
        }

        Info<<"Number of boundary maps:"<<bndMaps_.size()<<endl;
    }

    return;
}

label surfaceMap::patchId(const word &name)
{
    forAll(mesh_->boundaryMesh(),r)
    {
        if (mesh_->boundaryMesh()[r].name()==name)
        {
            return r;
        }
    }
    return -1;
}

void surfaceMap::getSymFields()
{
    Info<<"Obtain the symmetric fields..."<<endl;
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();

        const volScalarField &scfield=source().getScalarField(scname);
        word mscname=mirrorName(scname);
        volScalarField &mscfield=source().getScalarField(mscname);
        mscfield=0.5*(mscfield+scfield);
    }

}

void surfaceMap::getMirrorFields()
{
    std::vector<label> ulist;
    std::map<word,std::vector<label>> ublist;

    if (internalSym())
    {
        ulist.resize(mesh_->C().size());
        forAll(mesh_->C(),i)
        {
            const point &pt=mesh_->C()[i];
            point mpt=mirror().mirror(pt);
            scalar xc=source().xbar(mpt.x());
            scalar yc=source().ybar(mpt.y());
            scalar zc=source().zbar(mpt.z());
            std::vector<scalar> qvect(3);
            qvect[0]=xc;
            qvect[1]=yc;
            qvect[2]=zc;
            std::vector<label> nb;
            std::vector<scalar> dist;
            internalMap_.search
            (
                qvect,
                nb,
                dist
            );
            label ic=nb[0];
            ulist[i]=ic;
        }
    }

    if (bndSym())
    {
        for (unsigned r=0;r<sourcePatchNames_.size();r++)
        {
            const word &pname=sourcePatchNames_[r];
            label ir=patchId(pname);
            if (ir<0)
            {
                continue;
            }

            forAll(mesh_->boundaryMesh()[ir],j)
            {
                const point &pt=mesh_->Cf().boundaryField()[ir][j];
                point mpt=mirror().mirror(pt);
                std::vector<scalar> pts(3);
                pts[0]=mpt.x();
                pts[1]=mpt.y();
                pts[2]=mpt.z();
                std::vector<label> nb;
                std::vector<scalar> dist;
                bndMaps_[pname].search
                (
                    pts,
                    nb,
                    dist
                );
                label ib=nb[0];
                label jc=faceIdMaps_[pname][ib];

                ublist[pname].push_back(jc);
            }
        }
    }//if bnd
    wait();

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();
        const volScalarField &scfield=source().getScalarField(scname);
        word mscname=mirrorName(scname);
        volScalarField &mscfield=source().getScalarField(mscname);
        volScalarField dsc( scfield-mscfield );

        if (ulist.size()>0)
        {
            for (unsigned i=0;i<ulist.size();++i)
            {
                label ic=ulist[i];
                scalar msc;
                if (!Pstream::parRun())
                {
                    msc=scfield.primitiveField()[ic];
                }
                else
                {
                    msc=sourceScalarFields_[scname][ic];
                }
                mscfield.primitiveFieldRef()[i]=msc;
            }
        }

        if (ublist.size()>0)
        {
            Info<<"search boundary field..."<<endl;
            std::map<word,std::vector<label>>::iterator it;
            for (it=ublist.begin();it !=ublist.end();it++)
            {
                word pname=it->first;
                word pkey=pname+','+scname; //parallel
                label ir=patchId(pname);
                if (ir<0)
                {
                    Info<<"error occurred in boundary mapping."<<endl;
                    continue;
                }
                const std::vector<label> &blist=it->second;
                for (unsigned j=0;j<blist.size();j++)
                {
                    label jc=blist[j];

                    scalar sc;
                    if (!Pstream::parRun())
                    {
                        sc=scfield.boundaryField()[ir][jc];
                    }
                    else
                    {
                        sc=bsourceScalarFields_[pkey][jc];
                    }
                    mscfield.boundaryFieldRef()[ir][j]=sc;
                }
            }
        }

    }//for isc
}
