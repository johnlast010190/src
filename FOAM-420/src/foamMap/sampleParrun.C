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
    Foam::samplePar

Group
    grpFoamMap

Description
    Parallel version of foamSample class

SourceFile
    samplePar.C

\*---------------------------------------------------------------------------*/

#include "sampleParrun.H"

using namespace Foam;

samplePar::samplePar()
:
foamSample()
{
    sampleMethod_="nearest";
}

samplePar::samplePar
(
    const Time* runTime,
    const fvMesh *mesh
)
:
foamSample
(
runTime,mesh
)
{
    sampleMethod_="nearest";
}


template<class T>
void Foam::interProcTrans
(
    DynamicList<T> &field,
    const DynamicList<T> &source,
    const word &name
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    //  all send to master
    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        UOPstream toMaster(Pstream::masterNo(), pBufs);
        toMaster << source;
    }

    pBufs.finishedSends();

    if (Pstream::myProcNo() == Pstream::masterNo())
    {
        // Collect my own data
        field.append(source);

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            slave++
        )
        {
            DynamicList<T> posGlob;
            UIPstream fromSlave(slave, pBufs);
            fromSlave >> posGlob;
            field.append(posGlob);
        }
    }

    Info<<"field transferred for :"<<name<<endl;

    return;
}

void samplePar::setOptions
(
    const argList& args
)
{
    foamSample::setOptions(args);

    args.optionReadIfPresent("sampleMethod", sampleMethod_);
}

void samplePar::getDirNames
(
    wordvec &files,
    const word &pattern
)
{
    get_files(files,pattern);
}


void samplePar::getMapTime()
{
    if (mapTimeName_!="latest")
    {
        return;
    }

    word cdir=casedir_+"/processor"+std::to_string(myid_);
    wordvec dirs;

    word pattern=cdir+"/*";

    getDirNames
    (
        dirs,
        pattern
    );

    scalar latestTime=0;

    for (unsigned i=0;i<dirs.size();i++)
    {
        if (!is_dir(dirs[i].c_str()))
        {
            continue;
        }
        wordvec vtmp;
        stringtok(vtmp,dirs[i],"/");
        word wlast=vtmp[vtmp.size()-1];
        label itime=atoi(wlast.c_str());
        latestTime=max(latestTime,itime);
    }

    mapTimeName_=name(latestTime);
    return;
}

void samplePar::getBodySurface()
{
    if (!Pstream::parRun())
    {
        return;
    }
    getbodyPatchNames();
    List<vectorField> posLoc(Pstream::nProcs());
    label myid=Pstream::myProcNo();
    DynamicList<point> xyz;

    forAll(mesh->boundaryMesh(),r)
    {
        word type=mesh->boundaryMesh()[r].type();
        if (type == "empty" || type == "symmetryPlane")
        {
            continue;
        }
        word name=mesh->boundaryMesh()[r].name();
        if (!inlist(name,bodyPatchNames_))
        {
            continue;
        }
        forAll(mesh->boundaryMesh()[r],j)
        {
            const point &pt=mesh->Cf().boundaryField()[r][j];
            xyz.append(pt);
        }
    }

    posLoc[myid] =xyz;

    Pstream::allGatherList(posLoc);

    bodySurfCenter_=ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());

}

void samplePar::combineFieldsBK()
{
    if (!Pstream::parRun())
    {
        return;
    }
    List<vectorField> posLoc(Pstream::nProcs());
    label myid=Pstream::myProcNo();
    DynamicList<point> xyz;
    xyz.setSize(mesh->C().size());
    forAll(mesh->C(),i)
    {
        xyz[i]=mesh->C()[i];
    }
    posLoc[myid] =xyz;

    Pstream::allGatherList(posLoc); //make it available for other processors
    if (master())
    {
        sourceXyz_=ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());
    }

    wait();
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();
        const volScalarField &sfield=source().getScalarField(name);
        label sz=sfield.primitiveField().size();
        scalarField  scfldLoc(sz);
        for (label i=0;i<sz;i++)
        {
            scfldLoc[i]=sfield[i];
        }

        List<scalarField> scalar_field(Pstream::nProcs());
        scalar_field[myid]=scfldLoc;
        Pstream::allGatherList(scalar_field);
        scalarField fldGlobal=ListListOps::combine<scalarField>(scalar_field,accessOp<scalarField>());

        sourceScalarFields_[name].resize(0);
        if (master())
        {
            sourceScalarFields_[name]=fldGlobal;
        }
        wait();
        fldGlobal.resize(0);
    }

    //velocity fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const volVectorField &vfield=source().getVectorField(name);
        label sz=vfield.primitiveField().size();

        vectorField  vectLoc(sz);
        for (label i=0;i<sz;i++)
        {
            vectLoc[i]=vfield[i];
        }
        List<vectorField> vector_field(Pstream::nProcs());
        vector_field[myid]=vectLoc;
        Pstream::allGatherList(vector_field);
        vectorField vfldGlobal=ListListOps::combine<vectorField>(vector_field,accessOp<vectorField>());
        sourceVectorFields_[name].resize(0);
        if (master())
        {
            sourceVectorFields_[name]=vfldGlobal;
        }

        wait();
        vfldGlobal.resize(0);
    }
    return;
}


void samplePar::combineFields()
{
    if (!Pstream::parRun())
    {
        return;
    }

    // cell centre coordinates

    interProcTrans
    (
        sourceXyz_,
        source().xyz(),
        "cellVerts"
    );

    //scalar fields
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();
        const DynamicList<scalar> &scfld=source().scalarfields_[name];
        interProcTrans
        (
            sourceScalarFields_[name],
            scfld,
            name
        );
    }

      //vector fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();
        const DynamicList<point> &vecfld=source().vectorfields_[name];
        interProcTrans
        (
            sourceVectorFields_[name],
            vecfld,
            name
        );
    }

    return;
}

void samplePar::constructKnn()
{
    if (!master())
    {
        return;
    }

    label ncells=sourceXyz_.size();
    label cnt=0;
    farray2d trainData;
    for (label i=0;i<ncells;i++)
    {
        const point &pt=sourceXyz_[i];
        std::vector<scalar> pts(3);

        pts[0]=pt.x();
        pts[1]=pt.y();
        pts[2]=pt.z();
        trainData.push_back(pts);
        cellMap_[cnt]=i;
        cnt++;
    }

    Info<<"build knn..."<<trainData.size()<<endl;
    internalMap_.build(trainData);

    trainData.resize(0);
    for (label i=0;i<bodySurfCenter_.size();i++)
    {
        const point &pt=bodySurfCenter_[i];
        std::vector<scalar> pts(3);

        pts[0]=pt.x();
        pts[1]=pt.y();
        pts[2]=pt.z();
        trainData.push_back(pts);
    }

    Info<<"build object knn..."<<trainData.size()<<endl;
    bodyBoundaryTree_.build(trainData);
    return;
}

void samplePar::searchNodeField
(
    gridField &grid,
    const word &funct
)
{
    if (!master())
    {
        return;
    }

    grid.setFieldSize();//set placeholders

    label cnt=0;
    label cnt1=0;
    label cnt2=0;


    for (label i=0;i<grid.nx_+1;++i)
    {
        for (label j=0;j<grid.ny_+1;++j)
        {
            for (label k=0;k<grid.nzpt();++k)
            {
                const point &pt=grid.gridNodes_(i,j,k);
                std::vector<scalar> qvect(3);
                qvect[0]=pt.x();
                qvect[1]=pt.y();
                qvect[2]=pt.z();
                std::vector<label> nb;
                std::vector<scalar> dist;

                internalMap_.search
                (
                    qvect,
                    nb,
                    dist
                );

                //search for the nearest body surface face center
                std::vector<label> nb1;
                std::vector<scalar> dist1;
                if (bodyBoundaryTree_.active)
                {
                    bodyBoundaryTree_.search
                    (
                        qvect,
                        nb1,
                        dist1
                    );
                }

                // check if a grid node is in solid
                label insolid=1;

                if (nb1.size()==0)
                {
                    insolid=0;//cannot find boundary boundary match
                    Info<<"cannot find boundary boundary match..."<<endl;

                }
                else
                {
                    //Info<<"dist comp..."<<dist1[0]<<" "<<dist[0]<<endl;
                    if (dist1[0]>dist[0])
                    {
                        insolid=0;//nearer to fluid cell
                        cnt1++;
                    }
                    else
                    {
                        // dist1[1]<=dist[0] can also be a fluid cell if cell_alpha>1.0e-3
                        label bid=nb1[0];
                        const point &pf=bodySurfCenter_[bid];

                        label ic=nb[0];
                        point pc=sourceXyz_[ic];
                        scalar delta=(pc-pt)&(pf-pt);
                        if (delta>0)
                        {
                            insolid=1;
                        }
                        else
                        {
                            insolid=0;//delta <0, in fluid
                            cnt2++;
                        }

                    }
                }

                grid.insolids_(i,j,k)=insolid;
                cnt+=insolid;

                if (funct=="field")
                {
                    if (grid.scalarSamples_.size()>=1)
                    {
                        getGridField(grid,i,j,k,nb,dist,insolid); //get the scalar field at a nodal point
                    }
                    if (grid.vectorSamples_.size()>=1)
                    {
                        getGridVectorField(grid,i,j,k,nb,dist,insolid);
                    }
                }


            }
        }
    }

    Info<<"number of node in solids:"<<cnt<<" "<<cnt1<<" "<<cnt2<<endl;


    return;
}
