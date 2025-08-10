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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::foamSample

Group
    grpfoamMap

Description
    A set of classes designed to sample an arbitrary number of flow fields from an unstructured
    mesh region into a uniform grid, or another unstructured mesh region.

SourceFile
    blockSample.C

\*---------------------------------------------------------------------------*/

#include "blockSample.H"
#include <cstdlib>
#include <iomanip>

using namespace Foam;
// get files or sub-directories

void Foam::get_files
(
    wordvec& v,
    const word& root,
    const word &filter
)
{
    DIR* dirp = opendir(root.c_str());
    struct dirent * dp;
    struct stat s;
    while ((dp = readdir(dirp)) != nullptr)
    {
        word dname=dp->d_name;
        if (dname=="." || dname=="..") continue;
        if (stat(dname.c_str(),&s) == 0)
        {
            if (s.st_mode)
            {
                //Info<<"path is a file:"<<dname<<endl;
                word path=root+dname;
                if (!is_file(path.c_str())) continue;

                if (filter=="none")
                {
                    v.push_back(dname);
                }
                else
                {
                    label id=dname.find(filter);
                    if (id>=0)
                    {
                        v.push_back(dname);
                    }
                }
            }
        }
    }

    closedir(dirp);
    return;
}

void Foam::read_directory
(
    wordvec& v,
    const word& root,
    const word &filter
)
{
    DIR* dirp = opendir(root.c_str());
    struct dirent * dp;
    struct stat s;
    while ((dp = readdir(dirp)) != nullptr)
    {
        word dname=dp->d_name;

        if (dname=="." || dname=="..") continue;
        if (stat(dname.c_str(),&s) == 0)
        {
            if (s.st_mode & S_IFDIR)
            {

                if (filter=="none")
                {
                    v.push_back(dname);
                }
                else
                {
                    label id=dname.find(filter);
                    if (id>=0)
                    {
                        v.push_back(dname);
                    }
                }
             }
        }
    }
    closedir(dirp);

    return;
}

word foamSample::databaseDir_;
void foamSample::init()
{
    nx_=100;
    ny_=100;
    nz_=1;
    toDatabase_=false;
    visualize_=true;
    writeSamples_=true;
    bodySurfaceType_="bodySurface";

    mapBoundary_=true;

    mapCase_=cwd();
    sourceCase_=cwd();

    alphaMax_=0.5;
    byPhysicalType_=false;
    byname_=false;
    bytype_=false;
    byNameList_=true;
    nproc_=1;
    word homedir=getenv("HOME");
    databaseDir_=homedir+"/foamAI/database";
}

foamSample::foamSample()
:
foamMap()
{
    init();
}

foamSample::foamSample
(
    const Time* runTime,
    const fvMesh *_mesh
)
:
foamMap()
{
    runTime_=runTime;
    mesh=_mesh;
    init();
}

void foamSample::getMapTime()
{
}

void foamSample::combineFields()
{
}

void foamSample::constructKnn()
{
}

void foamSample::getBodySurface()
{
}


void foamSample::buildSearchTrees()
{

    source().construct_internal_knn(-0.1);

    farray2d trainData;

    forAll(mesh->boundaryMesh(),r)
    {
        word type=mesh->boundaryMesh()[r].type();
        if (type == "empty" || type == "symmetryPlane")
        {
            continue;
        }
        word ptype=mesh->boundaryMesh()[r].physicalType();
        word pname=mesh->boundaryMesh()[r].name();

        if (!byPhysicalType_)
        {
            ptype=mesh->boundaryMesh()[r].type();
        }

        //physical type as an indicator body surface
        if (bytype_ || byPhysicalType_)
        {
            if (ptype!=bodySurfaceType_)
            {
                continue;
            }
        }
        else if (byname_)
        {
            word key="none";
            if (bodyPatchNames_.size()>=1)
            {
                key=bodyPatchNames_[0];
            }
            if (pname!=key)
            {
                continue;
            }
        }
        else if (byNameList_)
        {
            bool found=false;
            for (unsigned ip=0;ip<bodyPatchNames_.size();ip++)
            {
                word filter=bodyPatchNames_[ip];
                if (pname==filter)
                {
                    found=true;
                    break;
                }
                label id=pname.find(filter);
                if (id>=0)
                {
                    found=true;
                    break;
                }
            }
            if (!found)
            {
                continue;
            }
        }

        forAll(mesh->boundaryMesh()[r],j)
        {
            point pt=mesh->Cf().boundaryField()[r][j];

            scalar xb=source().xbar(pt.x());
            scalar yb=source().ybar(pt.y());
            scalar zb=source().zbar(pt.z());
            std::vector<scalar> pts(3);
            pts[0]=xb;
            pts[1]=yb;
            pts[2]=zb;

            std::pair<label,label> faceid(r,j);
            bodyFaces_.append(faceid);
            trainData.push_back(pts);
        }
    }
    if (trainData.size()<2 && !parRun())
    {
        FatalErrorInFunction
            << "Error construct solid object search tree, zero number of solid\n"
            << "object boundary faces obtained. Please check whether bodySurfaceType name is set correctly."
            << exit(FatalError);
    }
    if (trainData.size()>=1)
    {
        Info<<"body surface faces:"<<trainData.size()<<endl;
        bodyBoundaryTree_.build(trainData);
    }
    else
    {
        bodyBoundaryTree_.deactivate();
    }

    return;
}


void foamSample::getGeomSample
(
    gridField &grid
)
{
    grid.setFieldSize();
    for (label i=0;i<grid.nx_+1;++i)
    {
        for (label j=0;j<grid.ny_+1;++j)
        {
            for (label k=0;k<grid.nzpt();++k)
            {
                point pt=grid.gridNodes_(i,j,k);

                scalar xb=source().xbar(pt.x());
                scalar yb=source().ybar(pt.y());
                scalar zb=source().zbar(pt.z());
                std::vector<scalar> qvect(3);
                qvect[0]=xb;
                qvect[1]=yb;
                qvect[2]=zb;
                std::vector<label> nb(2);
                std::vector<scalar> dist(2);
                source().internalMap_.search
                (
                    qvect,
                    nb,
                    dist
                );

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
                    insolid=0;
                }
                else
                {
                    if (dist1[0]>dist[0])
                    {
                        insolid=0;//nearer to fluid cell
                    }
                    else
                    {
                        // dist1[1]<=dist[0] can also be a fluid cell if cell_alpha>1.0e-3
                        label bfaceid=nb1[0];
                        label r=bodyFaces_[bfaceid].first;
                        label jf=bodyFaces_[bfaceid].second;
                        point pf=mesh->Cf().boundaryField()[r][jf];

                        label ic=nb[0];
                        point pc=mesh->C()[ic];
                        scalar delta=(pc-pt)&(pf-pt);
                        if (delta>0)
                        {
                            insolid=1;
                        }
                        else
                        {
                            insolid=0;//delta <0, in fluid
                        }

                    }
                }
                grid.insolids_(i,j,k)=insolid;

            }//for k
        }//for j
    }//for i

    return;
}

//grid map does not need wall distance

void foamSample::searchNodeField
(
    gridField &grid,
    const word &funct
)
{
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

                scalar xb=source().xbar(pt.x());
                scalar yb=source().ybar(pt.y());
                scalar zb=source().zbar(pt.z());
                std::vector<scalar> qvect(3);
                qvect[0]=xb;
                qvect[1]=yb;
                qvect[2]=zb;
                std::vector<label> nb(2);
                std::vector<scalar> dist(2);

                source().internalMap_.search
                (
                    qvect,
                    nb,
                    dist
                );

                //nb stores the nearest cell index
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
                    insolid=0;//cannot find body boundary match
                }
                else
                {
                    if (dist1[0]>dist[0])
                    {
                        insolid=0;//nearer to fluid cell
                        cnt1++;
                    }
                    else
                    {
                        // dist1[1]<=dist[0] can also be a fluid cell if cell_alpha>1.0e-3
                        label bfaceid=nb1[0];
                        label r=bodyFaces_[bfaceid].first;
                        label jf=bodyFaces_[bfaceid].second;
                        point pf=mesh->Cf().boundaryField()[r][jf];

                        label ic=nb[0];
                        point pc=mesh->C()[ic];
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
            }//k
        }//j
    }//i

    Info<<"number of solid node:"<<cnt<<" "<<cnt1<<" "<<cnt2<<endl;
    return;
}

void foamSample::getSampleData(const word &funct)
{
    std::map<word,gridField>::iterator it;


    for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
    {
        searchNodeField
        (
            it->second,
            funct
        );
    }
    return;
}

void foamSample::setOptions
(
    const argList& args
)
{
    foamMap::setOptions(args);

    args.optionReadIfPresent("nx", nx_);
    args.optionReadIfPresent("ny", ny_);
    args.optionReadIfPresent("nz", nz_);

    args.optionReadIfPresent("bodySurfaceType", bodySurfaceType_);
    args.optionReadIfPresent("bodySurfaceName", bodySurfaceType_);
}

void foamSample::getbodyPatchNames()
{
    forAll(mesh->boundaryMesh(),r)
    {
        word ptype=mesh->boundaryMesh()[r].physicalType();
        word pname=mesh->boundaryMesh()[r].name();

        if (!byPhysicalType_)
        {
            ptype=mesh->boundaryMesh()[r].type();
        }


        if (ptype==bodySurfaceType_)
        {
            bodyPatchNames_.push_back(pname);
        }
    }

}

void foamSample::writeVectorField
(
    gridField &grid,
    label k
)
{
    if (!visualize_)
    {
        return;
    }


    word casedir=cwd();

    word visudir=casedir+"/visualization/";
    word subdir=visudir+grid.name()+"_"+runTime_->timeName();
    word runtimdir=subdir+"/"+runTime_->timeName();
    for (unsigned iv=0;iv<grid.vectorSamples_.size();iv++)
    {
        word name=grid.vectorSamples_[iv];
        word fname=runtimdir+"/"+name;
        std::ofstream os(fname.c_str());
        if (!os)
        {
            FatalErrorInFunction
                << "Error write vector sample field, file: \n"<<fname<<" cannot open."
                << exit(FatalError);
        }

        scalar coef=source().scale(name);

        os<<"/*--------------------------------*- C++ -*----------------------------------*\\"<<std::endl;
        os<<"|       o        |                                                            |"<<std::endl;
        os<<"|    o     o     |  FOAM (R) : Open-source CFD for Enterprise                |"<<std::endl;
        os<<"|   o   O   o    |  Version : Dev                                             |"<<std::endl;
        os<<"|    o     o     |  ESI Ltd. <http://esi.com/>                            |"<<std::endl;
        os<<"|       o        |                                                            |"<<std::endl;
        os<<"\\*---------------------------------------------------------------------------*/"<<std::endl;
        os<<"FoamFile\n{\n";
        os<<"version     2.0;\n";
        os<<"format      ascii;\n";
        os<<"class       volVectorField;\n";
        os<<"object      "<<name<<";\n}\n";
        os<<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

        os<<"dimensions      [0 1 -1 0 0 0 0];\n\n\n";
        os<<"internalField   nonuniform List<vector>\n";
        os<<(grid.nx_+1)*(grid.ny_+1)*grid.nzpt()<<std::endl;
        os<<"(\n";

        DynamicList<point> uvw;
        uvw.setSize((grid.nx_+1)*(grid.ny_+1)*grid.nzpt());
        for (label i=0;i<grid.nx_+1;i++)
        {
            for (label j=0;j<grid.ny_+1;j++)
            {
                for (label s=0;s<grid.nzpt();s++)
                {
                    point v=grid.gridVectors_(i,j,s)[iv]*coef;
                    label index=grid.ijk(i,j,s);
                    uvw[index]=v;
                }
            }
        }

        for (label i=0;i<uvw.size();i++)
        {
                point v=uvw[i];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";

        }

        os<<")\n;\n";
        os<<"boundaryField\n{\n";
        os<<"    left\n    {\n";
        os<<"    type    inlet;\n";
        os<<"    value           nonuniform List<vector>\n";
        os<<(grid.ny_+1)*grid.nzpt()<<std::endl;
        os<<"(\n";
        for (label j=0;j<grid.ny_+1;j++)
        {
            for (label s=0;s<grid.nzpt();s++)
            {
                label idx=grid.ijk(0,j,s);
                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }
        }

        os<<")\n;\n}\n";

        os<<"    right\n    {\n";
        os<<"    type    outlet;\n";
        os<<"    value  nonuniform List<vector>\n";
        os<<(grid.ny_+1)*grid.nzpt()<<std::endl;
        os<<"(\n";
        for (label j=0;j<grid.ny_+1;j++)
        {
            for (label s=0;s<grid.nzpt();s++)
            {
                label idx=grid.ijk(grid.nx_,j,s);
                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }
        }
        os<<")\n;\n}\n";

        os<<"    top\n    {\n";
        os<<"    type    wall;\n";
        os<<"    value  nonuniform List<vector>\n";
        os<<(grid.nx_+1)*grid.nzpt()<<std::endl;
        os<<"(\n";
        for (label i=0;i<grid.nx_+1;i++)
        {
            for (label s=0;s<grid.nzpt();s++)
            {
                label idx=grid.ijk(i,0,s);
                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }
        }
        os<<")\n;\n}\n";

        os<<"    bottom\n    {\n";
        os<<"    type    wall;\n";
        os<<"    value   nonuniform List<vector>\n";
        os<<(grid.nx_+1)*grid.nzpt()<<std::endl;
        os<<"(\n";
        for (label i=0;i<grid.nx_+1;i++)
        {
            for (label s=0;s<grid.nzpt();s++)
            {
                label idx=grid.ijk(i,grid.ny_,s);

                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }
        }
        os<<")\n;\n}\n";

        os<<"    front\n{\n";
        os<<"    type    empty;\n}\n\n";

        os<<"    back\n{\n";
        os<<"    type    empty;\n}\n}\n";

        //
        os.close();

    }//iv
    return;

}

void foamSample::writeScalarField
(
    gridField &grid,
    label k
)
{
    if (!visualize_)
    {
        return;
    }

    word casedir=cwd();

    word visudir=casedir+"/visualization/";
    word subdir=visudir+grid.name()+"_"+runTime_->timeName();
    word runtimdir=subdir+"/"+runTime_->timeName();
    for (unsigned isc=0;isc<grid.scalarSamples_.size();isc++)
    {
        word name=grid.scalarSamples_[isc];
        word fname=runtimdir+"/"+name;
        std::ofstream os(fname.c_str());
        if (!os)
        {
            FatalErrorInFunction
                << "Error write scalar sample field, file: \n"<<fname<<" cannot open."
                << exit(FatalError);
        }

        scalar coef=source().scale(name);

        os<<"/*--------------------------------*- C++ -*----------------------------------*\\"<<std::endl;
        os<<"|       o        |                                                            |"<<std::endl;
        os<<"|    o     o     |  FOAM (R) : Open-source CFD for Enterprise                |"<<std::endl;
        os<<"|   o   O   o    |  Version : Dev                                             |"<<std::endl;
        os<<"|    o     o     |  ESI Ltd. <http://esi.com/>                            |"<<std::endl;
        os<<"|       o        |                                                            |"<<std::endl;
        os<<"\\*---------------------------------------------------------------------------*/"<<std::endl;
        os<<"FoamFile\n{\n";
        os<<"version     2.0;\n";
        os<<"format      ascii;\n";
        os<<"class       volScalarField;\n";
        os<<"object      "<<name<<";\n}\n";
        os<<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
        os<<"dimensions      [0 0 0 0 0 0 0];\n\n\n";
        os<<"internalField   nonuniform List<scalar>\n";
        os<<(grid.nx_+1)*(grid.ny_+1)*grid.nzpt()<<std::endl;
        os<<"(\n";
        DynamicList<scalar> scfield;
        scfield.setSize((grid.nx_+1)*(grid.ny_+1)*grid.nzpt());
        for (label i=0;i<grid.nx_+1;i++)
        {
            for (label j=0;j<grid.ny_+1;j++)
            {
                for (label s=0;s<grid.nzpt();s++)
                {
                    scalar d=grid.gridFields_(i,j,s)[isc]*coef;
                    label index=grid.ijk(i,j,s);
                    scfield[index]=d;
                }
            }
        }
        for (label i=0;i<scfield.size();i++)
        {
            os<<scfield[i]<<std::endl;
        }
        os<<")\n;\n";
        os<<"boundaryField\n{\n";
        os<<"    left\n    {\n";
        os<<"        type     zeroGradient;\n        }\n";
        os<<"    right\n    {\n";
        os<<"        type    zeroGradient;\n        }\n\n";

        os<<"    top\n    {\n";
        os<<"        type    zeroGradient;\n        }\n\n";

        os<<"    bottom\n    {\n";
        os<<"        type    zeroGradient;\n           }\n\n";

        os<<"    front\n{\n";
        os<<"    type    empty;\n}\n\n";

        os<<"    back\n{\n";
        os<<"    type    empty;\n}\n}\n";

        //
        os.close();
    }//isc
    return;
}

void foamSample::writeMeshDict
(
    gridField &grid
)
{
    if (!visualize_)
    {
        return;
    }

    scalar coef=1.0;
    scalar xmin=grid.xmin()-0.5*grid.dx_*coef;
    scalar xmax=grid.xmax()+0.5*grid.dx_*coef;

    scalar ymin=grid.ymin()-0.5*grid.dy_*coef;
    scalar ymax=grid.ymax()+0.5*grid.dy_*coef;

    scalar zmin=grid.zmin()-0.5*grid.dz_*coef;
    scalar zmax=grid.zmax()+0.5*grid.dz_*coef;

    word casedir=cwd();

    word visudir=casedir+"/visualization/";
    word subdir=visudir+grid.name()+"_"+runTime_->timeName();
    word sysdir=subdir+"/system";
    word gfile=sysdir+"/blockMeshDict";
    std::ofstream os(gfile.c_str());
    if (!os)
    {
        FatalErrorInFunction
        << "Error write blockDic file: \n"<<gfile<<" cannot open."
        << exit(FatalError);
    }

    os<<"/*--------------------------------*- C++ -*----------------------------------*\\"<<std::endl;
    os<<"|       o        |                                                            |"<<std::endl;
    os<<"|    o     o     |  FOAM (R) : Open-source CFD for Enterprise                |"<<std::endl;
    os<<"|   o   O   o    |  Version : Dev                                             |"<<std::endl;
    os<<"|    o     o     |  ESI Ltd. <http://esi.com/>                            |"<<std::endl;
    os<<"|       o        |                                                            |"<<std::endl;
    os<<"\\*---------------------------------------------------------------------------*/"<<std::endl;
    os<<"FoamFile\n{\n";
    os<<"version     2.0;\n";
    os<<"format      ascii;\n";
    os<<"class       dictionary;\n";
    os<<"object      blockMeshDict;\n}\n";
    os<<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    DynamicList<point> points;
    points.resize(8);
    points[0]=point(xmin,ymin,zmin);
    points[1]=point(xmax,ymin,zmin);
    points[2]=point(xmax,ymax,zmin);
    points[3]=point(xmin,ymax,zmin);
    points[4]=point(xmin,ymin,zmax);
    points[5]=point(xmax,ymin,zmax);
    points[6]=point(xmax,ymax,zmax);
    points[7]=point(xmin,ymax,zmax);
    //Info<<"number of points:"<<points.size()<<endl;
    os<<"convertToMeters 1;\n\n";
    os<<" vertices\n\(\n";
    for (label i=0;i<points.size();i++)
    {
        const point &pt=points[i];
        //Info<<"vrts:"<<pt<<endl;
        os<<"      ("<<pt.x()<<" "<<pt.y()<<" "<<pt.z()<<")\n";
    }
    os<<");\n\n";
    os<<"blocks\n\(\n";
    os<<"  hex (0 1 2 3 4 5 6 7) ("<<grid.nx_+1<<" "<<grid.ny_+1<<" "<<grid.nzpt()<<") simpleGrading (1 1 1)\n";
    os<<");\n"<<std::endl;
    os<<"edges\n(\n);\n";
    os<<"boundary\n\(\n\
    top\n\
    {\n\
        type patch;\n\
        faces\n\
        (\n\
        (3 7 6 2)\n\
        );\n\
    }\n\n";

    os<<"left\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
                (0 4 7 3  )\n\
            );\n\
        }\n\n";

    os<<"bottom\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
                (1 5 4 0)\n\
            );\n\
        }\n\n";

    os<<"right\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
                (2 6 5 1)\n\
            );\n\
        }\n\n";

    os<<"front\n\
        {\n\
            type empty;\n\
            faces\n\
            (\n\
                (4 5 6 7)\n\
            );\n\
        }\n\n";

    os<<"back\n\
        {\n\
        type empty;\n\
        faces\n\
        (\n\
          (0 3 2 1)\n\
        );\n\
        }\n\
\n\
    );\n";

    os.close();
    return;
}

void foamSample::writeVisualization()
{
    if (!visualize_)
    {
        Info<<"visualization option turned off, visualization files will not be written."<<endl;
        return;
    }

    Info<<"writing visualization file...."<<endl;

    word cmd="mkdir -p visualization";
    system(cmd.c_str());
    std::map<word,gridField>::iterator it;
    word casedir=cwd();
    rootcase_=casedir;
    word visudir=casedir+"/visualization/";
    for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
    {
        gridField &grid=it->second;
        word subdir=visudir+grid.name()+"_"+runTime_->timeName();
        cmd="mkdir -p "+subdir;
        system(cmd.c_str());

        word sysdir=subdir+"/system";
        cmd="mkdir -p "+sysdir;
        system(cmd.c_str());

        word constdir=subdir+"/constant";
        cmd="mkdir -p "+constdir;
        system(cmd.c_str());

        word meshdir=constdir+"/polyMesh";
        cmd="mkdir -p "+meshdir;
        system(cmd.c_str());

        word runtimdir=subdir+"/"+runTime_->timeName();
        cmd="mkdir -p "+runtimdir;
        system(cmd.c_str());
        writeMeshDict(it->second);
        cmd ="cp system/controlDict "+sysdir+"/";
        system(cmd.c_str());
        cmd ="cp system/fvSolution "+sysdir+"/";
        system(cmd.c_str());
        cmd ="cp system/fvSchemes "+sysdir+"/";
        system(cmd.c_str());

        writeScalarField(it->second);
        writeVectorField(it->second);
        chDir(subdir);
        cmd="blockMesh";
        system(cmd);
        chDir(rootcase_);
    }

    return;
}

void foamSample::writeNodeField()
{
    if (!writeSamples_)
    {
        Info<<"Write samples option set to false. Sample data will not be written into file."<<endl;
        return;
    }

    std::map<word,gridField>::iterator it;
    for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
    {
        writeNodeField(it->second);
    }

}


void foamSample::saveToDatabase(const word &funct)
{
    if (!toDatabase_)
    {
        Info<<"option to_database is false, data will not be saved to database"<<endl;
        return;
    }

    std::map<word,gridField>::iterator it;
    for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
    {
        saveToDatabase
        (
            it->second,
            funct
        );
    }
    return;
}

void foamSample::saveObjectMasks(const word &type)
{
    std::map<word,gridField>::iterator it;
    for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
    {
        saveObjectMask
        (
            it->second,
            type
        );
    }
    return;

}

void foamSample::saveObjectMask
(
    gridField &grid,
    const word &type
)
{

    word casedir=cwd();

    word dbsdir=databaseDir_+"/"+type;
    word dbsfile=dbsdir+"/objectMasks.dbs";
    databasePath_ = dbsfile;

    word fcaselist=dbsdir+"/caselist.dbs";
    std::vector<word> cases;
    std::ifstream is(fcaselist.c_str());
    if (is)
    {
        word str;
        while (!is.eof())
        {
            std::getline(is,str,'\n');
            if (str.length()<2)
            {
                continue;
            }
            cases.push_back(str);
        }
        is.close();
    }

    bool exists=false;
    if (cases.size()>=1)
    {
        std::vector<word>::iterator it;
        it=std::find (cases.begin(), cases.end(), casedir);
        if (it!=cases.end())
        {
            exists=true;
        }
    }
    if (exists)
    {
        Info<<"Warning: Current case already in database."<<endl;
        return;
    }

    std::ofstream os1(fcaselist.c_str(),std::ios::app);
    os1<<casedir<<std::endl;

    word cmd="mkdir -p "+dbsdir;
    system(cmd.c_str());
    std::ofstream os(dbsfile.c_str(),std::ios::app);
    word desc=type+","
             +casedir+","+grid.name()
             +","+std::to_string(grid.nx_+1)
             +","+std::to_string(grid.ny_+1)
             +","+std::to_string(grid.nzpt())
             +","+std::to_string(grid.xmin())
             +","+std::to_string(grid.ymin())
             +","+std::to_string(grid.zmin())
             +","+std::to_string(grid.xmax())
             +","+std::to_string(grid.ymax())
             +","+std::to_string(grid.zmax());
    os<<desc<<std::endl;
    word cnts="mask";
    for (label i=0;i<grid.nx_+1;++i)
    {
        for (label j=0;j<grid.ny_+1;++j)
        {
            for (label k=0;k<grid.nzpt();++k)
            {
                label insolid=grid.insolids_(i,j,k);
                cnts+=(","+std::to_string(insolid));
            }
        }
    }
    os<<cnts<<std::endl;
}

void foamSample::saveToDatabase
(
    gridField &grid,
    const word &funct
)
{
    if (parRun())
    {
        return;
    }

    word case_file_list=databaseDir_+"/caselist.dbs";
    //check whether case name already in database
    word bkname=caseName()+":"+grid.name();

    if (funct=="geometry")
    {
        //save case geometry
        Info<<"write geometric information into database."<<endl;
        word geofile=databaseDir_+"/geom.dbs";
        std::ofstream os(geofile.c_str(),std::ios::app);
        if (!os)
        {
            FatalErrorInFunction
                << "Error save sample geometry to database, database file \n"<<geofile<<" cannot open."
                << exit(FatalError);
        }

        os<<bkname<<std::endl;
        os<<grid.size()<<" "<<grid.nx_<<" "<<grid.ny_<<" "<<grid.nz_<<std::endl;
        std::vector<label> in_solid(grid.size());
        label cnt=0;
        for (label i=0;i<grid.nx_+1;++i)
        {
            for (label j=0;j<grid.ny_+1;++j)
            {
                for (label k=0;k<grid.nzpt();++k)
                {
                    label idx=grid.ijk(i,j,k);
                    in_solid[idx]=grid.insolids_(i,j,k);
                    cnt+=in_solid[idx];
                }
            }
        }
        Info<<"insolids:"<<cnt<<endl;
        for (unsigned i=0;i<in_solid.size();i++)
        {
            os<<in_solid[i]<<" ";
        }
        os<<std::endl;
        os.close();

        return;
    }//funct geometry


    std::ifstream is(case_file_list.c_str());
    std::map<word,word> saved_cases;
    if (is)
    {
        word str;
        while (!is.eof())
        {
            std::getline(is,str,'\n');
            if (str.length()<3) continue;
            std::vector<word> vtmp;
            stringtok(vtmp,str,"=");
            if (vtmp.size()!=2) continue;
            saved_cases[vtmp[0]]=vtmp[1];
        }
        is.close();
    }

    bool case_exists=true;
    if (saved_cases.size()==0)
    {
        case_exists=false;
    }
    else
    {
        std::map<word,word>::iterator it=saved_cases.find(bkname);
        if (it==saved_cases.end())
        {
            case_exists=false;
        }
    }

    if (!case_exists)
    {
        std::ofstream os(case_file_list.c_str(),std::ios::app);
        if (!os)
        {
            FatalErrorInFunction
            << "Error save sample data to database, database file \n"<<case_file_list<<" cannot open."
            << exit(FatalError);
        }

        os<<bkname<<"="<<source().rhoRef_<<" "<<source().Uref_.x()<<" "<<source().Uref_.y()<<" "
        <<source().Uref_.z()<<" "<<grid.xmin()<<" "<<grid.ymin()<<" "<<grid.zmin()<<" "
        <<grid.xmax()<<" "<<grid.ymax()<<" "<<grid.zmax()<<std::endl;
        os.close();
    }


    //write field data
    word caseDatadir=databaseDir_+"/caseData/";
    word cmd="mkdir -p "+caseDatadir;
    system(cmd.c_str());

    word datafile=caseDatadir+bkname+".dbs";
    std::ofstream os(datafile.c_str(),std::ios::out);
    if (!os)
    {
        FatalErrorInFunction
            << "Error save sample data to database, database file \n"<<datafile<<" cannot open."
            << exit(FatalError);
    }

    os<<grid.nx_<<" "<<grid.ny_<<" "<<grid.nz_<<" "<<grid.containSolidObject_<<" "
    <<grid.xmin()<<" "<<grid.ymin()<<" "<<grid.zmin()<<" "
    <<grid.xmax()<<" "<<grid.ymax()<<" "<<grid.zmax()<<std::endl;

    word header="insolid ";
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        //const word& scname = iter();
        header+=(iter() + " ");
    }
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& vname = iter();
        word vx=vname+"_x";
        word vy=vname+"_y";
        word vz=vname+"_z";
        header+=(" "+vx+" "+vy+" "+vz);
    }

    os<<header<<std::endl;

    DynamicList<std::vector<scalar>> data;
    label sz=grid.size();
    data.resize(sz);
    label ndata=1+mapScalarFields_.size()+3*mapVectorFields_.size();
    for (label i=0;i<grid.nx_+1;i++)
    {
        for (label j=0;j<grid.ny_+1;j++)
        {
            for (label k=0;k<grid.nzpt();k++)
            {
                label idx=grid.ijk(i,j,k);
                std::vector<scalar> vect(ndata);
                vect[0]=grid.insolids_(i,j,k);

                for (label isc=0;isc<mapScalarFields_.size();isc++)
                {
                    scalar d=grid.gridFields_(i,j,k)[isc];
                    vect[isc+1]=d;
                }//isc

                if (mapVectorFields_.size()>=1)
                {
                    for (label iv=0;iv<mapVectorFields_.size();iv++)
                    {
                        point vd=grid.gridVectors_(i,j,k)[iv];
                        label idata=mapScalarFields_.size()+3*iv+1;
                        vect[idata]=vd.x();
                        vect[idata+1]=vd.y();
                        vect[idata+2]=vd.z();
                    }
                }//if
                data[idx]=vect;
            }//k
        }//j
    }//i
    os<<data.size()<<std::endl;
    for (label ik=0;ik<data.size();ik++)
    {
        for (label id=0;id<ndata;id++)
        {
            os<<std::setprecision(4)<<data[ik][id];
            if (id<ndata-1)
            {
                os<<" ";
            }
            else
            {
                os<<std::endl;
            }
        }
    }

    os.close();
    return;
}

void foamSample::writeNodeField
(
    gridField &grid
)
{

    Info<<"mapScalarFields:"<<mapScalarFields_.size()<<endl;


    word dir=cwd()+"/postProcessing/";

    word filename=dir+grid.name()+"_"+runTime_->timeName();
    Info<<filename<<endl;
    word cmd="mkdir -p "+dir;
    system(cmd.c_str());
    label sz=grid.size();
    //write the distance file

    word dfile=filename+"_dists";
    std::ofstream osd(dfile.c_str());

    std::vector<word> ddata(sz);
    if (osd)
    {
        word ss;
        for (label i=0;i<grid.nx_+1;i++)
        {
            for (label j=0;j<grid.ny_+1;j++)
            {
                for (label k=0;k<grid.nzpt();k++)
                {
                    label idx=grid.ijk(i,j,k);
                    ddata[idx]=dbl2str(grid.dists_(i,j,k));
                }
            }
        }

        for (label s=0;s<sz;s++)
        {
            ss+=(ddata[s]+"\n");
        }

        osd<<ss;

    }

    std::ofstream os(filename.c_str());
    if (!os)
    {
        FatalErrorInFunction
            << "Error save sample data, sample file \n"<<filename<<" cannot open."
            << exit(FatalError);
    }
    word header="x y z";
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        //const word& scname = iter();
        header+=(" " + iter());
    }
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& vname = iter();
        word vx=vname+"_x";
        word vy=vname+"_y";
        word vz=vname+"_z";
        header+=(" "+vx+" "+vy+" "+vz);
    }

    header+=" in_solid_object";
    Info<<header<<endl;
    os<<header<<std::endl;
    DynamicList<std::vector<scalar>> data;
    //label sz=grid.size();
    data.resize(sz);
    label ndata=3+mapScalarFields_.size()+1+3*mapVectorFields_.size();
    for (label i=0;i<grid.nx_+1;i++)
    {
        for (label j=0;j<grid.ny_+1;j++)
        {
            for (label k=0;k<grid.nzpt();k++)
            {
                label idx=grid.ijk(i,j,k);
                std::vector<scalar> vect(ndata);

                vect[0]=grid.gridNodes_(i,j,k).x();
                vect[1]=grid.gridNodes_(i,j,k).y();
                vect[2]=grid.gridNodes_(i,j,k).z();

                forAllConstIter(HashSet<word>, mapScalarFields_, iter)
                {
                    const word& name = iter();
                    label isc = grid.scalarId(name);
                    scalar coef=source().scale(name);
                    scalar d=grid.gridFields_(i,j,k)[isc];
                    vect[3+isc]=d*coef;
                }//isc
                forAllConstIter(HashSet<word>, mapVectorFields_, iter)
                {
                    const word& name = iter();
                    label iv = grid.vectorId(name);
                    point vd=grid.gridVectors_(i,j,k)[iv];
                    scalar coef=source().scale(name);
                    vd*=coef;
                    label idata=3+mapScalarFields_.size()+3*iv;
                    vect[idata]=vd.x();
                    vect[idata+1]=vd.y();
                    vect[idata+2]=vd.z();
                }
                vect[ndata-1]=grid.insolids_(i,j,k);
                data[idx]=vect;
            }//k
        }//j
    }//i

    for (label ik=0;ik<data.size();ik++)
    {
        for (label id=0;id<ndata;id++)
        {
            os<<std::setprecision(4)<<data[ik][id];
            if (id<ndata-1)
            {
                os<<" ";
            }
            else
            {
                os<<std::endl;
            }
        }
    }

    os.close();
    return;
}



void foamSample::setInput
(
    const dictionary& dict,
    const argList& args
)
{
    foamMap::setInput(dict, args);

    getFieldTypes(*mesh, runTime_->timeName());

    visualize_=dict.lookupOrDefault<bool> ("visualize",false);

    bodySurfaceType_=dict.lookupOrDefault<word> ("bodySurfaceType","car");

    //byType byNameList byPhysicalType byName
    toDatabase_=dict.lookupOrDefault<bool> ("to_database",false);

    word identifyObject=dict.lookupOrDefault<word> ("identifyObject","byNameList");

    if (identifyObject=="byType")
    {
        byPhysicalType_=false;
        bytype_=true;
        byname_=false;
        byNameList_=false;
    }

    else if (identifyObject=="byPhysicalType")
    {
        byPhysicalType_=true;
        bytype_=false;
        byname_=false;
        byNameList_=false;
    }

    else if (identifyObject=="byName")
    {
        byPhysicalType_=false;
        bytype_=false;
        byname_=true;
        byNameList_=false;
    }
    else if (identifyObject=="byNameList")
    {
        byPhysicalType_=false;
        bytype_=false;
        byname_=false;
        byNameList_=true;
    }

    word blockName=dict.lookupOrDefault<word> ("blockName","block1");

    point p0(0,0,0);
    point p1(1,1,1);
    boundBox bbox0(p0,p1);
    boundBox bbox=dict.lookupOrDefault<boundBox> ("boundBox",bbox0);
    point dp=bbox.max()-bbox.min();
    scalar coef=0.05;
    point pmin=bbox.min()-coef*dp;
    point pmax=bbox.max()+coef*dp;
    samplebox_=new
    boundBox
    (
        pmin,
        pmax
    );

    vector ndivs0(127,127,1);
    vector ndivisions=dict.lookupOrDefault<vector> ("ndivisions",ndivs0);

    label nx=label(ndivisions[0]);
    label ny=label(ndivisions[1]);
    label nz=label(ndivisions[2]);
    bool contain_solid_object=true;
    writeSamples_=dict.lookupOrDefault<bool> ("write_samples",true);
    bool interp=dict.lookupOrDefault<bool> ("interpolation",false);

    gridField grid
    (
        bbox,
        nx,
        ny,
        nz,
        contain_solid_object
    );
    grid.interp_=interp;
    grid.name_=blockName;

    grid.setSampleNames
    (
        mapScalarFields_,
        mapVectorFields_
    );
    grid.generateNodes();
    grid.setFieldSize();
    sampleGrids_[blockName]=grid;

    const polyBoundaryMesh& bMesh = mesh->boundaryMesh();

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
        bodyPatchNames_.push_back
        (
            allPatchNames[*it]
        );
    }

    return;
}

