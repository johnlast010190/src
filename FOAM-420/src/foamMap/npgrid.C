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
    Foam::npGrid

Group
    grpFoamMap

Description
   A class designed to convert Numpy array into FOAM field for visualization
   in Paraview

SourceFile
    npgrid.C

\*---------------------------------------------------------------------------*/

#include "npgrid.H"
#include "foamMap.H"

using namespace Foam;
npGrid::npGrid()
:
gridField()
{;}

//read content of a text file
void npGrid::read_text_file
(
    word &input,
    const word &fname
)
{
    std::ifstream is(fname.c_str());

    if (is)
    {
        std::ostringstream ss;
        ss<<is.rdbuf();
        input=ss.str();
    }
    return;
}

npGrid::npGrid
(
    const word &blkdef,
    const word &npfiles
)
:
gridField()
{
    std::vector<word> vtmp;
    std::ifstream is(blkdef.c_str());
    if (!is)
    {
        FatalErrorInFunction<<"Error read np field, block definition file "<<blkdef<<" not found.\n"
            << exit(FatalError);
    }
    word s1;
    word line;
    while (!is.eof())
    {
        std::getline(is,s1,'\n');
        if (s1.length()<2) continue;
        if (s1[0]=='#') continue;
        line=s1;
        break;
    }
    foamMap::stringtok(vtmp,line," ");
    if (vtmp.size()!=9)
    {
        FatalErrorInFunction<<"Error read np field, incorrect block definition: "<<line<<"\n"
            << exit(FatalError);
    }
    //128 128 1 -0.5 -1 0.4 1.5 1 0.6
    nx_=atoi(vtmp[0].c_str())-1;
    ny_=atoi(vtmp[1].c_str())-1;
    nz_=atoi(vtmp[2].c_str())-1;
    if (nz_<1)
    {
        nz_=1;
    }
    xmin_=atof(vtmp[3].c_str());
    ymin_=atof(vtmp[4].c_str());
    zmin_=atof(vtmp[5].c_str());
    xmax_=atof(vtmp[6].c_str());
    ymax_=atof(vtmp[7].c_str());
    zmax_=atof(vtmp[8].c_str());

    std::vector<word> vfiles;
    foamMap::stringtok(vfiles,npfiles,",");


    for (unsigned i=0;i<vfiles.size();i++)
    {
        word fname=vfiles[i];
        foamMap::stringtok(vtmp,fname,":");
        if (vtmp.size()==1 && fname!=blkdef)
        {
            scalarSamples_.push_back(fname);
        }
        else if (vtmp.size()==2)
        {
            if (vectorSamples_.empty())
            {
                vectorSamples_.push_back(vtmp[0]);
            }
            else if (inlist(vtmp[0],vectorSamples_)<0)
            {
                vectorSamples_.push_back(vtmp[0]);
            }
        }
    }

    //blockdef:nxp nyp nzp xmin ymin zmin xmax ymax zmax

    for (unsigned f=0;f<scalarSamples_.size();f++)
    {
        word fname=scalarSamples_[f];
        DynamicList<scalar> scfield;
        scfield.setSize(nnodes());
        getScalarField(fname,scfield);
        scalarFields_[fname]=scfield;
    }


    for (unsigned i=0;i<vfiles.size();i++)
    {
        word fname=vfiles[i];
        foamMap::stringtok(vtmp,fname,":");
        if (vtmp.size()==2)
        {
            DynamicList<scalar> scfield;
            scfield.setSize(nnodes());
            getScalarField(fname,scfield);
            vectorComps_[fname]=scfield;
        }
    }

    combineVectorField();
}

void npGrid::combineVectorField()
{
    DynamicList<point> vecfield;
    vecfield.setSize(nnodes());
    for (unsigned ic=0;ic<vectorSamples_.size();ic++)
    {
        word vecnme=vectorSamples_[ic];
        word vnamex=vecnme+":0";
        word vnamey=vecnme+":1";
        std::map<word,nodeField>::iterator it=vectorComps_.find(vnamex);
        if (it==vectorComps_.end())
        {
            FatalErrorInFunction
                << "Error construct vector field, x-component field not found.\n"
                << exit(FatalError);
        }

        std::map<word,nodeField>::iterator it1=vectorComps_.find(vnamey);
        if (it1==vectorComps_.end())
        {
            FatalErrorInFunction
                << "Error construct vector field, y-component field not found.\n"
                << exit(FatalError);
        }

        const nodeField &Ux=it->second;
        const nodeField &Uy=it1->second;
        label dim=2;
        word vnamez=vecnme+":2";
        std::map<word,nodeField>::iterator it2=vectorComps_.find(vnamez);
        if (it2!=vectorComps_.end())
        {
            dim=3;
        }

        for (label ip=0;ip<vecfield.size();ip++)
        {
            point pt;
            pt[0]=Ux[ip];
            pt[1]=Uy[ip];
            pt[2]=0;
            vecfield[ip]=pt;
        }

        if (dim==3)
        {
            const nodeField &Uz=it2->second;
            for (label ip=0;ip<vecfield.size();ip++)
            {
                vecfield[ip][2]=Uz[ip];
            }
        }
        vectorFields_[vecnme]=vecfield;
    }//for ic

    Info<<"Number of vector fields constructed:"<<vectorFields_.size()<<endl;

    return;
}

void npGrid::writeMeshDict()
{
    scalar dx=(xmax()-xmin())/scalar(nx_);
    scalar dy=(ymax()-ymin())/scalar(ny_);
    scalar dz=(zmax()-zmin())/scalar(nz_);
    scalar xmin0=xmin()-0.5*dx;
    scalar xmax0=xmax()+0.5*dx;

    scalar ymin0=ymin()-0.5*dy;
    scalar ymax0=ymax()+0.5*dy;

    scalar zmin0=zmin()-0.5*dz;
    scalar zmax0=zmax()+0.5*dz;

    word casedir=cwd();

    word sysdir=casedir+"/system";

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
    points[0]=point(xmin0,ymin0,zmin0);
    points[1]=point(xmax0,ymin0,zmin0);
    points[2]=point(xmax0,ymax0,zmin0);
    points[3]=point(xmin0,ymax0,zmin0);
    points[4]=point(xmin0,ymin0,zmax0);
    points[5]=point(xmax0,ymin0,zmax0);
    points[6]=point(xmax0,ymax0,zmax0);
    points[7]=point(xmin0,ymax0,zmax0);
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
    os<<"  hex (0 1 2 3 4 5 6 7) ("<<nx_+1<<" "<<ny_+1<<" "<<nzpt()<<" )simpleGrading (1 1 1)\n";

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

void npGrid::writeVisualization()
{

    word casedir=cwd();

    word visudir=casedir;
    word sysdir=visudir+"/system";
    word cmd="mkdir -p "+sysdir;
    system(cmd.c_str());

    word constdir=casedir+"/constant";
    cmd="mkdir -p "+constdir;
    system(cmd.c_str());
    word meshdir=constdir+"/polyMesh";
    cmd="mkdir -p "+meshdir;
    system(cmd.c_str());

    word runtimdir=casedir+"/"+timeName();
    cmd="mkdir -p "+runtimdir;
    system(cmd.c_str());
    writeMeshDict();

    writeScalarField();
    writeVectorField();

    cmd="blockMesh";
    system(cmd.c_str());
    return;
}

void npGrid::writeVectorField
(
    label k
)
{

    word subdir=cwd();

    word runtimdir=subdir+"/1";
    word cmd="mkdir -p "+runtimdir;
    system(cmd.c_str());

    for (unsigned iv=0;iv<vectorSamples_.size();iv++)
    {
        word name=vectorSamples_[iv];
        word fname=runtimdir+"/"+name;
        std::ofstream os(fname.c_str());
        if (!os)
        {
            FatalErrorInFunction
                << "Error write vector sample field, file: \n"<<fname<<" cannot open."
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
        os<<"class       volVectorField;\n";
        os<<"object      "<<name<<";\n}\n";
        os<<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

        os<<"dimensions      [0 1 -1 0 0 0 0];\n\n\n";
        os<<"internalField   nonuniform List<vector>\n";

        os<<(nx_+1)*(ny_+1)*nzpt()<<std::endl;

        os<<"(\n";
        const DynamicList<vector> &uvw=vectorFields_[name];


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
        os<<(ny_+1)*nzpt()<<std::endl;
        os<<"(\n";
        for (label j=0;j<ny_+1;j++)
        {
            for (label s=0;s<nzpt();s++)
            {
                label idx=ijk(0,j,s);
                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }
        }

        os<<")\n;\n}\n";

        os<<"    right\n    {\n";
        os<<"    type    outlet;\n";
        os<<"    value  nonuniform List<vector>\n";
        os<<(ny_+1)*nzpt()<<std::endl;
        os<<"(\n";
        for (label j=0;j<ny_+1;j++)
        {
            for (label s=0;s<nzpt();s++)
            {
                label idx=ijk(nx_,j,s);
                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }

        }
        os<<")\n;\n}\n";

        os<<"    top\n    {\n";
        os<<"    type    wall;\n";
        os<<"    value  nonuniform List<vector>\n";
        os<<(nx_+1)*nzpt()<<std::endl;
        os<<"(\n";
        for (label i=0;i<nx_+1;i++)
        {
            for (label s=0;s<nzpt();s++)
            {
                label idx=ijk(i,0,s);
                point v=uvw[idx];
                os<<"("<<v.x()<<" "<<v.y()<<" "<<v.z()<<")\n";
            }
        }
        os<<")\n;\n}\n";

        os<<"    bottom\n    {\n";
        os<<"    type    wall;\n";
        os<<"    value   nonuniform List<vector>\n";
        os<<(nx_+1)*nzpt()<<std::endl;
        os<<"(\n";
        for (label i=0;i<nx_+1;i++)
        {
            for (label s=0;s<nzpt();s++)
            {
                label idx=ijk(i,ny_,s);
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

void npGrid::writeScalarField
(
    label k
)
{

    word subdir=cwd();


    word runtimdir=subdir+"/1";
    word cmd="mkdir -p "+runtimdir;
    system(cmd.c_str());

    for (unsigned isc=0;isc<scalarSamples_.size();isc++)
    {
        word name=scalarSamples_[isc];
        word fname=runtimdir+"/"+name;
        std::ofstream os(fname.c_str());
        if (!os)
        {
            FatalErrorInFunction
                << "Error write scalar sample field, file: \n"<<fname<<" cannot open."
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
        os<<"class       volScalarField;\n";
        os<<"object      "<<name<<";\n}\n";
        os<<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

        os<<"dimensions      [0 0 0 0 0 0 0];\n\n\n";
        os<<"internalField   nonuniform List<scalar>\n";

        os<<(nx_+1)*(ny_+1)*nzpt()<<std::endl;

        os<<"(\n";
        const DynamicList<scalar> &scfield=scalarFields_[name];

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

void npGrid::getScalarField
(
    const word& fname,
    DynamicList<scalar> &scfield
)
{
    word input;
    read_text_file(input,fname);
    std::vector<word> vtmp;
    foamMap::stringtok
    (
        vtmp,
        input,
        "\n"
    );
    DynamicList<DynamicList<scalar>>vdata;
    for (unsigned i=0;i<vtmp.size();i++)
    {
        word str=vtmp[i];
        if (str.size()<3)
        {
            continue;
        }

        foamMap::trimleft(str,"[");
        foamMap::trimright(str,"]");
        std::vector<word> vtmp1;
        foamMap::stringtok
        (
            vtmp1,
            str,
            ","
        );
        DynamicList<scalar> ddat;
        ddat.setSize(vtmp1.size());
        for (unsigned k=0;k<vtmp1.size();k++)
        {
            ddat[k]=atof(vtmp1[k].c_str());
        }
        vdata.append(ddat);
    }

    for (label i=0;i<vdata.size();i++)
    {
        for (label j=0;j<vdata[i].size();j++)
        {
            label idx=ijk(i,j,0);
            scfield[idx]=vdata[i][j];
        }
    }

    Info<<"number of data items for scalar "<<fname<<":"<<vdata.size()<<endl;
    return;
}
