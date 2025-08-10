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
    Foam::foamSample

Group
    grpfoamMap

Description
    A set of classes designed to sample an arbitrary number of flow fields from an unstructured
    mesh region into a uniform grid, or another unstructured mesh region.

SourceFile
    blockSample.H

\*---------------------------------------------------------------------------*/


#pragma once

#include <unistd.h>
#include "foamMap.H"
#include "array3d.H"
#include "gridField.H"
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

namespace Foam
{
typedef std::vector<word> wordvec;

inline bool is_file
(
    const char* path
)
{
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

inline bool is_dir
(
    const char* path
)
{
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

void get_files
(
    wordvec& v,
    const word& root="./",
    const word &contains="none"
);

void read_directory
(
    wordvec& v,
    const word& root="./",
    const word &filter="none"
);



class foamSample
:
    public foamMap
{
    protected:
        autoPtr<boundBox> samplebox_;
        word databasePath_;
    public:
        KNN bodyBoundaryTree_; //boundary map for the bodies
        // KNN interfaceTree_; //process interface search

        boundBox &sampleBox()
        {
            return samplebox_();
        }

        const fvMesh *mesh;
        const Time* runTime_;
        std::map<word,gridField> sampleGrids_;
        static word databaseDir_;
        bool toDatabase_;
        word rootcase_;
        bool visualize_;
        bool writeSamples_;
        label nproc_;
        DynamicList<std::pair<label,label>> bodyFaces_;
        DynamicList<std::pair<label,label>> procIfaces_;

        word bodySurfaceType_;

        bool byPhysicalType_;
        bool bytype_;
        bool byname_;
        bool byNameList_;
        label nx_;
        label ny_;
        label nz_;
        std::vector<word> bodyPatchNames_;

        virtual void setOptions
        (
            const argList& args
        );

        bool parRun()
        {
            if (nproc_>1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool isProcInterface(const word &name)
        {
            word wd="procBoundary";
            label j=name.find(wd);
            if (j>=0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        virtual void getMapTime();
        virtual void combineFields();
        virtual void constructKnn();
        virtual void getBodySurface();

        void getbodyPatchNames();

        virtual void setInput
        (
            const dictionary& dict,
            const argList& args
        );


        word databasePath() const
        {
            return databasePath_;
        }

        void getSampleData(const word &funct="field");

        void buildSearchTrees();

        virtual void searchNodeField
        (
            gridField &grid,
            const word &funct="field"
        );

        void getGeomSample
        (
            gridField &grid
        );

        void writeNodeField
        (
             gridField &grid
        );

        void saveObjectMask(gridField &grid,const word &type);
        void saveObjectMasks(const word &type);

        void saveToDatabase(const word &funct="field"); //save block samples to database
        void saveToDatabase
        (
            gridField &grid,
            const word &funct="field"
        );

        void writeNodeField();
        void writeVisualization();
        void writeMeshDict
        (
            gridField &grid
        );

        void writeScalarField
        (
            gridField &grid,
            label k=0
        );

        void writeVectorField
        (
            gridField &grid,
            label k=0
        );

        inline void getGridField
        (
            gridField &grid,
            label i,
            label j,
            label k,
            const std::vector<label> &nb,
            const std::vector<scalar> &dist,
            bool insolid
        );

        inline void outDomainGrid
        (
            gridField &grid,
            label i,
            label j,
            label k
        );

        inline void getGridVectorField
        (
            gridField &grid,
            label i,
            label j,
            label k,
            const std::vector<label> &nb,
            const std::vector<scalar> &dist,
            bool insolid
        );

        inline scalar interp
        (
            const word &name,
            const std::vector<label> & nb,
            const std::vector<scalar> &dist
        );

        inline point vinterp
        (
            const word &name,
            const std::vector<label> & nb,
            const std::vector<scalar> &dist
        );

        void init();

        foamSample();

        foamSample
        (
            const Time* runTime,
            const fvMesh* mesh
        );
        virtual ~foamSample(){;}

};

scalar foamSample::interp
(
    const word& name,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    scalar sum1=0;
    scalar sum2=0;
    scalar eps=1.0e-25;
    std::vector<scalar> prop(nb.size());
    for (unsigned i=0;i<nb.size();i++)
    {
        label ic=nb[i];
        if (!Pstream::parRun())
        {
            prop[i]=source().scalarfields_[name][ic];
        }
        else
        {
            prop[i]=sourceScalarFields_[name][ic];
        }
    }

    for (unsigned i=0; i<nb.size();i++)
    {
        scalar disti=dist[i]+eps;
        sum1+=prop[i]/disti;
        sum2+=1.0/disti;
    }

    return sum1/sum2;
}

point foamSample::vinterp
(
    const word& name,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    point sum1(0,0,0);
    scalar sum2=0;
    scalar eps=1.0e-25;
    std::vector<point> prop(nb.size());

    for (unsigned i=0;i<nb.size();i++)
    {
        label ic=nb[i];
        if (!Pstream::parRun())
        {
            prop[i]=source().vectorfields_[name][ic];
        }
        else
        {
            prop[i]=sourceVectorFields_[name][ic];
        }
    }

    for (unsigned i=0; i<nb.size();i++)
    {
        scalar disti=dist[i]+eps;
        sum1+=prop[i]/disti;
        sum2+=1.0/disti;
    }
    return sum1/sum2;
}

void foamSample::outDomainGrid
(
    gridField& grid,
    label i,
    label j,
    label k
)
{
    label nsc=grid.scalarSamples_.size();
    if (nsc>0)
    {
        nodeField node(nsc,1.0e20);
        grid.gridFields_(i,j,k)=node;
    }

    label nvc=grid.vectorSamples_.size();

    if (nvc>0)
    {
        point bigVector(1.0e20,1.0e20,1.0e20);
        nodeVectors vnode(nvc,bigVector);
        grid.gridVectors_(i,j,k)=vnode;
    }
    return;
}

//get the scalar field at a nodal point by kdTree search
void foamSample::getGridField
(
    gridField& grid,
    label i,
    label j,
    label k,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist,
    bool insolid
)
{
    label nsc=grid.scalarSamples_.size();

    grid.dists_(i,j,k)=dist[0];
    if (nsc>0)
    {
        nodeField node(nsc,0.0);
        if (insolid)
        {
            grid.gridFields_(i,j,k)=node;
        }
        else
        {
            for (label isc=0;isc<nsc;isc++)
            {
                word name=grid.scalarSamples_[isc];
                label ic=nb[0];
                if (!grid.interp_)
                {
                    scalar sval;
                    if (!Pstream::parRun())
                    {
                        sval=source().scalarfields_[name][ic];
                    }
                    else
                    {
                        sval=sourceScalarFields_[name][ic];
                    }
                    node[isc]=sval;
                }
                else
                {
                    scalar sval=interp(name,nb,dist);
                    node[isc]=sval;
                }
            }

            grid.gridFields_(i,j,k)=node;
        }

    }//nsc>0

    return;
}

//get the vector field at a nodal point
void foamSample::getGridVectorField
(
    gridField& grid,
    label i,
    label j,
    label k,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist,
    bool insolid
)
{
    label nvc=grid.vectorSamples_.size();

    if (nvc>0)
    {
        point zeroVector(0,0,0);
        nodeVectors node(nvc,zeroVector);

        if (insolid)
        {
            grid.gridVectors_(i,j,k)=node;
        }
        else
        {
            for (label ivc=0;ivc<nvc;ivc++)
            {
                word name=grid.vectorSamples_[ivc];
                label ic=nb[0];
                point svec;
                if (!grid.interp_)
                {
                    if (!Pstream::parRun())
                    {
                        svec=source().vectorfields_[name][ic];
                    }
                    else
                    {
                        svec=sourceVectorFields_[name][ic];
                    }

                    node[ivc]=svec;
                }
                else
                {
                    point svec=vinterp(name,nb,dist);
                    node[ivc]=svec;
                }
            }
            grid.gridVectors_(i,j,k)=node;
        }

    }//nvc>0

    return;
}


}//Foam
