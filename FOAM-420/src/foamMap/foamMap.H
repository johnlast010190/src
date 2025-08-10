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
    grpFoamMap

Description
    A set of classes designed to map an arbitrary number of flow fields from an unstructured
    mesh region into a uniform grid, or another unstructured mesh region.

SourceFile
    FoamMap.C

\*---------------------------------------------------------------------------*/
#pragma once
#include "array3d.H"
#include "readFields.H"
#include <sstream>
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/Time/Time.H"
#include <set>
#include "smBoundBox.H"


namespace Foam
{

template<class T>
void interProcTrans
(
    DynamicList<T>& field,
    const DynamicList<point>& verts,
    const List<boundBox>& tboxs,
    const List<boundBox>& sboxs,
    const DynamicList<T>& source,
    const word& name
);

class foamMap
{
    protected:
        bool norm_;
        label parmodel_;
        scalar extcoef_;
        labelHashSet wallDistPatchs_;
        wordReList wallPatchNames_;
        word bboxType_;
        label boxCellNum_;
        bool reportMem_; //report memory usage for debug
        bool scaleSource_;
    public:
        static label trimleft(word& str,const word& del = " ");
        static label trimright(word& str, const word& del = " ");
        static void trim(word& str, const word& del = " ");
        static void replace_str(word& str);

        typedef Tuple2<label,label> index;

        word casedir_;
        word caseType_;
        const fvMesh* mesh;

        label nprocs_;
        label myid_;
        label masterid_;
        scalar deltaBox_;
        DynamicList<label> procs_;
        Array3d <label> gblocks_;

        bool wdistMap_;
        DynamicList<bool> cellInDomain_;

        DynamicList<point> sourceXyz_;
        std::map<word,DynamicList<scalar>> sourceScalarFields_;
        std::map<word,DynamicList<point>> sourceVectorFields_;
        DynamicList<scalar> sourceWdists_;
        std::map<label,label> cellMap_;
        std::vector<KNN> kdTrees_;
        KNN internalMap_;
        scalar deltaY_;

        scalar sourceAlphaMax_;
        Array2d<label> alphaMap_;//parallel run alphaMap

        scalar xmin_;
        scalar xmax_;
        scalar ymin_;
        scalar ymax_;
        scalar zmin_;
        scalar zmax_;

        DynamicList<label> toProcs_;
        DynamicList<label> fromProcs_;

        DynamicList<index> interProcs_;
        scalar alphaMax_;
        scalar errorBound_;

        scalar rhoRefSource_;
        scalar rhoRefTarget_;
        scalar UrotDegreeFromSource_; //rotating degree frpm source to target, ani-clockwise positive
        vector UrefSource_;
        vector UrefTarget_;
        bool interp_;
        bool mapBoundary_;
        bool mapFixedBC_;
        word function_;
        word mapTimeName_;
        word tgtTimeName_;
        HashTable<word> fieldTypes_;

        word sourceCase_;
        HashSet<word> mapScalarFields_;
        HashSet<word> mapVectorFields_;
        wordList mapFieldNames_;

        HashSet<word> mapSurfaceScalarFields_;
        label nwdist_;
        word mapCase_;

        wordList sourceRegions_;
        wordList targetRegions_;
        wordList targetRegionDirs_;
        bool allRegions_;

        autoPtr<class readFields> source_;
        autoPtr<class readFields> target_;

        const labelHashSet& wallDistPatchs()
        {
            return wallDistPatchs_;
        }

        bool reportMemUsage()
        {
            return reportMem_;
        }

        bool scaleSource()
        {
            return scaleSource_;
        }

        void reportMemory
        (
            memInfo& mem,
            const word& desc,
            scalar& maxUsage
        );


        void getWallDistPatchs
        (
            const fvMesh& mesh
        );


        label parModel()
        {
            return parmodel_;
        }

        scalar extcoef()
        {
            return extcoef_;
        }

        bool inlist
        (
            label k,
            const DynamicList<label>& list
        );

        bool inlist
        (
            const word& key,
            const std::vector<word>& list
        );


        virtual void buildpKdTrees();//build parallel run kdTrees

        virtual void construct_internal_knn
        (
            scalar cutoff = 0.8
        );


        scalar xbar(scalar x)
        {
            return target().xbar(x);
        }

        scalar ybar(scalar y)
        {
            return target().ybar(y);
        }

        scalar zbar(scalar z)
        {
            return target().zbar(z);
        }


        label inGroup(label i)
        {
            if (i > sourceWdists_.size() - 1)
            {
                i = sourceWdists_.size() - 1;
            }
            return label(sourceWdists_[i]/deltaY_);
        }

        bool master()
        {
            if (myid_ == masterid_)
            {
                return true;
            }
            return false;
        }

        label wait()
        {
            if (!parRun())
            {
                return 0;
            }
            else
            {
                return gMax(procs_);
            }
        }

        static word int2str(label i)
        {
            std::stringstream ss;
            ss << i;
            return ss.str();
        }

        static word dbl2str(scalar d)
        {
            std::stringstream ss;
            ss << d;
            return ss.str();
        }

        word caseName()
        {
            if (!parRun())
            {
                std::vector<word> vtmp;
                stringtok(vtmp, casedir_, "/");

                word cname = "fields";
                if (vtmp.size() == 1)
                {
                    cname += casedir_;
                }
                else
                {
                    for (unsigned i = 0; i < vtmp.size(); i++)
                    {
                        cname += ("_" + vtmp[i]);
                    }
                }
                return cname;
            }
            word ss = "processor";
            ss += int2str(myid_);
            return ss;
        }

        void parInit(label np);

        void getDomainFields(); //for parallel run, get domain fields from source case
        void getDomainFields0();//single-field parallel model

        void createSource
        (
            const fvMesh* mesh,
            const Time* runTime
        );

        void createTarget
        (
            const fvMesh* mesh,
            const Time* runTime
        );

        void clearSourceTarget();

        static void stringtok
        (
            std::vector<word>& container,
            const word& in,
            const word& delimiters
        );

        label search_error();
        scalar map_error(const word& scname); //mapping error based on a scalar
        scalar Umap_error(const word& vecname);

        void mapScalarField(const word& scname);
        void mapVectorField(const word& vecname);
        void mapScalarFieldPar(const word& scname);
        void mapVectorFieldPar(const word& vecname);

        void getFieldTypes
        (
            const fvMesh& mesh,
            const word& timeName
        );

        void scalePressure
        (
            const dimensionSet& sourceDim,
            const dimensionSet& targetDim
        );

        void checkFieldDimensions();

        void mapFields();

        void setMapTime
        (
            Time& runTime
        );

        virtual void setInput
        (
            const dictionary& dict
        );

        virtual void setInput
        (
            const dictionary& dict,
            const argList& args
        );


        bool parRun()
        {
            if (nprocs_ > 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        class readFields& source()
        {
            return source_();
        }

        class readFields& target()
        {
            return target_();
        }

        virtual void generateNodes(const fvMesh* mesh);

        virtual void setOptions
        (
            const argList& args
        );

        label cellIndex
        (
            scalar alphai,
            label kdTreeId
        )
        {
            if (alphai > sourceAlphaMax_ || !wdistMap_)
            {
                return cellMap_[kdTreeId];
            }
            else
            {
                label ialpha = int(alphai/deltaY_);
                return alphaMap_(ialpha, kdTreeId);
            }

        }

    inline void search_nearest
    (
        std::vector<label>& nb,
        std::vector<scalar>& dist,
        scalar x,
        scalar y,
        scalar z
    );

    bool indomain
    (
        const point& pt
    )
    {
        return target().indomain(pt);
    }

    inline void search_nearest
    (
        std::vector<label>& nb,
        std::vector<scalar>& dist,
        scalar alphai,
        scalar x,
        scalar y,
        scalar z
    );

    foamMap();

    virtual ~foamMap(){;}
};


void foamMap::search_nearest
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar wdisti,
    scalar x,
    scalar y,
    scalar z
)
{
    if (!wdistMap_)
    {
        return search_nearest
        (
            nb,
            dist,
            x,
            y,
            z
        );
    }

    unsigned int ic = label(wdisti/deltaY_);
    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;
    if (ic > kdTrees_.size() - 1) ic = kdTrees_.size() - 1;
    if (!kdTrees_[ic].active)
    {
        return internalMap_.search
        (
            qvect,
            nb,
            dist
        );
    }

    if (wdisti > sourceAlphaMax_ || ic > kdTrees_.size() - 1)
    {
        return internalMap_.search
        (
            qvect,
            nb,
            dist
        );
    }
    else
    {
        if (kdTrees_[ic].active)
        {
            kdTrees_[ic].search
            (
                qvect,
                nb,
                dist
            );
        }
        else
        {
            internalMap_.search
            (
                qvect,
                nb,
                dist
            );
        }
    }

    return;
}

void foamMap::search_nearest
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar x,
    scalar y,
    scalar z
)
{
    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;
    internalMap_.search
    (
        qvect,
        nb,
        dist
    );
}

}//Foam

