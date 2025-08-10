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

Description:
    Base class for mapping two flow fields with different meshes

SourceFile
    readFieldsPar.C
\*---------------------------------------------------------------------------*/

#pragma once

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "src/knn.h"
#include "array2d.H"
#include "smBoundBox.H"

namespace Foam
{

typedef std::vector<std::vector<label>> iarray2d;
typedef std::vector<word> words;

typedef std::vector<std::vector<scalar>> farray2d;
typedef unsigned short ushort;

typedef std::map<word, DynamicList<scalar>> ScalarFields;
typedef std::map<word, DynamicList<point>> VectorFields;


class foamMap;

class readFields
{
    protected:
        // autoPtr<boundBox> bbox_;

        autoPtr<smBoundBox> bbox_;
        autoPtr<boundBox> gbox_;
        DynamicList<point> xyz_;

        word type_; //source,target
    public:
        ScalarFields scalarfields_;
        VectorFields vectorfields_;

        word mapTime_;
        const fvMesh* mesh;
        const Time* runTime_;
        foamMap* fieldmap_;
        bool norm_;
        std::vector<KNN> kdTrees_;

        KNN boundaryMap_;

        KNN internalMap_;
        scalar alphaMax_;
        std::map<label,label> cellMap_;//map from kdTree id to cell id
        label ncells_;

        Array2d<label> alphaMap_;

        DynamicList<scalar> wdists_;
        scalar deltaY_;
        scalar xmin_;
        scalar xmax_;
        scalar ymin_;
        scalar ymax_;
        scalar zmin_;
        scalar zmax_;
        HashSet<word> scfields_;
        HashSet<word> vecfields_;
        scalar rhoRef_;
        scalar UrotDegreeFromSource_;
        vector Uref_;
        bool interp_;
        scalar maxWallDist_;

        DynamicList<point>& xyz()
        {
            return xyz_;
        }

        const DynamicList<point>& xyz() const
        {
            return xyz_;
        }

        word myType()
        {
            return type_;
        }

        readFields
        (
            const fvMesh* _mesh,
            const Time* _runTime,
            foamMap *_map,
            const word& _type = "source"
        );

        readFields();

        void setInputs();

        void setNorm(bool norm)
        {
            norm_ = norm;
        }

        void constructBoundBox
        (
            scalar,
            const word,
            label
        );

        void transformPoints
        (
            const boundBox& tbox,
            DynamicList<point>& xyz
        );

        const smBoundBox& bbox() const
        {
            return bbox_();
        }

        const boundBox& gbox() const
        {
            return gbox_();
        }

        point pbar(const point& pt)
        {
            point pp(xbar(pt.x()), ybar(pt.y()), zbar(pt.z()));
            return pp;
        }

        scalar U0()
        {
            return mag(Uref_);
        }
        scalar pref()
        {
            return U0()*U0()*rhoRef_;
        }
        scalar scale(const word& fldname);

        scalar xbar(scalar x)
        {
            if (norm_)
            {
                return (x - xmin_)/(xmax_ - xmin_);
            }
            else
            {
                return x;
            }
        }

        scalar ybar(scalar y)
        {
            if (norm_)
            {
                return (y - ymin_)/(ymax_ - ymin_);
            }
            else
            {
                return y;
            }
        }

        scalar zbar(scalar z)
        {
            if (norm_)
            {
                return (z - zmin_)/(zmax_ - zmin_);
            }
            else
            {
                return z;
            }
        }

        label bfaceIndex(label j)
        {
            return ncells_ + j;
        }

        inline label cellIndex
        (
            scalar alphai,
            label kdTreeId
        );

        void construct_boundary_knn();

        void construct_internal_knn
        (
            scalar cutoff = 0.8
        );

        label numBoundFaces();

        void getBoundBox
        (
            const boundBox& box
        );

        void createFields(const word& timeName = "");

        void buildKdTrees
        (
            const volScalarField& y,
            scalar cutoff = 0.8,
            const word& excludeBndType = "none"
        );

        label inGroup(label i)
        {
            return label(wdists_[i]/deltaY_);
        }

        void setWallDist
        (
            const volScalarField& y,
            scalar maxwdist
        );

        volScalarField& getScalarField(const word& name)
        {
            return mesh->lookupObjectRef<volScalarField>(name);
        }

        surfaceScalarField& getFlux(const word& name)
        {
            return mesh->lookupObjectRef<surfaceScalarField>(name);
        }

        const surfaceScalarField& getFlux(const word& name) const
        {
            return mesh->lookupObjectRef<surfaceScalarField>(name);
        }

        scalar getScalarField(const word& name, label i)
        {
            return mesh->lookupObjectRef<volScalarField>(name).primitiveField()[i];
        }

        const volScalarField& getScalarField(const word& name) const
        {
            return mesh->lookupObjectRef<volScalarField>(name);
        }

        volVectorField& getVectorField(const word& name)
        {
            return mesh->lookupObjectRef<volVectorField>(name);
        }

        const volVectorField& getVectorField
        (
            const word& name
        )
        const
        {
            return mesh->lookupObjectRef<volVectorField>(name);
        }

        const point& getVectorField
        (
            const word& name,
            label i
        )
        const
        {
            return mesh->lookupObjectRef<volVectorField>(name).primitiveField()[i];
        }

        void storeFields(); //store source fields
        void storeFields
        (
            const boundBox& sampleBox
        );

        inline void search_nearest
        (
            //search is always based on relative coordinates
            std::vector<label>& nb,
            std::vector<scalar>& dist,
            scalar alphai,
            scalar x,
            scalar y,
            scalar z
        );

        inline void search_nearest_bface
        (
            std::vector<label>& nb,
            std::vector<scalar>& dist,
            scalar x,
            scalar y,
            scalar z
        );

        inline bool indomain(const point& pt);
        virtual ~readFields();

};

bool readFields::indomain(const point& pt)
{
    return bbox().contains(pt);
}


void readFields::search_nearest_bface
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar x,
    scalar y,
    scalar z
)
{
    if (type_ != "source") return;
    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;
    boundaryMap_.search
    (
        qvect,
        nb,
        dist
    );
    return;
}


void readFields::search_nearest
(
    std::vector<label>& nb,
    std::vector<scalar>& dist,
    scalar wdisti,
    scalar x,
    scalar y,
    scalar z
)
{
    if (type_ != "source") return;

    std::vector<scalar> qvect(3);
    qvect[0] = x;
    qvect[1] = y;
    qvect[2] = z;

    if (wdisti > alphaMax_)
    {
        internalMap_.search
        (
            qvect,
            nb,
            dist
        );
    }
    else
    {
        label ic = label(wdisti/(deltaY_ + VSMALL));
        kdTrees_[ic].search
        (
            qvect,
            nb,
            dist
        );
    }

    return;
}

label readFields::cellIndex
(
    scalar alphai,
    label kdTreeId
)
{
    if (alphai > alphaMax_)
    {
        return cellMap_[kdTreeId];
    }
    else
    {
        label ialpha = int(alphai/deltaY_);

        if (kdTreeId < label(alphaMap_[ialpha].size()))
        {
            return alphaMap_(ialpha, kdTreeId);
        }
        else
        {
            return -1;
        }
    }
}


}//Foam
