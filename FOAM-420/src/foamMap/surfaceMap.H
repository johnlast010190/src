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

#pragma once
#include <unistd.h>
#include "foamMap.H"
#include "blockSample.H"
#include "meshes/primitiveShapes/plane/plane.H"

namespace Foam
{
class surfaceMap: public foamMap
{
    protected:
        autoPtr<plane> mirror_;
        const fvMesh *mesh_;
        const Time *runTime_;
        bool imap_;
        bool bmap_;
        bool byType_;
        bool byPhysicalType_;
        bool byNameList_;
        word bodySurfaceType_;
    public:
        void defaults();

        std::vector<word> sourcePatchNames_;//patch name separated by ,

        std::map<word,KNN> bndMaps_;
        std::map<word,std::vector<label>> faceIdMaps_;//bnd tree id => boundary region face id maps
        std::map<word,DynamicList<point>> bsource_xyz_;//boundary face centres
        std::map<word,DynamicList<scalar>> bsourceScalarFields_;//boundary scalar field, bndName,scalarName: field

        bool internalSym()
        {
            return imap_;
        }

        bool bndSym()
        {
            return bmap_;
        }

        void setInput
        (
            const dictionary& dict
        );

        void setInput
        (
            const dictionary& dict,
            const argList& args
        );

        void createMirrorFields
        (
            const word &mapTime
        );

        plane &mirror()
        {
            return mirror_();
        }

        label patchId(const word &name);

        word mirrorName(const word &name)
        {
            return name+"_symp";
        }

        void getParallelFields();
        void getMirrorFields();
        void getSymFields();
        void buildSearchTree();
        void buildParSearchTree();

        surfaceMap
        (
            const fvMesh* mesh,
            const Time* runTime
        );

        surfaceMap
        (
            const fvMesh& mesh
        );

        surfaceMap();
        virtual ~surfaceMap(){;}
};



}//Foam
