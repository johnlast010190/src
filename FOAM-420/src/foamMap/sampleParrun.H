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
    samplePar.H

\*---------------------------------------------------------------------------*/

#pragma once
#include "blockSample.H"

//#include <glob.h>


namespace Foam
{

//transfer a field to a specific processor
template<class T>
void interProcTrans
(
    DynamicList<T> &field,
    const DynamicList<T> &source,
    const word &name
);

class samplePar : public foamSample
{
    public:
        word sampleMethod_;
        DynamicList<point> bodySurfCenter_;

        void setOptions
        (
            const argList& args
        );

        void getDirNames
        (
            wordvec &dnames,
            const word &pattern
        );

        void getMapTime();
        void combineFields();
        void combineFieldsBK();
        void constructKnn();
        void getBodySurface();

        void searchNodeField
        (
            gridField &grid,
            const word &funct="field"
        );

        samplePar();
        samplePar
        (
            const Time* runTime,
            const fvMesh *_mesh
        );
        virtual ~samplePar() {;}
};


}//Foam
