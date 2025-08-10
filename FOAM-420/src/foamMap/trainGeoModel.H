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
    (c) 2009-2021 Esi Ltd.

Description
    Train a machine-learning model to recognize the closest flow topology
    in the database to the one in consideration.

SourceFiles
    trainGeoModel.C
\*---------------------------------------------------------------------------*/

#pragma once

#include "dataItem.H"
#include "triSurface/triSurface.H"

namespace Foam
{
class trainGeoModel
:
    public foamSample
{
	private:
		word task_;
		label K_;

	public:
		std::map<word,scalar> caseLabels_; //enumeriate case name with number for training
		std::map<word,dataItem> caseItems_;
		std::map<label,label> idmap_; //map between search-returned id and the id in vcaseItems_

		std::vector<dataItem> vcaseItems_;//vector version of the data items
		std::vector<word> caseNames_;
		std::vector<scalar> geoDists_;

		void task(const word& tsk)
		{
			task_=tsk;
		}

		void clearGeoData(); //clear mask array once database search tree is build

		autoPtr<KNN> model_;

		trainGeoModel();
        trainGeoModel
        (
            const Time* runTime,
            const dictionary& dict
        );

		trainGeoModel
		(
			const Time* runTime,
			const fvMesh* mesh
		);

		void setOptions
		(
			const DynamicList<word>& argNames,
			const DynamicList<word>& argValues
		);

		void readDatabase1(const word& dbsfile);
		void readDatabase(const word& dbsfile);
		void readCaseList();
		void getCaseLabels();
		void readGeometryFile();
		void readVector
		(
			std::vector<scalar>& vect,
            std::ifstream& is
		);

		void build_model();
		void build_model
        (
            label m,
            label n,
            label k
        );//build search model for specific sample resolution

		void save_model();
		void getSolidLabels
		(
			gridField& grid,
			std::vector<label>& in_solids
		);

		void load_model();

		void get_prediction();
		scalar getDist(label i, label j);

		void getGeomSample(); //get sampling data for geometry topology only.
        virtual ~trainGeoModel(){;}

};


}// End namespace Foam


// ************************************************************************* //



