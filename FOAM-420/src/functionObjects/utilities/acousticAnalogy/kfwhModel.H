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
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::kfwhModel

Group
    grpUtilityFunctionObjects

Description
	A founction object which implements an extended aero-acoustics model
	as defined in : Journal of Sound and Vibration (1997), 202 (4), 491-509
	A New Boundary Integral Formulation for the Prediction of Sound Radiation
	By P. di Francescantonio

SourceFiles
    kfwhModel.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "acousticAnalogy/FfowcsWilliamsHawkings.H"
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "triSurface/triSurface.H"
#include "array2d.H"
#include "array3d.H"
#include "src/knn.h"
#include "acousticAnalogy/solidObject.H"

namespace Foam
{
namespace functionObjects
{
//extended fwh model

class kfwhModel: public FfowcsWilliamsHawkings
{

protected:
	point prefLoc_;
	label icPref_; //reference pressure cell
	bool viscousEffect_;
	label hpcMode_;
	label phaseNo_;

	scalar veloCoef_;
	scalar presCoef_;
	scalar tkeCoef_;
	scalar term3Coef_;

	bool meshChange_;
	word prefOpt_;
	// word modelType_;
	List<List<scalarList>> pLoad_;
	List<List<scalarList>> pNonlinear_;

	Array2d<scalar> L_;//Lighthill stress
	Array2d<scalar> L0_;
	Array2d<scalar> L00_;
	Array2d<scalar> dotL_;
	Array2d<scalar> Tij_;

	vectorField Uf_; //rho/rho0*U: at surface faces
	vectorField Uf0_;
	vectorField Uf00_;

	labelList outCells_;

	Array3d<scalar> ovterm1_;
	Array3d<scalar> ovterm2_;
	Array3d<scalar> dvfield10_;

	scalar pref_;
	word boundSurfaceFile_;
	autoPtr<triSurface> surfaces_;
	autoPtr<boundBox> samplebox_;
	KNN surfaceKnn_;

	vectorField faceNormals_;
	vectorField searchZoneXyz_;
	labelList cellIds_;
	labelList procIds_;
	//surface face id to cell id map
	labelList faceToCell_;//index is face id, cell id is the index of cellIds_
	labelList insolid_; //whether a sample surface face is in a solid
	//labelHashSet solidPatchs_;
	solidObject sobject_;


public:
	 //- Runtime type information
    TypeName("kfwhModel");
	const label tensorIdx[9]={0,1,2,1,3,4,2,4,5};
	kfwhModel
	(
		const word& name,
		const Time& runTime,
		const dictionary& dict
	);

	void setTensorField();

	scalarField rho();
	void calculate();
	const tmp<volSymmTensorField> R() ;

	tmp<volScalarField>  turbEng();

	virtual ~kfwhModel(){;}
	void initialise();

	inline scalar delta(label i, label j);
	inline label ij(label i,label j);

	inline scalar Lij(label k, label i, label j);
	inline scalar Tij(label k,label i,label j);
	inline scalar Tii(label k);
	inline scalar Trr(label k,const point &r);
	inline point rvect(label obsSetI, label obsI, const point &fldPt);
	inline const vectorField& Sf() const;
	inline const vectorField& Cf() const;
	inline const scalarField& magSf() const;
	inline bool inSurfObject(const point &pt); //test whether a point is inside object enclosed by surface
	inline bool	inSearchZone(const point &pt); //test if a field point is inside the search zone
	inline const point &nf(label jf);
	inline scalar Fi(label jf,label i);
	inline scalar dotFi(label jf,label i);

	inline scalar Firi(label jf,const point &rv);
	inline scalar dFiri(label jf,const point &rv);

	void getSampleBox();
	void writeFfowcsWilliamsHawkings();

	void setSurfaceState(); //exclude the surface which cut into solid objects

	bool read(const dictionary&);
	word surfaceFile();
	void buildSurfaceKnn();
	void buildCellMapKnn(KNN &knn);
	void getFaceNormals();
	void getCellList();
	void getFaceCellMap(KNN &knn);
	void getVolTerms(bool hasOldField);
};


scalar kfwhModel::dFiri(label jf,const point &rv)
{
	scalar fsum = 0;
	for (label i=0;i<3;i++)
	{
		scalar r = mag(rv);
		scalar rhati = rv[i]/r;
		fsum += dotFi(jf,i)*rhati;
	}
	return fsum;
}


scalar kfwhModel::Firi(label jf,const point &rv)
{
	scalar fsum = 0;
	for (label i=0;i<3;i++)
	{
		scalar r = mag(rv);
		scalar rhati = rv[i]/r;
		fsum += Fi(jf,i)*rhati;
	}
	return fsum;
}


scalar kfwhModel::dotFi(label jf,label i)
{
	scalar fdotSum = 0;
	for (label j=0;j<3;j++)
	{
		fdotSum += dotL_(jf,ij(i,j))*faceNormals_[jf][j];
	}
	return fdotSum;
}


scalar kfwhModel::Fi(label jf,label i)
{
	scalar fsum = 0;
	for (label j=0;j<3;j++)
	{
		fsum += L_(jf,ij(i,j))*faceNormals_[jf][j];
	}
	return fsum;
}


const point & kfwhModel::nf(label jf)
{
	return faceNormals_[jf];
}


scalar kfwhModel::Lij
(
	label k,
	label i,
	label j
)
{
	return L_(k,ij(i,j));
}


bool kfwhModel::inSearchZone
(
	const point &pt
)
{
	bool contains = samplebox_().contains(pt);
	bool inSurfobj = inSurfObject(pt);
	if (contains && !inSurfobj)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool kfwhModel::inSurfObject
(
	const point &pt
)
{
	std::vector<label> nb;
	std::vector<scalar> dist;
	std::vector<scalar> qvect(3);
	qvect[0] = pt.x();
	qvect[1] = pt.y();
	qvect[2] = pt.z();
	surfaceKnn_.search
	(
		qvect,
		nb,
		dist
	);
	label j0 = nb[0];
	const point& nsf = faceNormals_[j0];
	point dpt = Cf()[j0] - pt;
	scalar d = nsf&dpt;
	if (d>=0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


const vectorField& kfwhModel::Sf() const
{
	return surfaces_().Sf();
}


const vectorField& kfwhModel::Cf() const
{
	return surfaces_().Cf();
}


const scalarField& kfwhModel::magSf() const
{
	return surfaces_().magSf();
}


scalar kfwhModel::Tij(label k,label i,label j)
{
	return Tij_(k,ij(i,j));
}


point kfwhModel::rvect
(
	label obsSetI,
	label obsI,
	const point &fldPt
)
{
	const point& pt = observerSets_[obsSetI].positions()[obsI];
	point dpt = pt - fldPt;
	return dpt;
}


scalar kfwhModel::Trr(label k,const point &r)
{
	scalar tsum = 0;
	scalar rr = mag(r);
	for (label i=0;i<3;i++)
	{
		for (label j=0;j<3;j++)
		{
			tsum += Tij(k,i,j)*r[i]*r[j]/(rr*rr);
		}
	}
	return tsum;
}


scalar kfwhModel::delta(label i, label j)
{
	if (i==j)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


label kfwhModel::ij(label i,label j)
{
	label i0 = 3*i+j;
	return tensorIdx[i0];
}


scalar kfwhModel::Tii(label k)
{
	return Tij(k,0,0) + Tij(k,1,1) + Tij(k,2,2);
}

}//functionObjects
}//Foam
