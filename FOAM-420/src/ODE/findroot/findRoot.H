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
    (c) held by original author

Application
    findRoot

Description
    root-finding functions

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "findroot/BisectionRoot.H"
#include "findroot/RiddersRoot.H"
#include "findroot/NewtonRoot.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class rootFunction
{
public:

    rootFunction();

    rootFunction(double a1, double a2, double a3, double a4 );

    void modify(double x,double y,double z, double w);

    double operator()(const double& x) const; //returns the value of the function

    double d(const double& x) const; // returns the value of the differential of the function

private:
//  This function has the formulation for the transition model.
    double coef1;
    double coef2;
    double coef3;
    double coef4;
};







// ************************************************************************* //
