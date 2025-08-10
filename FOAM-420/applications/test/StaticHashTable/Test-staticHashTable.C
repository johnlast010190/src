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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "containers/HashTables/StaticHashTable/StaticHashTable.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOstreams/StringStreams/IStringStream.H"
#include "db/IOstreams/StringStreams/OStringStream.H"

using namespace Foam;

// use define so we can easily test other implementations
#define HASHTABLE_CLASS StaticHashTable

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    HASHTABLE_CLASS<double> table1(13);

    table1.insert("aaa", 1.0);
    table1.insert("aba", 2.0);
    table1.insert("aca", 3.0);
    table1.insert("ada", 4.0);
    table1.insert("aeq", 5.0);
    table1.insert("aaw", 6.0);
    table1.insert("abs", 7.0);
    table1.insert("acr", 8.0);
    table1.insert("adx", 9.0);
    table1.insert("aec", 10.0);

    table1.erase("aaw");
    table1.erase("abs");

    Info<< "\ntable1 toc: " << table1.toc() << endl;
    table1.printInfo(Info)
        << "table1 [" << table1.size() << "] " << endl;
    forAllIter(HASHTABLE_CLASS<double>, table1, iter)
    {
        Info<< iter.key() << " => " << iter() << nl;
    }

    table1.set("acr", 108);
    table1.set("adx", 109);
    table1.set("aec", 100);
    table1("aaw") -= 1000;
    table1("aeq") += 1000;

    Info<< "\noverwrote some values table1: " << table1 << endl;

    Info<< "\ntest find:" << endl;
    Info<< table1.find("aaa")() << nl
        << table1.find("aba")() << nl
        << table1.find("aca")() << nl
        << table1.find("ada")() << nl
        << table1.find("aeq")() << nl
        << table1.find("acr")() << nl
        << table1.find("adx")() << nl
        << table1.find("aec")() << nl
        << table1["aaa"] << nl;

    {
        OStringStream os;
        os  << table1;
        HASHTABLE_CLASS<double> readTable(IStringStream(os.str())(), 100);

        Info<< "Istream constructor:" << readTable << endl;
    }


    HASHTABLE_CLASS<double> table2(table1);
    HASHTABLE_CLASS<double> table3(table1.xfer());

    Info<< "\ncopy table1 -> table2" << nl
        << "transfer table1 -> table3 via the xfer() method" << nl;

    Info<< "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl
        << "\ntable3" << table3 << nl;

    Info<< "\nerase table2 by iterator" << nl;
    forAllIter(HASHTABLE_CLASS<double>, table2, iter)
    {
        Info<< "erasing " << iter.key() << " => " << iter() << " ... ";
        table2.erase(iter);
        Info<< "erased" << endl;
    }

    Info<< "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl
        << "\ntable3" << table3 << nl;

    table3.resize(1);
    Info<< "\nresize(1) table3" << nl;
    table3.printInfo(Info)
        << table3 << nl;

    table3.resize(10000);
    Info<< "\nresize(10000) table3" << nl;
    table3.printInfo(Info)
        << table3 << nl;

    HASHTABLE_CLASS<double> table4;

    table4 = table3;
    Info<< "\ncopy table3 -> table4 " << table4 << nl;

    Info<< "\nclear table4 ... ";
    table4.clear();
    Info<< "[" << table4.size() << "] " << table4 << nl;

    table1 = table3;
    Info<< "\ncopy table3 -> table1 (previously transferred)" << table1 << nl;

    Info<< "test table1 == table3 : " << (table1 == table3) << nl;
    table1.erase(table1.begin());
    Info<< "removed an element - test table1 != table3 : "
        << (table1 != table3) << nl;

    // insert a few things into table2
    table2.set("ada", 14.0);
    table2.set("aeq", 15.0);
    table2.set("aaw", 16.0);
    table2.set("abs", 17.0);
    table2.set("adx", 20.0);

    Info<< "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl;

    label nErased = table1.erase(table2);

    Info<< "\nerase table2 keys from table1 (removed "
        << nErased << " elements)" << nl
        << "\ntable1" << table1 << nl
        << "\ntable2" << table2 << nl;


    Info<< "\ntable3" << table3
        << "\nclearStorage table3 ... ";
    table3.clearStorage();
    Info<< table3 << nl;

    Info<< "\nDone\n";

    return 0;
}


// ************************************************************************* //
