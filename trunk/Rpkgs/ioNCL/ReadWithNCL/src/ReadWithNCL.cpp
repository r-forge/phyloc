// RcppExample.cpp: Part of the R/C++ interface class library, Version 5.0
//
// Copyright (C) 2005-2006 Dominick Samperi
//
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Lesser General Public License as published by 
// the Free Software Foundation; either version 2.1 of the License, or (at 
// your option) any later version.
//
// This library is distributed in the hope that it will be useful, but 
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License 
// along with this library; if not, write to the Free Software Foundation, 
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 

#include <fstream>

#include "Rcpp.hpp"

//TODO - figure out where/why length is getting defined as Rf_length so that
//this isn't necessary to compile
#define  Rf_length length

#include "ncl.h"
#include "NCLInterface.h"

using namespace std;

//This function receives a list of parameters from
//R, which can then be extracted below
RcppExport SEXP ReadWithNCL(SEXP params) {

    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {


		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
		string filename = rparam.getStringValue("filename");
	
		BASICCMDLINE reader;
		
		//this is where the reader would be passed the filename to read
		//reader.Run(filename.c_str()); //Will not compile
		//reader.Run(NULL);
		reader.Initialize(const_cast < char* > (filename.c_str()));
	
		string filenameString = "I was told to read a file named ";
		filenameString += filename;
	
		//Various calls to the reader can be made here to fill the various strings
		//string treeString = "This was passed in the tree string";
		NxsString treeStringNxs;
		reader.RReturnTrees(treeStringNxs);
		string treeString=treeStringNxs.c_str();
	
		string otherString = "";
	
		//string discreteString = "This was passed in the discrete string";
		
		NxsString characterStringNxs;
		reader.RReturnCharacters(characterStringNxs,false, true, false);
		//string characterString=characterStringNxs.c_str();
		
		string continuousString = "This was passed in the continuous string";
		//reader.GetContinuousString();
	
		//Just to test to see if reader is running
		string testString = reader.TestRunning();
	
		// Build result set to be returned as a list to R.
		RcppResultSet rs;
		
		//if the various strings are nonempty, add them to the ResultSet
		if(filenameString.length() > 0)
			rs.add("filenamestring", filenameString);
		if(treeString.length() > 0)
			rs.add("treestring", treeString);
		if(otherString.length() > 0)
			rs.add("otherstring", otherString);
		if(characterStringNxs.length() > 0)
			rs.add("characterstring", characterStringNxs.c_str());
		if(continuousString.length() > 0)
			rs.add("continuousstring", continuousString);
		if(testString.length() > 0)
			rs.add("teststring", testString);

		
		// Get the list to be returned to R.
		rl = rs.getReturnList();
	
    	} catch(std::exception& ex) {
			exceptionMesg = copyMessageToR(ex.what());
    	} catch(...) {
			exceptionMesg = copyMessageToR("unknown reason");
    	}
    
    if(exceptionMesg != NULL)
	error(exceptionMesg);

    return rl;

}

RcppExport SEXP ReadTreesWithNCL(SEXP params) {
	
    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {
		
		
		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
		string filename = rparam.getStringValue("filename");
		
		BASICCMDLINE reader;
		
		//this is where the reader would be passed the filename to read
		//reader.Run(filename.c_str()); //Will not compile
		//reader.Run(NULL);
		reader.Initialize(const_cast < char* > (filename.c_str()));
		
		NxsString treeStringNxs;
		reader.RReturnTrees(treeStringNxs);
		string treeString=treeStringNxs.c_str();
		
		// Build result set to be returned as a list to R.
		RcppResultSet rs;
		
		//if the various strings are nonempty, add them to the ResultSet
		if(treeString.length() > 0)
			rs.add("treestring", treeString);
			
		
		// Get the list to be returned to R.
		rl = rs.getReturnList();
		
	} catch(std::exception& ex) {
		exceptionMesg = copyMessageToR(ex.what());
	} catch(...) {
		exceptionMesg = copyMessageToR("unknown reason");
	}
    
    if(exceptionMesg != NULL)
		error(exceptionMesg);
	
    return rl;
	
}

RcppExport SEXP ReadCharsWithNCL(SEXP params) {
	
    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {
		
		
		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
		string filename = rparam.getStringValue("filename");
		bool allchar = rparam.getBoolValue("allchar");
		bool levelsall=rparam.getBoolValue("levelsall");
		bool polymorphictomissing=rparam.getBoolValue("polymorphictomissing");
		
		BASICCMDLINE reader;
		
		//this is where the reader would be passed the filename to read
		//reader.Run(filename.c_str()); //Will not compile
		//reader.Run(NULL);
		reader.Initialize(const_cast < char* > (filename.c_str()));
		
		NxsString charStringNxs;
		reader.RReturnCharacters(charStringNxs,allchar, polymorphictomissing, levelsall);
		string charString=charStringNxs.c_str();
		
		// Build result set to be returned as a list to R.
		RcppResultSet rs;
		
		//if the various strings are nonempty, add them to the ResultSet
		if(charString.length() > 0)
			rs.add("charstring", charString);
		
		
		// Get the list to be returned to R.
		rl = rs.getReturnList();
		
	} catch(std::exception& ex) {
		exceptionMesg = copyMessageToR(ex.what());
	} catch(...) {
		exceptionMesg = copyMessageToR("unknown reason");
	}
    
    if(exceptionMesg != NULL)
		error(exceptionMesg);
	
    return rl;
	
}


/*
 * Sample function illustrates how to use the Rcpp R/C++ interface library.
 */
RcppExport SEXP Rcpp_Example(SEXP params, SEXP nlist, 
			     SEXP numvec, SEXP nummat,
			     SEXP df, SEXP datevec, SEXP stringvec,
			     SEXP fnvec, SEXP fnlist) {

    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
/*
    try {

	int i=0, j=0;

	// Get parameters in params.
	RcppParams rparam(params);
	string method = rparam.getStringValue("method");
	double tolerance = rparam.getDoubleValue("tolerance");
	int    maxIter = rparam.getIntValue("maxIter");
	RcppDate startDate = rparam.getDateValue("startDate");

	// The output of commands like this may not appear under Windows.
	Rprintf("Start Date: %d/%d/%d\n", 
		startDate.getMonth(),
		startDate.getDay(),
		startDate.getYear());

	// QuantLib note: an RcppDate is automatically converted to QuantLib
	// Date when the context calls for this, provided 
	// USING_QUANTLIB is set.

	RcppDateVector dateVec(datevec);
	//dateVec(0) = RcppDate(12, 15, 1989); // update one element.

	RcppStringVector stringVec(stringvec);
	//stringVec(1) = string("New York"); // update one element.

	// and nl.getValue(i) to fetch data.
	RcppNumList nl(nlist);

	// numvec parameter viewed as vector of ints (with truncation).
	//RcppVector<int> vecI(numvec);

	// mat parameter viewed as matrix of ints (with truncation).
	//RcppMatrix<int> matI(nummat);

	// vec parameter viewed as vector of doubles.
	RcppVector<double> vecD(numvec);

	// mat parameter viewed as matrix of doubles.
	RcppMatrix<double> matD(nummat);

	// Do some computations with the matrices.
	int nrows = matD.getDim1();
	int ncols = matD.getDim2();
	for(i = 0; i < nrows; i++)
	    for(j = 0; j < ncols; j++)
		matD(i,j) = 2 * matD(i,j);

	int len = vecD.size();
	for(i = 0; i < len; i++)
	    vecD(i) = 3 * vecD(i);

	// Get copy of matrix/vector in standard (unchecked) C/C++ format.
	// May be useful when these vectors need to be passed to
	// C/C++ code that does not know about Rcpp classes...
	double **a = matD.cMatrix();
	double  *v = vecD.cVector();

	// ...or we might want to use an STL container...
	vector<double> stlvec(vecD.stlVector());
	nrows = (int)stlvec.size();
	for(i = 0; i < nrows; i++)
	    stlvec[i] += 1;

	// ...or perhaps a container of containers.
	vector<vector<double> > stlmat(matD.stlMatrix());
	nrows = (int)stlmat.size();
	ncols = (int)stlmat[0].size();
	for(i = 0; i < nrows; i++)
	    for(j = 0; j < ncols; j++)
		stlmat[i][j] += 2;

	// Get a zero matrix the same size as matD.
	//RcppMatrix<double> matZ(nrows, ncols);

	// Make a vector of strings
 	vector<string> svec(2);
        svec[0] = "hello";
	svec[1] = "world";

	// Process the input data frame and show factors and dates.
	RcppFrame inframe(df);
*/
	/*
	Rprintf("\nFactors and Dates in frame...");
	vector<vector<ColDatum> > table = inframe.getTableData();
	int nrow = table.size();
	int ncol = table[0].size();
	for(int row=0; row < nrow; row++) {
	    for(int col=0; col < ncol; col++) {
		RcppDate d;
		string name;
		int level;
		switch(table[row][col].getType()) {
		case COLTYPE_FACTOR:
		    level = table[row][col].getFactorLevel();
		    name  = table[row][col].getFactorLevelName();
		    Rprintf("Level, name: %d, %s\n",
		            level, name.c_str());
		    break;
		case COLTYPE_DATE:
		    d = table[row][col].getDateValue();
		    Rprintf("Start Date: %d/%d/%d\n", 
		            d.getMonth(),
			    d.getDay(),
			    d.getYear());
		    break;
		default:
		    ; // Ignore other types.
		}
	    }
	}
	*/
/*
	// Make a pre-data frame, that is, a list object that when passed
	// the the R function data.frame() will return a data frame with
	// the specified column names and data types. The first row added
	// determines the types for all columns.
	int numCol=4;
	vector<string> colNames(numCol);
	colNames[0] = "alpha"; // column of strings
	colNames[1] = "beta";  // column of reals
	colNames[2] = "gamma"; // factor column
	colNames[3] = "delta"; // column of Dates
	RcppFrame frame(colNames);

	// Third column will be a factor. In the current implementation the
	// level names are copied to every factor value (and factors
	// in the same column must have the same level names). The level names
	// for a particular column will be factored out (pardon the pun) in
	// a future release.
	int numLevels = 2;
	string *levelNames = new string[2];
	levelNames[0] = string("pass"); // level 1
	levelNames[1] = string("fail"); // level 2

	// First row (this one determines column types).
	vector<ColDatum> row1(numCol);
	row1[0].setStringValue("a");
	row1[1].setDoubleValue(3.14);
	row1[2].setFactorValue(levelNames, numLevels, 1);
	row1[3].setDateValue(RcppDate(7,4,2006));
	frame.addRow(row1);

	// Second row.
	vector<ColDatum> row2(numCol);
	row2[0].setStringValue("b");
	row2[1].setDoubleValue(6.28);
	row2[2].setFactorValue(levelNames, numLevels, 1);
	row2[3].setDateValue(RcppDate(12,25,2006));
	frame.addRow(row2);

	// Done with levelNames.
	delete [] levelNames;

	// Test MyRVectorFunction defined above...
	MyRVectorFunc vfunc(fnvec);
	int n = 10;
	vector<double> vecInput(n);
	for(int i=0; i < n; i++)
	    vecInput[i] = i;
	double vecSum = vfunc.getSum(vecInput);
	Rprintf("vecSum = %lf\n", vecSum);
	
	// Test MyRListFunction defined above...
	MyRListFunc lfunc(fnlist);
	double alpha=1, beta=2, gamma=3;
	vector<double> vecOut = lfunc.addOne(alpha, beta, gamma);
	Rprintf("vecOut: %lf, %lf, %lf\n", vecOut[0], vecOut[1], vecOut[2]);

	RcppDate aDate(12, 25, 1999);

	// Build result set to be returned as a list to R.
	RcppResultSet rs;

	rs.add("date", aDate);
	rs.add("dateVec", dateVec);
	rs.add("method", method);
	rs.add("tolerance", tolerance);
	rs.add("maxIter", maxIter);
	rs.add("nlFirstName", nl.getName(0));
	rs.add("nlFirstValue", nl.getValue(0));
	rs.add("matD", matD);
	rs.add("stlvec", stlvec);
	rs.add("stlmat", stlmat);
	rs.add("a", a, nrows, ncols);
	rs.add("v", v, len);
	rs.add("stringVec", stringVec);
	rs.add("strings", svec);
	rs.add("InputDF", inframe);
	rs.add("PreDF", frame);


	// Instead of returning selected input parameters as we did in
	// the last three statements, the entire input parameter list
	// can be returned like this:
	rs.add("params", params, false);


	// Get the list to be returned to R.
	rl = rs.getReturnList();
	
    } catch(std::exception& ex) {
	exceptionMesg = copyMessageToR(ex.what());
    } catch(...) {
	exceptionMesg = copyMessageToR("unknown reason");
    }
    
    if(exceptionMesg != NULL)
	error(exceptionMesg);

    return rl;
*/
}



/*
 * The following class definitions employ advanced features of the Rcpp
 * library and R, permitting the C++ programmer to call user-defined functions
 * on the R side. They should be skipped on first reading.
 */

/*
 * Define a class that can be used to call an R function that expects a
 * real vector argument and returns a scalar. The R function is defined in
 * the example section of the documentation page for RcppExample (see
 * RcppExample.Rd).
 */
class MyRVectorFunc : public RcppFunction {
public:
    MyRVectorFunc(SEXP fn) : RcppFunction(fn) {}

    // This trivial function will use an R function to compute the
    // sum of the elements of v!
    double getSum(vector<double>& v) {

	// Turn vector into a SEXP that can be passed to R as an argument.
	setRVector(v);

	// Call the R function that was passed in as the paramter fn, with
	// the SEXP vector that was just set as its argument.
	SEXP result = vectorCall();

	// Assuming that the R function simply returns a real number we
	// pass it back to the C++ user as follows. If the R function returns
	// something more complicated transform result into a C++ object to
	// be returned, and  clear the part of the protection stack due to
	// this object before returning (to prevent protection stack overflow).
	// Note that it is unsafe to do this if the returned result depends
	// on PROTECT-ed SEXP's. For example, result should not be 
	// wrapped in a class like RcppParams where objects hold onto the
	// the PROTECT-ed SEXP that was used to construct them.

	double value = REAL(result)[0];

	// Safe now to clear the contribution of this function to the
	// protection stack.
	clearProtectionStack();

	return value;
    }
};

/*
 * Define a class that can be used to call an R function that expects a
 * heterogeneous list argument, and returns a vector of the same length
 * with 1 added to each component (no names). The R function is defined in
 * the example section of the documentation page for RcppExample (see
 * RcppExample.Rd).
 */
/*
class MyRListFunc : public RcppFunction {
public:
    MyRListFunc(SEXP fn) : RcppFunction(fn) {}
    vector<double> addOne(double alpha, double beta, double gamma) {

	// Build argument list.
	setRListSize(3);
	appendToRList("alpha", alpha);
	appendToRList("beta",  beta);
	appendToRList("gamma", gamma);

	// Call the R function passed in as fn with the list argument just
	// constructed.
	SEXP result = listCall();

	// Turn returned R vector into a C++ vector, clear protection stack,
	// and return.
	vector<double> vec(length(result));
	for(int i=0; i < length(result); i++)
	    vec[i] = REAL(result)[i];

	// See comments in previous class definition on the purpose of this.
	clearProtectionStack();

	return vec;
    }
};
*/
