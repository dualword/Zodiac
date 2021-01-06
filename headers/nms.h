/*   Zodiac - molecular modelling environment. www.zeden.org
 *class implementation of Nelder Mead Simplex Search
 *Adam Gurson College of William & Mary 1999
 *
 * modified slightly by Anne Shepherd (pls), 8/00


 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#ifndef NMS_H
#define NMS_H

#include "optimiser.h"
#include "subscrpt.h"
#include "cmat.h"
#include "vec.h"



class NMS : public Optimiser {
	public:
	NMS ()  : Optimiser () 
	{
	};	
	NMS (Function* f) : Optimiser (f) 
	{
		init_parameters ();

		vector<Variable *> vars = optFunc ->access_variables ();
		init (optFunc ->access_variables ());

	};	
	
	~NMS ()  {
	if(simplex != NULL) delete simplex;
   if(simplexValues != NULL) delete [] simplexValues;
   delete centroid;
   delete reflectionPt;
   delete expansionPt;
   delete contractionPt;
   delete scratch;
   delete scratch2;
   //NOTE: Matrix and Vector classes have their own destructors
	};
	
void init (vector<Variable*>& vars) {
		for (unsigned int i = 0; i < vars.size (); i++) {
			(*scratch)[i] = (double) *(vars[i]->value);
		}
		InitRegularTriangularSimplex(scratch, 2);	
}


void set_iteration_limit (int iter) {
	maxCalls = iter;
}
void init_parameters () {
   dimensions = optFunc ->access_variables ().size ();
   functionCalls = 0;
   simplex = NULL;
   simplexValues = NULL;
   centroid = new Vector<double>(dimensions,0.0);
   reflectionPt = new Vector<double>(dimensions,0.0);
   expansionPt = new Vector<double>(dimensions,0.0);
   contractionPt = new Vector<double>(dimensions,0.0);
   alpha = 1.0;
   beta = 0.5;
   gamma = 2.0;
   sigma = 0.5;
   maxCalls = 100000;
   scratch = new Vector<double>(dimensions,0.0);
   scratch2 = new Vector<double>(dimensions,0.0);
   stoppingStepLength = 0.02;
}

	float run () {
		double secondHighestPtValue; // used for contraction/reflection decision
		toleranceHit = 0;
		
		FindMinMaxIndices();
		do {
			FindCentroid();
			secondHighestPtValue = simplexValues[SecondHighestPtIndex()];
			// reflection step
			FindReflectionPt();
			
			// stop if at maximum function calls and update the simplex
			/*changed  8/8/00 to fix the problem of maxCalls == -1  --pls
			 formerly read if(functionCalls <= maxCalls)
			 */
			if (maxCalls > -1 && functionCalls >= maxCalls) {
				FindMinMaxIndices();
				ReplaceSimplexPoint(maxIndex, *reflectionPt);
				simplexValues[maxIndex] = reflectionPtValue;
				FindMinMaxIndices(); 
				return 0.;
			} // if using call budget
			
			// possibility 1
			if(simplexValues[minIndex] > reflectionPtValue) {
				FindExpansionPt(); // expansion step
				
				if (reflectionPtValue > expansionPtValue) {
					ReplaceSimplexPoint(maxIndex, *expansionPt);
					simplexValues[maxIndex] = expansionPtValue;
				} // inner if
				else {
					ReplaceSimplexPoint(maxIndex, *reflectionPt);
					simplexValues[maxIndex] = reflectionPtValue;
				} // else         
			} // if for possibility 1
			
			// possibility 2
			
			else if( (secondHighestPtValue > reflectionPtValue        ) &&	(   reflectionPtValue >= simplexValues[minIndex]) ) {
				ReplaceSimplexPoint(maxIndex, *reflectionPt);
				simplexValues[maxIndex] = reflectionPtValue;
			} // else if for possibility 2
			
			// possibility 3
			else if( reflectionPtValue >= secondHighestPtValue ) {
				FindContractionPt(); // contraction step
				if(maxPrimePtId == 0) {
					if( contractionPtValue > maxPrimePtValue ) {
						ShrinkSimplex();
					} // inner if
					else {
						ReplaceSimplexPoint(maxIndex, *contractionPt);
						simplexValues[maxIndex] = contractionPtValue;
					} // inner else
				} // maxPrimePtId == 0
				else if(maxPrimePtId == 1) {
					if( contractionPtValue >= maxPrimePtValue ) {
						ShrinkSimplex();
					} // inner if
					else {
						ReplaceSimplexPoint(maxIndex, *contractionPt);
						simplexValues[maxIndex] = contractionPtValue;
					} // inner else
				} // maxPrimePtId == 1
			} // else if for possibility 3
			
			// if we haven't taken care of the current simplex, something's wrong
			else {
				cerr << "Error in ExploratoryMoves() - "
				<< "Unaccounted for case.\nTerminating.\n";
				return 0.;
			}
			FindMinMaxIndices();
		} while (!Stop());   // while stopping criteria is not satisfied
		return 0.f;
	} // ExploratoryMoves()



	void ReplaceSimplexPoint(int index, const Vector<double>& newPoint) 
	{  
		for( int i = 0; i < dimensions; i++ ) {
			(*simplex)[index][i] = newPoint[i];
		} // for;
	}

   void CalculateFunctionValue(int index) {
      *scratch = (*simplex).row(index);
	int success;
	fcnCall(dimensions, (*scratch).begin(), simplexValues[index], success);
	if(!success) cerr<<"Error calculating point in CalculateFunctionValue().\n";
   }
   // finds the f(x) value for the simplex point indexed at index and
   // replaces the proper value in simplexValues

void SetAlpha(double newAlpha)
{
   alpha = newAlpha;
} // SetAlpha()

void SetBeta(double newBeta)
{
   beta = newBeta;
} // SetBeta()

void SetGamma(double newGamma)
{
   gamma = newGamma;
} // SetGamma()

void SetSigma(double newSigma)
{
   sigma = newSigma;
} // SetGamma()
   // these functions allow the user to set values of the reflection,
   // contraction, expansion, and shrinking coefficients
   

bool Stop()
{
   if(maxCalls > -1) {
      if(functionCalls >= maxCalls)
         return true;
   }

   double mean = 0.0;

   for( int i = 0; i <= dimensions; i++) {
      if( i != minIndex ) {
         mean += simplexValues[i];
      } // if
   } //for 

   mean /= (double)dimensions;

   // Test for the suggested Nelder-Mead stopping criteria
   double total = 0.0;
   for( int i = 0; i <= dimensions; i++ ) {
      total += pow( simplexValues[i] - mean ,2);
   } //for
   total /= ((double)dimensions + 1.0);
   total = sqrt(total);
   
   
 //   printSimplex();
   if(total < stoppingStepLength) {
      toleranceHit = 1;
      return true;
   }
   else
      return false;
} // Stop()
   // returns true if the stopping criteria have been satisfied

void fcnCall(int n, double *x, double& f, int& flag)
{
	for (unsigned int i = 0; i < n; i++) {
		vector <Variable*>& func_vars = optFunc ->access_variables();
	//	func_vars[i] ->value = &x[i];
		*func_vars[i] ->value =  x[i];
	}
	f =  (double) optFunc -> evaluate ();
   functionCalls++;
   flag = 1;
} // fcnCall()

// Simplex-altering functions

   // indirection of function call for purposes of keeping an accurate
   // tally of the number of function calls

 // Simplex-altering functions

void InitRegularTriangularSimplex(const Vector<double> *basePoint,
                                            const double edgeLength)
{
  //  This routine constructs a regular simplex (i.e., one in which all of 
  //  the edges are of equal length) following an algorithm given by Jacoby,
  //  Kowalik, and Pizzo in "Iterative Methods for Nonlinear Optimization 
  //  Problems," Prentice-Hall (1972).  This algorithm also appears in 
  //  Spendley, Hext, and Himsworth, "Sequential Application of Simplex 
  //  Designs in Optimisation and Evolutionary Operation," Technometrics, 
  //  Vol. 4, No. 4, November 1962, pages 441--461.

   int i,j;
   double p, q, temp;
   Matrix<double> *plex = new Matrix<double>(dimensions+1,dimensions,0.0);
   for( int col = 0; col < dimensions; col++ ) {
      (*plex)[0][col] = (*basePoint)[col];
   }

   temp = dimensions + 1.0;
   q = ((sqrt(temp) - 1.0) / (dimensions * sqrt(2.0))) * edgeLength;
   p = q + ((1.0 / sqrt(2.0)) * edgeLength);

   for(i = 1; i <= dimensions; i++) { 
      for(j = 0; j <= i-2; j++) {
         (*plex)[i][j] = (*plex)[0][j] + q;
      } // inner for 1
      j = i - 1;
      (*plex)[i][j] = (*plex)[0][j] + p;
      for(j = i; j < dimensions; j++) {
            (*plex)[i][j] = (*plex)[0][j] + q;
      } // inner for 2
   } // outer for

   InitGeneralSimplex(plex);
   delete plex;
} // InitRegularTriangularSimplex()

   // deletes any existing simplex and replaces it with a regular
   // triangular simplex in the following manner:
   //
   // basePoint points to a point that will be the "origin" of the
   //    simplex points (it will be a part of the simplex)
   // edgeLength is the length of each edge of the "triangle"
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: basePoint is assumed to be of proper dimension

void InitFixedLengthRightSimplex(const Vector<double> *basePoint,
                                           const double edgeLength)
{
  // to take advantage of code reuse, this function simply turns
  // edgeLength into a vector of dimensions length, and then
  // calls InitVariableLengthRightSimplex()

   double* edgeLengths = new double[dimensions];
   for( int i = 0; i < dimensions; i++ ) {
      edgeLengths[i] = edgeLength;
   }
   InitVariableLengthRightSimplex(basePoint,edgeLengths);
   delete [] edgeLengths;
} // InitFixedLengthRightSimplex()
   // deletes any existing simplex and replaces it with a right-angle
   // simplex in the following manner:
   //
   // basePoint points to a point that will be the "origin" of the
   //    simplex points (it will be a part of the simplex)
   // edgeLength is to be the length of each simplex side extending
   //    from the basePoint along each positive coordinate direction.
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: basePoint is assumed to be of proper dimension

 void InitVariableLengthRightSimplex(const Vector<double> *basePoint,
                                              const double* edgeLengths)
{
   Matrix<double> *plex = new Matrix<double>(dimensions+1,dimensions,0.0);
   for( int i = 0; i < dimensions; i++ ) {
      // we're building the basePoint component-by-component into
      // the (n+1)st row
      (*plex)[dimensions][i] = (*basePoint)[i];

      // now fill in the ith row with the proper point
      for( int j = 0; j < dimensions; j++ ) {
         (*plex)[i][j] = (*basePoint)[j];
         if( i == j )
            (*plex)[i][j] += edgeLengths[i];
      }
   } // for
   InitGeneralSimplex(plex);
   delete plex;
} // InitVariableLengthRightSimplex()
   // deletes any existing simplex and replaces it with a right-angle
   // simplex in the following manner:
   //
   // basePoint points to a point that will be the "origin" of the
   //    simplex points (it will be a part of the simplex)
   // edgeLengths points to an array of n doubles, where n is the
   //    dimension of the given search. x_1 will then be located
   //    a distance of edgeLengths[0] away from the basepoint along the
   //    the x_1 axis, x_2 is edgeLengths[1] away on the x_2 axis, etc.
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: basePoint and edgeLengths are assumed to be of proper dimension

void InitGeneralSimplex(const Matrix<double> *plex)
{
   functionCalls = 0;
   if( simplex != NULL ) { delete simplex; }
   if( simplexValues != NULL ) { delete [] simplexValues;}
   simplex = new Matrix<double>((*plex));
   simplexValues = new double[dimensions+1];

   int success;
   for( int i = 0; i <= dimensions; i++ ) {
      *scratch = (*plex).row(i);
      fcnCall(dimensions, (*scratch).begin(), simplexValues[i], success);
      if(!success) cerr<<"Error with point #"<<i<<" in initial simplex.\n";
   } // for
   FindMinMaxIndices();
} // InitGeneralSimplex()
   // deletes any existing simplex and replaces it with the one
   // pointed to by plex
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: THIS ASSUMES THAT plex IS OF PROPER DIMENSION

void ReadSimplexFile(istream& fp)
{
   if(fp == NULL) {
      cerr<<"No Input Stream in ReadSimplexFile()!\n";
      return; // There's no file handle!!
   }

   Matrix<double> *plex = new Matrix<double>(dimensions+1,dimensions);
   for( int i = 0; i <= dimensions; i++ ) {
      for ( int j = 0; j < dimensions; j++ ) {
         fp >> (*plex)[i][j];
      } // inner for
   } // outer for
   InitGeneralSimplex(plex);
   delete plex;
} // ReadSimplexFile()
   // may also pass cin as input stream if desired
   // input the values of each trial point
   // NOTE: THIS FUNCTION WILL ONLY ACCEPT n+1 POINTS
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.

 // Query functions

int GetFunctionCalls() const
{
   return functionCalls;
} // GetFunctionCalls()
   // number of objective function evaluations

void GetMinPoint(Vector<double>* &minimum) const
{
   minimum = new Vector<double>((*simplex).row(minIndex));
} // GetMinPoint()
   // simplex point which generates the best objective function
   // value found thus far
   // USER SHOULD PASS JUST A NULL POINTER, WITHOUT PREALLOCATED MEMORY

double GetMinVal() const
{
   return simplexValues[minIndex];
} // GetMinVal()
   // best objective function value found thus far

void GetCurrentSimplex(Matrix<double>* &plex) const
{
   plex = new Matrix<double>((*simplex));
} // GetCurrentSimplex()
   // performs a deep copy of the simplex to a Matrix pointer
   // points to a newly allocated chunk of memory upon return
   // USER SHOULD PASS JUST A NULL POINTER, WITHOUT PREALLOCATED MEMORY

void GetCurrentSimplexValues(double* &simValues) const
{
   simValues = new double[dimensions+1];
   for( int i = 0; i <= dimensions; i++ ) {
      simValues[i] = simplexValues[i];
   } // for
} // GetCurrentSimplexValues()
   // performs a deep copy of the simplexValues array to a double pointer
   // points to a newly allocated chunk of memory upon return
   // USER SHOULD PASS JUST A NULL POINTER, WITHOUT PREALLOCATED MEMORY

int GetVarNo() const
{
   return dimensions;
} // GetVarNo()
   // returns the dimension of the problem

int GetTolHit() const
{
   return toleranceHit;
} // GetTolHit()

   // returns toleranceHit

void printSimplex() const
{
  for( int i = 0; i <= dimensions; i++ ) {
     cout << "   Point:";
     for ( int j = 0; j < dimensions; j++ ) {
        cout << (*simplex)[i][j] << "\t";
     } // inner for
     cout << "Value:" << simplexValues[i] << "\n";
  } // outer for

  cout << "\nFCalls: " << functionCalls << endl << endl;
}
   // prints out the simplex points by row, their corresponding f(x)
   // values, and the number of function calls thus far

 private:

void FindMinMaxIndices()
{
   if(simplexValues == NULL) {
      cerr << "Error in FindMinMaxIndices() - "
           << "The vector of simplexValues is NULL!!\n";
      return;
   }
   minIndex = 0;
   maxIndex = dimensions;
   double min = simplexValues[0];
   double max = simplexValues[dimensions];
   for( int i = 1; i <= dimensions; i++ ) {
      if( simplexValues[i] < min ) {
         min = simplexValues[i];
         minIndex = i;
      } // if
      if( simplexValues[dimensions-i] > max ) {
         max = simplexValues[dimensions-i];
         maxIndex = dimensions - i;
      } // if
   } // for
} // FindMinMaxIndices()

   // sets minIndex to the simplex index of the point which generates
   // the lowest value of f(x)
   // sets maxIndex to the simplex index of the point which generates
   // the highest value of f(x)

int SecondHighestPtIndex()
{
   if(simplexValues == NULL) {
      cerr << "Error in SecondHighestPtValue() - "
           << "The vector of simplexValues is NULL!!\n";
      return -1;
   }
   int secondMaxIndex = minIndex;
   double secondMax = simplexValues[minIndex];
   for( int i = 0; i <= dimensions; i++ ) {
      if(i != maxIndex) {
         if( simplexValues[i] > secondMax ) {
            secondMax = simplexValues[i];
            secondMaxIndex = i;
         } // inner if
      } // outer if
   } // for
   return secondMaxIndex;
} // SecondHighestPtValue()

   // returns simplex index of the point which
   // generates the second highest value of f(x)

void FindCentroid()
{
   (*centroid) = 0.0;
   for( int i = 0; i <= dimensions; i++ ) {
      if( i != maxIndex ) {
         (*centroid) = (*centroid) + (*simplex).row(i);
      } // if
   } // for
   (*centroid) = (*centroid) * ( 1.0 / (double)dimensions );
} // FindCentroid()

   // finds the centroid

void FindReflectionPt()
{ 
   (*reflectionPt) = 0.0;
   (*reflectionPt) = ( (*centroid) * (1.0 + alpha) ) -
                     ( alpha * (*simplex).row(maxIndex) );
   int success;
   fcnCall(dimensions, (*reflectionPt).begin(), reflectionPtValue, success);
   if(!success) {
      cerr << "Error finding f(x) for reflection point at"
           << "function call #" << functionCalls << ".\n";
   } // if
} // FindReflectionPt()

   // finds the reflection point and sets its f(x) value

void FindExpansionPt()
{
   (*expansionPt) = 0.0;
   (*expansionPt) = ( (*centroid) * (1.0 - gamma) ) +
                    ( gamma * (*reflectionPt) );
   int success;
   fcnCall(dimensions, (*expansionPt).begin(), expansionPtValue, success);
   if(!success) {
      cerr << "Error finding f(x) for expansion point at"
           << "function call #" << functionCalls << ".\n";
   } // if
} // FindExpansionPt()
   // finds the expansion point and sets its f(x) value

void FindContractionPt()
{
   // need to first define maxPrimePt
   Vector<double> *maxPrimePt = scratch;
   if(simplexValues[maxIndex] <= reflectionPtValue) {
      *maxPrimePt = (*simplex).row(maxIndex);
      maxPrimePtValue = simplexValues[maxIndex];
      maxPrimePtId = 1;
   } // if
   else {
      maxPrimePt = reflectionPt;
      maxPrimePtValue = reflectionPtValue;
      maxPrimePtId = 0;
   } // else

   (*contractionPt) = ( (*centroid) * (1.0 - beta) ) +
                      ( beta * (*maxPrimePt) );
   int success;
   fcnCall(dimensions, (*contractionPt).begin(), contractionPtValue, success);
   if(!success) {
      cerr << "Error finding f(x) for contraction point at"
           << "function call #" << functionCalls << ".\n";
   } // if
} // FindContractionPt()

   // finds the contraction point and sets its f(x) value

void ShrinkSimplex()
{
   // stop if at maximum function calls
  // changed 5/01 to reflect maxcalls = -1 possibility ---pls
  if ( (maxCalls != (-1)) 
       && (functionCalls >= maxCalls) ) {return;}

   Vector<double> *lowestPt = scratch;
   *lowestPt = (*simplex).row(minIndex);
   Vector<double> *tempPt = scratch2;
   int success;
   for( int i = 0; i <= dimensions; i++ ) {
      if( i != minIndex ) {
         *tempPt = (*simplex).row(i);
         (*tempPt) = (*tempPt) + ( sigma * ( (*lowestPt)-(*tempPt) ) );
         for( int j = 0; j < dimensions; j++ ) {
            (*simplex)[i][j] = (*tempPt)[j];
         } // inner for
         fcnCall(dimensions,(*tempPt).begin(),simplexValues[i],success);
         if (!success) cerr << "Error shrinking the simplex.\n";
         
         // stop if at maximum function calls 
	 // changed 5/01 to reflect maxcalls = -1 possibility ---pls
	 if ( (maxCalls != (-1)) 
	      && (functionCalls >= maxCalls) ) {return;}
	 
      } // if
   } // outer for
} // ShrinkSimplex()

   // this function goes through the simplex and reduces the
   // lengths of the edges adjacent to the best vertex


   int dimensions;                // the number of dimensions
                                  //    (the dimension of the problem)
   Matrix<double> *simplex;       // the current simplex
   double *simplexValues;         // their corresponding f(x) values
   double alpha;                  // reflection coefficient
   double beta;                   // contraction coefficient
   double gamma;                  // expansion coefficient
   double sigma;                  // shrinking coefficient
   int maxCalls;
   double stoppingStepLength;
   int minIndex;                  // index of point generating min f(x)
   int maxIndex;                  // index of point generating max f(x)
   Vector<double> *centroid;      // the current centroid
   Vector<double> *reflectionPt;  // the reflection point
   double reflectionPtValue;      // the value of f(reflectionPt)
   Vector<double> *expansionPt;   // the expansion point
   double expansionPtValue;       // the value of f(expansionPt)
   Vector<double> *contractionPt; // the contraction point
   double contractionPtValue;     // the value of f(contractionPt)
   double maxPrimePtValue;        // min(f(maxIndexPoint),reflectionPtValue)
   long functionCalls;            // tally of number of function calls
   int toleranceHit;              // 1 if stop due to tolerance, 0 if funcCalls
   int maxPrimePtId;              // set by FindContractionPt() and used in
                                  // ExploratoryMoves() to branch in pos. 3

   // the following vectors are simply extra storage space
   // that are used by functions that require a vector of
   // size = dimensions
   Vector<double> *scratch, *scratch2;



};

#endif //NMS_H