/*****************************  parameterloop.h   *****************************
* Author:        Agner Fog
* Date created:  1995-02-10
* Last modified: 2024-10-13
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines parameter loops that let you run multiple simulations
* and search for parameter ranges that lead to a certain result
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#pragma once

// ParameterLoop type values
const int loopUnused = 0;              // loop unused
const int loopSweep  = 1;              // parameter sweep
const int loopSearch = 2;              // search for limiting parameter values
const int loopJob    = 3;              // job loop. not implemented

// ParameterLoop state values
const int loopStateUnused   = 0;       // unused
const int loopStateIdle     = 1;       // finished or about to start
const int loopStateRunning  = 2;       // running
const int loopStateFound1   = 3;       // loopSearch first limit found,
const int loopStateLast     = 4;       // last run
const int loopStateFinished = 5;       // last run finished
const int loopStateAbort    = 6;       // aborted

class Worker;                          // defined in run.cpp
struct AltruData;                      // defined in altuist.h


// record in list of parameters than can be swept
struct SweepParameter {
    int offset;                        // offset into AltruData
    int varType;                       // varInt32 or varFloat
    int arraysize;                     // 1 if scalar, >1 if array
    int bbGeographyParametersUsed;     // enabling bit in bGeographyParametersUsed
    int bbGroupPropertiesUsed;         // enabling bit in bGroupPropertiesUsed
    const char * name;                 // name of parameter
};

extern SweepParameter sweepParameterList[]; // list of parameters than can be swept


// record in list of simulation results used for parameter search
struct ResultSet {
    float par;                         // swept parameter
    float gFraction;                   // fraction of mutant gene
    int32_t generations;               // generations
    int32_t result;                    // simulation result
    int32_t order;                     // simulation number
};
const int resultListLength = 128;      // length of list of simulation results

// Record in list of search results used for 2-dimensional map of parameter limits.
// This list contains limiting x values for each y value
struct ParameterLimits {
    float sValue;                      // value of swept parameter in loop 1
    float limit1;                      // y value limit between altruism and polymorphism
    float limit2;                      // y value limit between polymorphism and egoism
    int result;                        // bit 0: 0 = low y gives altruism, 1 = high y gives altruism (= direction)
                                       // bit 1: results out of order
                                       // bit 2: all died
                                       // bit 3: no limits found. Don't know if direction
};
const int limitListLength = 1024;      // length of list of search results


// define a parameter loop or parameter search
class ParameterLoop {
public:
    ParameterLoop() {                  // constructor
        reset();  worker = 0;
    }
    void reset(void);                  // clear all
    void setup(int iSweep, Worker * worker); // set all parameters
    void init(int iSweep);             // start sweep or search
    void step(void);                   // get next
    void abort(void);                  // abort sweep or search
    void initSweep(void);              // start sweep
    void stepSweep(void);              // get next
    void initSearch(void);             // start search
    void stepSearch(void);             // get next in search
    void checkLimitList(void);         // check for entries in limitlist with unknown direction
    void changeParameters(int iSweep, Worker * worker); // modify parameters on unfinished loop

    int parIndex;                      // index to sweepParameterList
    float increment;                   // for linear sweeps: increment,
                                       // for logarithmic sweeps: steps per decade
    float parameterValue;              // current value
    float logValue;                    // log10(parameterValue) if logarithmic
    float lastValue;                   // parameterValue for the last finished simulation, corresponds with Worker::lastResult
    bool direction;                    // loopSweep: false: low value gives altruism, true: high value gives altruism
    bool logarithmic;                  // true: logarithmic
    bool isFitExpo;                    // variable is group fitness exponent
    int type;                          // 1: parameter sweep, 2: search
    int state;                         // state of loop
    int varType;                       // 3: 32 bit integer, 8: float
    float limit1;                      // limit between altruism and polymorphism (initial guess or result)
    float limit2;                      // limit between polymorphism and egoism (initial guess or result)
    int limitsFound;                   // 1 = limit1 valid, 2 = limit2 valid
    int lastFound;                     // limitsFound from previous search, used for initial guess
    float startValue;                  // start value of sweep parameter
    float endValue;                    // end value of sweep parameter
protected:
    Worker * worker;                   // access functions in worker thread
    AltruData * d;                     // point to common data
    int arraySize;                     // size of array. 1 if scalar
    int varOffset;                     // offset to variable in AltruData
    float logStartValue;               // log10(startValue)
    float logEndValue;                 // log10(endValue)
    int8_t * pSweepVar;                // pointer to sweep variable
    // variables and functions for loopSearch:
    int lastParIndex;                  // parameter in last search
    int iSweep;                        // index into AltruData::sweepType
    int32_t num;                       // run number
    void setValue(float x);            // set model parameter to value
    int regression(int, int, int, float*, float*);// linear regression on resultList
    int checkPoints(float*p, float deltax, int i);// check distance between search points
public:
    int allOutcomes;                   // search results: 1=aborted, 2=egoism found, 4=polymorphism found, 8=altruism found, 16=all died
    float x2p(float x);                // convert search parameter to model parameter
    float p2x(float p);                // convert model parameter to search parameter
};