/**************************  parameterloop.cpp   ******************************
* Author:        Agner Fog
* Date created:  1995-02-10
* Last modified: 2024-06-22
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file controls parameter loops that let you run multiple simulations
* and search for parameter ranges that lead to a certain result
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

// list of parameters that can be used in loop
/*
struct SweepParameter {
    int offset;                        // offset into AltruData
    int varType;                       // varInt32 or varFloat
    int arraysize;                     // 1 if scalar, >1 if array
    int bbGeographyParametersUsed      // enabling bit in bGeographyParametersUsed
    int bbGroupPropertiesUsed          // enabling bit in bGroupPropertiesUsed
    const char * name;                 // name of parameter
};*/

SweepParameter sweepParameterList[] = {
    {altruDataOffset(seed), varInt32, 1, -1, -1, "random seed"},
    {altruDataOffset(migrationRate[0]), varFloat, 2, 0x1000, 0, "migration rate"},
    {altruDataOffset(fitExpo), varFloat, 1, 0, 0x1000, "group fitness exponent"},
    {altruDataOffset(fit[2]), varFloat, 2, -1, -1, "altruist fitness"},
    {altruDataOffset(fit[0]), varFloat, 8, 0, -1, "individual fitness ratio"},  // fitness ratio of altruists over egoists
    {altruDataOffset(extinctionRate[0]), varFloat, 2, 0, 2, "extinction rate for egoist groups"},
    {altruDataOffset(warIntensity), varFloat, 1, 0, 0x80, "war intensity"},
    {altruDataOffset(nMaxPerDeme), varInt32, 1, 0x04, 0, "max group size"},    
    {altruDataOffset(carryingCapacity), varFloat, 2, 0x10, 0, "carrying capacity"},
    {altruDataOffset(colonySize), varInt32, 1, 0x200, 0, "colony size"},
    {altruDataOffset(surviv), varFloat, 1, 0, 0x10, "survival rate"},
    {altruDataOffset(haystackPeriod), varInt32, 1, 0, 0x100, "haystack period"},
    {altruDataOffset(fit[0]), varFloat, 9, 0, 0x1140, "fitness among altruists/egoists"},   // fitness ratio of all among altruists over egoists
    // the arraysize is actually 8, but the value 9 is used for distinguishing this from "individual fitness ratio" above
    {altruDataOffset(leaderAdvantage), varFloat, 1, 0, 0x10000, "leader advantage"},    
    {0, 0, 0, 0}                       // list must end with zeroes
};

const uint32_t sweepParlistLength = ARR_LEN(sweepParameterList) - 1;  // length of list


void ParameterLoop::reset(void) {
    // clear all
    type = 0;  state = 0;  varType = 0;  num = 0;
    parIndex = 0;  varOffset = 0;
    direction = false; logarithmic = false;
    startValue = 0;  endValue = 0;  increment = 0;
    pSweepVar = 0;  parameterValue = 0;  logValue = 0;  lastValue = 0;
    limit1 = 0;  limit2 = 0;
    limitsFound = 0;  lastFound = 0;  allOutcomes = 0;
}

    
void ParameterLoop::setup(int iSweep, Worker * worker) {
    // set all parameters
    // iSweep is index into AltruData::sweepType
    this->iSweep = iSweep;
    this->worker = worker;
    this->d = worker->d;
    if (iSweep >= maxSweeps) {
        errors.reportError("sweep index out of range");  
        return;
    }
    parIndex = d->sweepParameter[iSweep]; // index into sweepParameterList
    if ((unsigned)parIndex >= sweepParlistLength) {
        errors.reportError("parameter index out of range");  
        return;
    }
    state = loopStateIdle;
    switch (d->sweepType[iSweep]) {
    case 0:  // unused
        type = 0;
        state = loopStateUnused;
        break;
    case 1:  // linear sweep
        type = 1;  logarithmic = false;
        break;
    case 2:  // logarithmic sweep
        type = 1;  logarithmic = true;
        break;
    case 3:  // linear search
        type = 2;  logarithmic = false;
        break;
    case 4:  // logarithmic search
        type = 2;  logarithmic = true;
        break;
    default:
        type = 0;
    }
    varType = sweepParameterList[parIndex].varType;
    varOffset = sweepParameterList[parIndex].offset;
    pSweepVar = (int8_t*)d + varOffset;
    arraySize = sweepParameterList[parIndex].arraysize;
    if (arraySize < 1) arraySize = 1;
    startValue = d->sweepStartValue[iSweep];
    endValue = d->sweepEndValue[iSweep];
    direction = startValue > endValue;
    increment = d->sweepStep[iSweep];
    if (type == 1 && varType <= varInt64 && fabs(increment) < 1.f) {
        increment = 1.f; // integer increment must be at least 1
        d->sweepStep[iSweep] = increment;
    }
    if (direction) increment = -fabs(increment);
    isFitExpo = varOffset == altruDataOffset(fitExpo);  // fitExpo requires special treatment
    // check if valid
    if (type == 2 && varType < varFloat) {
        errors.reportError("Search not possible for integer variable");
        type = 0;    
    }
    // don't use previous guesses
    lastFound = 0; lastParIndex = 0;
}

void ParameterLoop::setValue(float x) {
    // set model parameter to value
    int i;  // loop counter

    switch (varType) {
    case varInt32:
        for (i = 0; i < arraySize; i++) {        
            *((int32_t*)pSweepVar + i) = (int32_t)lround(x);
        }
        break;
    case varFloat:
        if (varOffset == altruDataOffset(fit[0])) {
            if (arraySize == 8) {
                // special case: fitness ratio of altruist/egoists
                float * p = (float*)pSweepVar;
                if (iSweep == 1 && worker->loops[0].varOffset == varOffset && worker->loops[0].arraySize == 9) {
                    // double special case. Nested loops both use fit[] array
                    float x2 = worker->loops[0].startValue;
                    p[1] = p[0] * x2;  p[5] = p[4] * x2;
                }
                p[2] = p[0] * x;
                p[3] = p[1] * x;
                p[6] = p[4] * x;
                p[7] = p[5] * x;
                break;
            }
            if (arraySize == 9) {
                // special case: fitness ratio of all among altruist/egoists
                float * p = (float*)pSweepVar;
                if (iSweep == 1 && worker->loops[0].varOffset == varOffset && worker->loops[0].arraySize == 8) {
                    // double special case. Nested loops both use fit[] array
                    float x2 = worker->loops[0].startValue;
                    p[2] = p[0] * x2;  p[6] = p[4] * x2;
                }
                p[1] = p[0] * x;
                p[3] = p[2] * x;
                p[5] = p[4] * x;
                p[7] = p[6] * x;
                break;
            }
        }
        // all other cases
        for (i = 0; i < arraySize; i++) {
            *((float*)pSweepVar + i) = x;
        }
        break;
    default: // other types not supported yet
        break;
    }
    if (isFitExpo) { // variable is group fitness exponent
        if (x > 0.999f && x < 1.001f) {
             d->fitfunc = 2;
             *((float*)pSweepVar) = 1.0f;
        }
        else if (x <= 1.E-6f) d->fitfunc = 0;
        else if (x < 1.f) d->fitfunc = 1;
        else if (x <= 1.E6f) d->fitfunc = 3;
        else d->fitfunc = 4;
    }
}

void ParameterLoop::changeParameters(int iSweep, Worker * worker) {
    // modify parameters on unfinished loop
    setup(iSweep, worker);
    if (((startValue > parameterValue) ^ direction) != 0) {
        parameterValue = startValue;
        setValue(parameterValue);    
    }
}

void ParameterLoop::init(int iSweep) {
    // start sweep or search
    if (state == loopStateUnused) return;
    state = loopStateRunning;
    lastValue = 0.f;
    this->iSweep = iSweep;
    switch (type) {     
    case 1:
        initSweep();  break;
    case 2:
        initSearch();  break;
    }
    //worker->waitForScreen();
    d->requestUpdate = true;
}

void ParameterLoop::step(void) {
    // increment loop parameter
    if (state == loopStateUnused) return;
    switch (type) {     
    case 1:
        stepSweep();  break;
    case 2:
        stepSearch();  break;
    }
}

void ParameterLoop::abort(void) { 
    if (state == loopStateUnused) return;
    // abort sweep
    state = loopStateAbort;
    worker->loopsPaused = false;
    worker->doWorkSlot(state_stop);
    d->requestUpdate = true;
}

void ParameterLoop::initSweep(void) {
    // set parameter for start or restart of sweep
    parameterValue = logValue = startValue;
    lastValue = 0.f;
    if (logarithmic) {
        if (parameterValue > 0 && startValue > 0 && endValue > 0) {
            logValue = log10f(parameterValue);
            logStartValue = log10f(startValue);
            logEndValue = log10f(endValue);
        }
        else {
            errors.reportError("logarithm of non-positive value");
        }
    }
    if (type == 1 && iSweep == 1 && d->sweepType[0] >= 3) {
        // make x-y map of limiting parameter values
        // This ParameterLoop sweeps the y values, while an 
        // underlying search finds the limiting x values
        worker->limitList[0].sValue = parameterValue;   
        worker->limitListNum = 0;
        //d->graphicsType = graphicsLimits;
    }
    num = 0;  allOutcomes = 0;
    d->nSum = 0;
    for (int j = 0; j < d->nLoci; j++) d->fgsum[j] = 0;
    setValue(parameterValue);
    //num++;  // count number of elements in limitList
}


void ParameterLoop::stepSweep(void) {
    // increment loop parameter
    if (state == loopStateLast || state == loopStateFinished) {
        state = loopStateFinished;               // this loop is finished
        return;    
    }
    float lastValue = parameterValue;
    // get next
    if (!logarithmic) {
        if (!direction) { // ascending
            parameterValue += increment;
            if (parameterValue >= endValue - 0.01f * increment) {
                parameterValue = endValue;
                state = loopStateLast;
                //return;
            }
        }
        else { // direction
            parameterValue -= fabs(increment);
            if (parameterValue <= endValue + 0.01f * increment) {
                parameterValue = endValue;
                state = loopStateLast;
                //return;
            }
        }
    }
    else {
        if (!direction) { // ascending
            logValue += 1.f/increment;
            parameterValue = powf(10.f, logValue);
            if (logValue >= logEndValue - 0.01f/increment) {
                parameterValue = endValue;  logValue = logEndValue;
                state = loopStateLast;
                //return;
            }
        }
        else { // direction
            logValue -= 1.f/increment;
            parameterValue = powf(10.f, logValue);
            if (logValue <= logEndValue + 0.01f/increment) {
                parameterValue = endValue;  logValue = logEndValue;
                state = loopStateLast;
                //return;
            }
        }
    }
    if (state != loopStateFinished) {
        if (parameterValue == lastValue) {                 // make sure new value is not rounded to same as last value
            if (varType <= varInt64) {
                // integer variable has been rounded to same value as before. make sure it is incremented
                parameterValue += direction ? -1.f : 1.f;
                if (parameterValue > 0.f) {
                    logValue = log10(parameterValue);
                }
                else if (logarithmic) state = loopStateAbort; // logarithmic sweep must end if value is zero
            }
            else state = loopStateAbort;                   // unexpected error
        }
        setValue(parameterValue);                          // set next value
    }
}

void Worker::doMultipleSimulations() {
    // start or resume parameter sweeps or search involving multiple simulations
    int  iLoop;                                            // index for up to three loops
    bool finished = false;
    if (d->sweepState <= state_start) {    
        lastResult = resultUnknown;
    }
    if (loopsPaused) {
        if (d->sweepState == state_run || d->sweepState == state_single_step) {
            // resuming from pause
            loopsPaused = false;
        }
    }
    else if (d->sweepState == state_start || d->sweepState == state_run) {
        // initialize all loops
        for (iLoop = maxSweeps - 1; iLoop >= 0; iLoop--) {
            loops[iLoop].setup(iLoop, this);               // set up parameters
            loops[iLoop].init(iLoop);                      // set first value
            lastParameter[iLoop] = loops[iLoop].parameterValue; // save parameters for file output
        }
        if (d->sweepState == state_start) {
            // do first simulation before incrementing parameters
            d->sweepState = state_run;
            waitForScreen();                               // wait for screen to be updated before making new islands
            doWholeSimulation();                           // do simulation with first value of all parameters
            d->timeUsedLoops = d->timeUsed;                // add time used for this simulation        
            int i = 0;    
            //waitForScreen();
            d->requestUpdate = true;
            emit resultReadySignal(state_run);             // update status line
        }
    }

    // run all loops, innermost first
    while (!finished) {
        // check for stop or pause
        switch (d->sweepState) {
        case state_pause: case state_pause_next: 
            loopsPaused = true;
            d->sweepState = state_pause;
            emit resultReadySignal(state_pause);
            //return;
            while (d->sweepState == state_pause) {
                thread()->msleep(10);  // wait for state to change          
            }
            loopsPaused = false;
            break;

        case state_single_step:
            // pause after next simulation
            d->sweepState = state_pause_next;
            break;

        case state_stop:  case state_error:
            loopsPaused = false;
            d->sweepState = d->runState = state_stop;
            emit resultReadySignal(state_stop);
            return;
        }

        if (d->parametersChanged) {
            // user has changed parameters during run
            for (iLoop = 0; iLoop < maxSweeps; iLoop++) {
                loops[iLoop].changeParameters(iLoop, this);
            }
            d->parametersChanged = false;
        }

        loopsFinished = 0;
        // increment one or more loops
        for (iLoop = 0; iLoop < maxSweeps; iLoop++) {
            // save last parameters for file output
            lastParameter[iLoop] = loops[iLoop].parameterValue;

            // increment parameter
            loops[iLoop].step();

            if (loops[iLoop].state >= loopStateRunning && loops[iLoop].state < loopStateFinished) {
                // this loop is running and not finished
                // restart all lower loops
                for (int lower = 0; lower < iLoop; lower++) {
                    loops[lower].init(lower);              // restart lower loop, set first value
                }
                break;  // this loop is not finished. Don't increment higher loop
            }
            // this loop is finished or unused
            loopsFinished |= 1 << iLoop;
            if (iLoop == maxSweeps - 1) {
                finished = true;                           // highest loop finished
            }
        }
        // write file output
        if (loopsFinished) fileOutSweepFinished();         // at least one sweep finished
        else fileOutSweepNext();                           // write last parameters

        if (finished) break;                               // all sweeps finished

        waitForScreen();                                   // wait for screen to be updated before making new islands
        doWholeSimulation();                               // do simulation with new parameter values
        d->timeUsedLoops += d->timeUsed;                   // add time used for this simulation
        
        // steady state statistics
        uint32_t nSteady = d->nSteady;                     // number of generations before steady state
        if (!(d->bOutOptions & 0x1000)) {                  // steady state not specified by user
            nSteady = d->maximumGenerations / 4u;        
        }
        if (d->generations > nSteady) {
            for (int j = 0; j < d->nLoci; j++) {
                if (d->locusUsed[j] && d->totalPopulation > 0) {                
                    d->fgsum[j] += float(d->genePool[j][1]) / float(d->totalPopulation * 2);
                }
            }
            d->nSum++;
        }
        d->requestUpdate = true;
        emit resultReadySignal(d->sweepState);             // update status line
    }

    // tell GUI thread that we have finished all loops
    d->sweepState = state_stop;
    d->runState =  state_stop;
    loopsPaused = false;
    d->requestUpdate = true;
    emit resultReadySignal(state_stop);
}


// Functions for parameter search:

float ParameterLoop::x2p(float x) {
    // convert search parameter to model parameter.
    // The search parameter is linearized for use in regression and extrapolation
    float p = x;
    if (logarithmic) p = powf(10.f, p);
    if (varType <= varInt32) p = round(p);
    return p;
}

float ParameterLoop::p2x(float p) {
    // convert model parameter to search parameter
    float x = p;
    if (logarithmic) {
        if (p < 1.E-6f) x = -6.f;                // zero or close to zero
        else if (p > 1.E6f) x = 6.f;             // infinity or very high
        else x = log10(p);
    }
    return x;
}

int ParameterLoop::regression(int i1, int i2, int returnType, float* intercept, float* slope) {
    // linear regression on resultList.
    // performs a linear regression of some result parameter in the list results versus x.
    // i1 and i2 are indexes to the first and last parameter set to include.
    // returnType indicates the dependent parameter:
    // returnType = 0:  gFraction, 1:  generations
    // return value = 1 if success, 0 if slope is near infinite

    float x, y;                                  // x and y variables
    float sx, sxx, sy, sxy;                      // sums and square sums
    int i;                                       // index into resultList
    int n;                                       // length of resultList
    ResultSet * resultList = worker->resultList; // list of ResultSet records

    *intercept = *slope = 0;
    n = i2 - i1 + 1;
    if (n < 2) return 0;
    // compute sums
    sx = sxx = sy = sxy = y = 0.f;
    for (i = i1; i <= i2; i++) {
        x = p2x(resultList[i].par);
        switch (returnType) {
        case 0:  y = resultList[i].gFraction;  break;
        case 1:  y = resultList[i].generations;  break;
        }
        sx += x; sxx += x * x; sy += y; sxy += x * y;
    }
    y = sxy - sx * sy / n;  // SAP
    x = sxx - sx * sx / n;  // SAKx
    if (fabs(x) > 1.E-20f) {
        *slope = y / x;
        *intercept = (sy - *slope * sx) / n;
        return 1;
    }
    else {
        return 0;
    }
}

int ParameterLoop::checkPoints(float*p, float deltax, int i) {
    // check distance between search points
    // checks the parameter value supposed to lie between points # i and i+1.
    // The minimum distance from these points is deltax (in x domain).
    // return value:  0:  OK
    //                1:  distance < deltax.  p modified to make distance = deltax
    //                2:  distance between points i and i+1 <= 2*deltax. finished
    //                3:  p < point # i  or  p > point i+1.
    //                4:  p < startValue or  p > endValue. p changed to limit

    float a, b, ax, bx, x;
    a = worker->resultList[i].par;  
    b = worker->resultList[i + 1].par;
    ax = p2x(a);
    bx = p2x(b);
    x  = p2x(*p);

    if (*p < startValue || *p > endValue) {
        // out of range
        *p = x2p((ax + bx) / 2);
        if (*p < startValue) *p = startValue;
        if (*p > endValue)   *p = endValue;
        return 4;
    }
    if (x < ax || x > bx) {
        return 3;
    }
    if (fabs(bx - ax) <= 2 * deltax || (!logarithmic && fabs(b - a) <= 2 * deltax)) {
        // distance between points is small. accuracy reached
        *p = x2p((ax + bx) / 2);
        if (*p == a) (*p)++;
        return 2;
    }
    if (x < ax + deltax) {
        *p = x2p(ax + deltax);
        return 1;
    }
    if (x > bx - deltax) {
        *p = x2p(bx - deltax);
        return 1;
    }
    return 0;
}

// functions used in paremeter search
void ParameterLoop::initSearch(void) {
    // set parameter for start or restart of search
    initSweep();  // most initializations are the same as for sweep      
    limitsFound = 0;
    if (d->maximumGenerations == 0)  d->maximumGenerations = 1000;
    if (d->stopCriterionDegree == 0) d->stopCriterionDegree = 0.99f;
    if (startValue == endValue) {
        state = loopStateAbort;  return;
    }
    if (startValue > endValue) {                           // make sure startValue < endValue
        float temp = startValue;  startValue = endValue;  endValue = temp;
    }
    if (logarithmic) {
        if (startValue < 0.f) {
            state = loopStateAbort;  return;
        }
        if (startValue == 0.f && endValue >= 1.E6f) {      // 0 to infinity
            startValue = 1.E-3f;  endValue = 1.E3f;
        }
        if (startValue == 0.f) startValue = 1.E-4f * endValue;
        if (endValue >= 1.E6f) endValue = 1.E4f * startValue;
        increment = 1.f / ((logEndValue - logStartValue) * d->sweepStep[0]);
        if (d->sweepStep[0] < 10.f) {        
            increment = (logEndValue - logStartValue) * 0.001f;
        }
    }
    else {  // linear
        if (endValue >= 1.E6f) {
            if (startValue == 0.f) endValue = 1.E4f; else endValue = 1000.f * startValue;
        }
        if (increment < 0) increment = -increment;
        if (increment == 0) increment = (endValue - startValue) * 0.001;
    }
    setValue(startValue);

    // make sure steady state statistics are calculated
    if (d->nSteady == 0 || d->nSteady >= d->maximumGenerations) {
        d->nSteady = (uint32_t)d->maximumGenerations * 2u / 3u;
    }
    // check if initial guesses are available
    if (lastParIndex != parIndex ||
        limit1 <= startValue || limit1 >= endValue ||
        limit2 <= startValue || limit2 >= endValue) {
        lastFound = 0;
    }
    lastParIndex = parIndex;
}

void ParameterLoop::stepSearch(void) {
    // get next in search
    // This function is searching for limits between the parameter values that 
    // lead to fixation of wild gene, stable polymorphism, or fixation of mutant gene

    int   i, j;                                  // loop counters or temporary values
    int   h;                                     // temporary index
    int   g0, g1, g2, g3;                        // number of results of each type
    int   ord1 = 0, ord2 = 0;                    // position of out-of-order points in resultList
    float p;                                     // value of the search parameter
    float x;                                     // x = log(p) if logarithmic scale, otherwise x = p
    float xmin, xmax;                            // limits for x
    float fa;                                    // fraction of mutant gene
    float a, b;                                  // interval in which a limit is searched for, in x domain
    float deltax;                                // minimum difference between x points
    float dd;                                    // increment in search
    float slope, intercept;                      // slope and intercept in linear regression
    bool  died = false;                          // all died

    ResultSet * resultList = worker->resultList; // list of simulation results
    if (state >= loopStateFinished) goto FINISHED; // aborted
    allOutcomes |= 1 << worker->lastResult;      // register outcome of last simulation

    p = lastValue;

    if (p < startValue) errors.reportError("p<start"); // should never occur

    if (state == loopStateFinished) {
        checkLimitList();                        // check for entries in limitlist with unknown direction
        return;    
    }

    // get fraction of altruism (mutant) gene
    if (d->nSum > 1 && d->stopCause == stopCauseMaxGen) {
        // use steady state value of gene fraction if stop criterion not reached
        fa = float(d->fgsum[1]) / float(d->nSum);
    }
    else {
        fa = d->geneFraction[0];                 // use final gene fraction if stop criterion reached
    }

    // put simulation result in list ordered by parameter p
    for (i = 0; i < num; i++) if (p < resultList[i].par) break;
    if (i < num) {
        memmove(resultList + i + 1, resultList + i, (num - i) * (int)sizeof(ResultSet));
    }

    resultList[i].par = p;
    resultList[i].gFraction = fa;
    resultList[i].generations = int32_t(d->generations);
    resultList[i].result = worker->lastResult;
    resultList[i].order = ++num;

    if (worker->lastResult >= resultDied) {      // all died or simulation aborted
        died = true;
        goto FINISHED;
    }

    if (num >= resultListLength) {               // result list full
        goto FINISHED;
    }

    // have we had p=EndValue?
    if (num == 1) {
        p = endValue;  goto NEXTP;               // go with endValue
    }

    // direction ascending or descending?
    i = resultList[0].result;
    j = resultList[num-1].result;
    if (i == j) {  // all simulations give the same reult. stop searching
        goto FINISHED;
    }
    direction = i < j;                          // true if result is a descenting function of p
    xmin = p2x(startValue);
    xmax = p2x(endValue);
    deltax = (xmax - xmin) / 1000.f;             // minimum difference between x points

    if (!direction) {
        if (i < resultAltruism && j < resultAltruism) {
            state = loopStateFound1;             // altruism fixation not possible
            limit1 = startValue;                 // limit between altruism and polymorphism
        }
    }
    else {                                       // direction
        if (i > resultEgoism && j > resultEgoism) {
            state = loopStateFound1;             // egoism fixation not possible
            limit1 = startValue;                 // limit egoism and polymorphism
        }    
    }

    // have we had guess 1 ?
    if (lastFound & 1) {                         // bit 0: limit1 found, bit 1: limit2 found
        a = limit1;                              // limit1 from previous simulation
        lastFound &= ~1;                         // don't use the same value again
        i = checkPoints(&a, deltax, 0);
        if (i < 3) {
            p = a;                               // try this value
            goto NEXTP;
        }
    }

    // have we had guess 2 ?
    if ((lastFound & 2) && limit1 != limit2) {
        a = limit2;                              // limit2 from previous simulation
        lastFound &= ~2;                         // don't use the same value again
        i = checkPoints(&a, deltax, 0);
        if (i == 3 && num > 2) i = checkPoints(&a, deltax, 1);
        if (i < 3) {
            p = a;                               // try this value
            goto NEXTP; 
        }
    }

    // count points of each type and check if they are in order
    g0 = g1 = g2 = g3 = 0;
    ord1 = ord2 = 0;
    for (i = 0; i < num; i++) {
        int result2 = resultList[i].result;      // result from list
        if (direction) {                        // swap results if direction
            if (result2 == resultAltruism) result2 = resultEgoism;
            else if (result2 == resultEgoism) result2 = resultAltruism;
        }
        switch (result2) {
        case resultAltruism:                     // altruism
            g0++;
            if (g1 | g2 | g3) ord1 = i;
            break;
        case resultPolymorphism:                 // polymorphism
            if (direction) {
                if (resultList[i].gFraction <= 0.9f) g1++; else g2++;
            }
            else {
                if (resultList[i].gFraction >= 0.1f) g1++; else g2++;
            }
            if (g3) ord2 = i;
            break;
        case resultEgoism:                       // egoism
            g3++;  break;
        default:
            break;
        }
    }

    // search for limit1 between altruism (mutant gene) and polymorphism
    if (ord1) {
        // results out of order. Set limit1 to the middle between the last two points
        limit1 = x2p((p2x(resultList[ord1-1].par) + p2x(resultList[ord1].par)) / 2);
        state = loopStateFound1;  limitsFound |= 1;         
    }
    else {
        // get last altruism point and first non-altruism point, or opposite if direction
        if (g0) {
            a = p2x(resultList[g0 - 1].par);
            b = p2x(resultList[g0].par);
            if (b - a <= deltax) {               // found limit1
                limit1 = x2p((a + b) / 2.);
                limitsFound |= 1;
                state = loopStateFound1;
                if (g1 + g2 == 0) {
                    // found double limit
                    limit2 = limit1;
                    limitsFound |= 3;
                    goto FINISHED;
                }
            }
        }
        else {
            // no altruism
            state = loopStateFound1;
        }
    }

    // search for polymorphism points
    if (state < loopStateFound1) {
        if (g1 + g2 == 0) {
            // no polymorphism points
            if (b - a > (xmax - xmin) / 10.) {
                // try to find polymorphism point
                x = (a + b) / 2;  p = x2p(x);
                checkPoints(&p, deltax, g0 - 1);  
                goto NEXTP;
            }
        }
        else {
            if (g1 == 0 || g1 + g2 < 2) {
                // not enough polymorphism points
                dd = (b - a) / 3;  if (dd < deltax) dd = deltax;
                x = a + dd;  p = x2p(x); 
                checkPoints(&p, deltax, g0 - 1);
                goto NEXTP;
            }
            // enough polymorphism points for linear regression
            i = g1;  if (i < 2) i = 2;  if (i > 4) i = 4;
            // linear regression on polymorphism points
            j = regression(g0, g0 + i - 1, 0, &intercept, &slope);
            if (j == 0 || slope == 0) {
                // regression error
                x = (a + b) * 0.5f;
            }
            else {
                if (direction) x = -intercept / slope;    // extrapolate to gf = 0
                else x = (1.f - intercept) / slope;        // extrapolate to gf = 1
                // check if extrapolation result is useful
                if (x <= a || x >= b) x = (a + b) / 2;
                if (g1 + g2 < 5) dd = (b - a) / 8.f; else dd = (b - a) / 4.f;
                if (dd < deltax) dd = deltax;
                if (x < a + dd) x = a + dd;
                if (x > b - dd) x = b - dd;
            }
            p = x2p(x);
            checkPoints(&p, deltax, g0 - 1);
            goto NEXTP;
        }
    }

    if (state < loopStateFound1) {
        // give up searching for polymorphism points
        if (g0 < 2) {
            // not enough altruism points, or egoism points if direction
            x = (a + b) / 2.f;
            p = x2p(x);
            checkPoints(&p, deltax, g0 - 1);
            goto NEXTP;
        }
        // linear regression of altruism points, or egoism points if direction
        i = g0;  if (i > 4) i = 4;
        j = regression(g0, g0 + i - 1, 1, &intercept, &slope);
        if (j == 0 || slope == 0) {
            // regression error
            x = (a + b) / 2;
            p = x2p(x);
            checkPoints(&p, deltax, g0 - 1);
        }
        else {
            // extrapolate these points to gn = limit
            x = (d->maximumGenerations - intercept) / slope;
            if (x <= a || x >= b) x = (a + b) / 2;
            dd = (b - a) / 4;
            if (dd < deltax) dd = deltax;
            p = x2p(x);
            checkPoints(&p, dd, g0 - 1);
        }
        goto NEXTP;
    }

    // state = loopStateFound1. search for limit2 between polymorphism and egoism (or altruism if direction)
    if (!(allOutcomes & (1 << resultPolymorphism))) {
        // no polymorphism found. finished
        if (limitsFound & 1) {
            limit2 = limit1;  limitsFound |= 3;
        }
        goto FINISHED;
    }

    if (!(allOutcomes & (1 << (direction ? resultAltruism : resultEgoism)))) {
        // no egoism (altruism) found. finished
        limit2 = endValue;
        goto FINISHED;
    }

    if (ord2) {
        // results out of order. finished
        limit2 = x2p((p2x(resultList[ord2-1].par) + p2x(resultList[ord2].par)) / 2.f);
        limitsFound |= 2;
        goto FINISHED;
    }

    h = g0 + g1 + g2;
    a = p2x(resultList[h-1].par);
    b = p2x(resultList[h].par);

    if (b - a <= deltax) {  // found limit2
        limit2 = x2p((a + b) / 2.);
        limitsFound |= 2;
        goto FINISHED;
    }

    if (g1 + g2 < 2) {
        // not enough polymorphism points
        x = (a + b) / 2;  p = x2p(x);
        checkPoints(&p, deltax, h - 1);
        goto NEXTP;
    }

    if (g2 < 4 && g1 + g2 < 8) {
        // linear regression of polymorphism points
        i = g1 + g2;
        if (i < 2) i = 2;
        if (i > 4) i = 4;
        j = regression(h - i, h - 1, 0, &intercept, &slope);
        if (j == 0 || slope == 0) {
            // regression error
            x = (a + b) / 2;
        }
        else {
            if (direction) x = (1.f - intercept) / slope; // extrapolate to gf = 1
            else x = -intercept / slope;                   // extrapolate to gf = 0
            if (x <= a || x >= b) x = (a + b) / 2;
            dd = (b - a) / 6.f;
            if (dd < deltax) dd = deltax;
            if (x < a + dd) x = a + dd;
            if (x > b - dd) x = b - dd;
        }
        p = x2p(x);
        checkPoints(&p, deltax, h - 1);
        goto NEXTP;
    }

    // extrapolation from polymorphism points not successful
    if (g3 < 2) {
        // not enough egoism points (or altruism if direction)
        x = (a + b + b) / 3;
        p = x2p(x);
        checkPoints(&p, deltax, h - 1);
        goto NEXTP;
    }

    // linear regression of egoism points, or altruism points if direction
    i = g3;
    if (i > 4) i = 4;
    j = regression(h, h + i - 1, 1, &intercept, &slope);
    if (j == 0 || slope == 0) {
        // regression error
        x = (a + b) / 2;
        p = x2p(x);
        checkPoints(&p, deltax, h - 1);
    }
    else {
        // extrapolate to gn = limit
        x = (d->maximumGenerations - intercept) / slope;
        if (x <= a || x >= b) x = (a + b) / 2;
        dd = (b - a) / 4.f;
        if (dd < deltax) dd = deltax;
        p = x2p(x);
        checkPoints(&p, dd, h - 1);
    }
    goto NEXTP;

NEXTP: // Make a new simulation with the calculated guess

    if (p < startValue) {
        errors.reportError("p < start C");
    }

    parameterValue = p;
    setValue(p);
    return;

FINISHED: // Search for limits is finished

    parameterValue = endValue;
    setValue(endValue);
    if (!(limitsFound & 1)) limit1 = startValue;
    if (!(limitsFound & 2)) limit2 = endValue;
    if (!limitsFound) {
        // all simulations give same result
        switch (resultList[0].result) {
        case resultEgoism:
            limit1 = limit2 = startValue;
            break;
        case resultPolymorphism: 
            limit1 = startValue;
            limit2 = endValue;
            break;
        case resultAltruism:
            limit1 = limit2 = endValue;
            break;        
        }
    }
    lastParIndex = parIndex;
    lastFound = limitsFound;
    allOutcomes = 0;
    state = loopStateFinished;

    // put search results in limitList
    if (worker->limitListNum < limitListLength) {
        worker->limitList[worker->limitListNum].limit1 = limit1;
        worker->limitList[worker->limitListNum].limit2 = limit2;
        worker->limitList[worker->limitListNum].result = (int)direction | (limitsFound ? 0 : 8);
        // get y value from next higher loop
        worker->limitList[worker->limitListNum].sValue = worker->loops[1].parameterValue;
        if (ord1 | ord2) worker->limitList[worker->limitListNum].result |= 2;
        if (died) worker->limitList[worker->limitListNum].result |= 4;
        worker->limitListNum++;
    }
    checkLimitList();
}

void ParameterLoop::checkLimitList(void) {
    // check for entries in limitlist with unknown direction
    ParameterLimits * limitList = worker->limitList;
    int limitListNum = worker->limitListNum;
    if (limitListNum < 2) return;
    // find dominating direction
    int countDirection0 = 0, countDirection1 = 0;
    for (int i = 0; i < limitListNum; i++) {
        if (limitList[i].result & 8) continue;   // direction not known
        if (limitList[i].result & 1) countDirection1++; else countDirection0++;
    }
    if (countDirection1 > countDirection0) {      // majority is direction 1
        for (int i = 0; i < limitListNum; i++) {
            if (limitList[i].result & 8) {
                // fix record with unknown direction
                limitList[i].result |= 1;        // mark as direction 1
                limitList[i].result &= ~8;       // avoid correcting twice
                if (limitList[i].limit1 == limitList[i].limit2) {
                    if (limitList[i].limit1 == startValue) {
                        limitList[i].limit1 = limitList[i].limit2 = endValue;
                    }
                    else if (limitList[i].limit1 == endValue) {
                        limitList[i].limit1 = limitList[i].limit2 = startValue;
                    }                 
                }
            }
        }
    }
}
                                  