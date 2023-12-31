/******************************  run.cpp   ************************************
* Author:        Agner Fog
* Date created:  2023-03-27
* Last modified: 2023-12-31
* Version:       3.00.00
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the overhead of running a model, including generation
* loop, parameter loops, and run statistics.
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"


void Altruist::run() {
    // start button pressed
    int state = d.sweepsUsed ? d.sweepState : d.runState;
    if (state == state_pause || state == state_pause_next) {    
        resume();
        return;
    }
    worker->loopsPaused = false;

    if (state == state_idle || state == state_stop || d.parametersChanged) {
        d.currentModel = (models.getModel(d.iModel));      // get model descriptor
        // calculate necessary memory
        d.demeStructureSize = d.modelDemeStructureSize;    // redundant
        if (d.currentModel->demeFields[0].type != 1) {
            errorMessage("Model code fails to define deme size"); return;
        } 
        else {
            d.demeStructureSize = d.currentModel->demeFields[0].offset;
        }
        // calculate required memory for maxIslands. Add 1 for migrant pool or unused record 0
        uint64_t requiredMemory = (d.maxIslands + 1) * d.demeStructureSize;
        if (d.bufferSize < requiredMemory) {
            // must allocate memory
            if (d.bufferSize > 0) delete[] d.demeData;     // free previously allocated memory
            d.demeData = new int8_t[requiredMemory]();     // () means initialize to zero (C++ 11)
            if (d.demeData == 0) {
                errorMessage("Memory allocation failed"); return;
            }
            else {            
                d.bufferSize = requiredMemory;
            }
        }
    }
    // other initializations
    d.rowLength = int32_t(sqrt(d.maxIslands));
    if (d.sweepsUsed) d.sweepState = state_start;
    if (d.runState == state_idle || d.runState == state_stop || d.runState == state_start) {
        // seed random number generator
        d.runState = state_start;
        d.generations = 0;
        if (d.sweepsUsed) d.sweepState = state_start;
        emit doWorkSignal(state_start);
    }
    if (d.sweepType[1] > 0 && d.sweepType[1] <= 2 && d.sweepType[0] >= 3) {
        // make x-y map of limiting parameter values
        // This ParameterLoop sweeps the y values, while an 
        // underlying search finds the limiting x values
        d.graphicsType = graphicsLimits;
    }
    else {
        d.graphicsType = d.graphicsTypeForModel;
    }
    d.parametersChanged = false;
}

void Altruist::pause() {
    // pause run
    int volatile &state = d.sweepsUsed ? d.sweepState : d.runState;
    if (state == state_run || state == state_start) {
        state = state_pause_next;
        setStateTitle();
        emit doWorkSignal(state_pause);
        setStateTitle();
    }
}

void Altruist::resume() {
    // resume after pause
    int volatile &state = d.sweepsUsed ? d.sweepState : d.runState;
    if (state >= state_run && state <= state_pause) {
        state = state_run;
        setStateTitle();
        checkIfParametersChanged();
        emit doWorkSignal(state_run);
    }
}

void Altruist::singleStep() {
    if (d.sweepsUsed == 0) {  
        // single generation step
        if (d.runState == state_stop) {
            return;                              // stop single stepping until resat
        }
        if (d.runState == state_idle) {
            run();
            d.runState = state_single_step;
            checkIfParametersChanged();
            emit doWorkSignal(state_single_step);
            setStateTitle();
        }
        else if (d.runState >= state_run && d.runState <= state_pause) {
            d.runState = state_single_step;
            checkIfParametersChanged();
            emit doWorkSignal(state_single_step);
            setStateTitle();
        }
    }
    else {
        // single simulation in parameter sweep
        if (d.sweepState == state_stop) {
            d.sweepState = state_idle;
            return;                              // stop single stepping until resat
        }
        if (d.sweepState == state_idle) {
            run();
            d.sweepState = state_single_step;
            checkIfParametersChanged();
            emit doWorkSignal(state_single_step);
            setStateTitle();
        }
        else if (d.sweepState >= state_run && d.sweepState <= state_pause) {
            d.sweepState = state_single_step;
            checkIfParametersChanged();
            emit doWorkSignal(state_single_step);
            setStateTitle();
        }    
    }
}

void Altruist::stop() {
    // stop running
    int volatile &state = d.sweepsUsed ? d.sweepState : d.runState;
    if (state == state_stop) {
        d.runState = state_idle;
        d.sweepState = state_idle;
        d.requestUpdate = true;
        setStateTitle();
    }
    else if (state > state_idle) {
        d.runState = state_stop;
        d.sweepState = state_stop;
        d.requestUpdate = true;
        setStateTitle();
        emit doWorkSignal(state_stop);
    }
}

void Altruist::checkIfParametersChanged() {
    // check if parameters changed while running
    int volatile &state = d.sweepsUsed ? d.sweepState : d.runState;
    if (d.parametersChanged && state > state_idle && state < state_stop) {
        // parameters have changed while running
        // check memory requirement
        uint64_t requiredMemory = (d.maxIslands + 1) * d.demeStructureSize;
        if (requiredMemory > d.bufferSize) {
            // must allocate more memory
            int8_t * newBuffer = new int8_t[requiredMemory]();
            if (newBuffer == 0) {
                errorMessage("Memory allocation failed"); return;
                d.maxIslands = d.bufferSize / d.demeStructureSize;
                return;
            }
            int8_t * oldBuffer = d.demeData;
            memcpy (newBuffer, oldBuffer, d.bufferSize);   // copy old data to new buffer
            delete[] oldBuffer;                            // delete old buffer
            d.bufferSize = requiredMemory;                 // set new buffer size
        }
        (*d.currentModel->initFunction)(&d, model_parameters_changed); // model init function parameters changed 
        d.requestUpdate = true;
    }
}

void Altruist::resultReadySlot(int state) {
    // signal received from worker
    if (d.sweepsUsed == 0) { // a simulation is finished
        setStateTitle();
    }
    else { // a simulation in a sweep is finished 
        setStateTitle();
    }
}

void Altruist::setStateTitle() {
    // update the main window title to indicate run state
    int state;
    static int oldState = -1;
    const char * text = "?";
    if (d.sweepsUsed == 0) {  // single simulation
        state = d.runState;   // save d.runState, because it may change during writing
        if (state != oldState) {
            switch (d.runState) {
            case state_idle:
                text = "Altruist - idle"; break;
            case state_start:
                text = "Altruist - start"; break;
            case state_run:
                text = "Altruist - running";
                if (d.sweepsUsed) text = "parameter loop running";
                break;
            case state_single_step:
                text = "Altruist - single step"; break;
            case state_pause_next:
                text = "Altruist - pause next"; break;
            case state_pause:
                text = "Altruist - pause"; break;
            case state_stop:
                text = "Altruist - stop"; break;
            default:
                text = "Altruist - error"; break;
            }
            this->setWindowTitle(text);
        }
        oldState = state;
    }
    else {  // parameter sweep
        state = d.sweepState;  // save d.sweepState, because it may change during writing
        if (state != oldState) {
            switch (d.sweepState) {
            case state_idle:
                text = "Altruist - idle"; break;
            case state_start:
                text = "Altruist - start"; break;
            case state_run:
                text = "Altruist - parameter loop running";
                break;
            case state_single_step:
                text = "Altruist - single step"; break;
            case state_pause_next:
                text = "Altruist - pause next"; break;
            case state_pause:
                text = "Altruist - parameter loop paused"; break;
            case state_stop:
                text = "Altruist - stop"; break;
            default:
                text = "Altruist - error"; break;
            }
            this->setWindowTitle(text);
        }    
        oldState = state;
    }
}

void Worker::doWorkSlot(int state) {
    // time-consuming simulation activated here
    if (d->sweepsUsed) {
        doMultipleSimulations();                 // start or resume parameter sweeps or search
        return;
    }

    if (state == state_start) {
        waitForScreen();
        randomGenerator.randomInit(d->seed);     // seed random number generator
        timer.start();
        d->generations = 0;
        bool singleStep = d->runState == state_single_step;
        (*(models.getModel(d->iModel)->generationFunction))(d, state_start);
        if (singleStep) d->runState = state_single_step;
        else d->runState = state_run;
        fileOutStart();
    }

    if (d->runState == state_pause_next) {
        // call generation function to indicate state
        (*(models.getModel(d->iModel)->generationFunction))(d, d->runState);
        d->runState = state_pause;
        d->requestUpdate = true;
    }
    else if (d->runState == state_single_step) {
        // call generation function
        (*(models.getModel(d->iModel)->generationFunction))(d, state_run);
        fileOutGeneration();
        d->runState = state_pause;
        d->requestUpdate = true;
    }
    else if (d->runState == state_run) {
        // running: call model generation function repeatedly
        do {
            waitForScreen();
            // call generation function
            (*(models.getModel(d->iModel)->generationFunction))(d, d->runState);
            if (d->runState == state_pause_next) {
                d->runState = state_pause;
            }
            d->requestUpdate = true;
            fileOutGeneration();
            if (d->delayms && d->runState == state_run) {
                thread()->msleep(d->delayms);    // delay after screen update
            }
        } while (d->runState == state_run);
    }

    if (d->runState == state_stop) {
        d->timeUsed = timer.elapsed();           // measure elapsed time in milliseconds
        fileOutSimulationFinished();             // write results
        d->runState = state_idle;                // idle
    }
    emit resultReadySignal(d->runState);
}

void Worker::doWholeSimulation() {
    // do a simulation run as part of a parameter sweep or search
    // initialize
    d->runState = state_start;
    randomGenerator.randomInit(d->seed);         // restart and seed random number generator
    timer.start();                               // measure time consumption of this simulation
    d->generations = 0;
    (*(models.getModel(d->iModel)->generationFunction))(d, state_start);
    d->runState = state_run;
    // run: call model generation function repeatedly
    do {
        // call generation function
        (*(models.getModel(d->iModel)->generationFunction))(d, d->runState);
        if (d->runState == state_pause_next) {
            d->runState = state_pause;
        }

        // save simulation result and the corresponding parameter values
        lastResult = d->simulationResult;
        for (int iLoop = 0; iLoop < maxSweeps; iLoop++) {
            lastParameter[iLoop] = loops[iLoop].lastValue = loops[iLoop].parameterValue;
        }
        // file output for generation not used in parameter loops
        // if (d->runState != state_stop) fileOutGeneration();

    } while (d->runState == state_run);

    // save simulation result and the corresponding parameter values
    lastResult = d->simulationResult;
    for (int iLoop = 0; iLoop < maxSweeps; iLoop++) {
        lastParameter[iLoop] = loops[iLoop].parameterValue;
    }

    if (d->runState == state_stop) {
        d->timeUsed = timer.elapsed();           // measure elapsed time in milliseconds
        // file output for simulation is redundant in parameter loops because results are written by fileOutSweepNext
        //fileOutSimulationFinished();
        d->runState = state_idle;                // idle
    }

    int32_t delay1 = d->delayms;
    if (delay1 > 0 && d->sweepState == state_run) {
        if (delay1 < 200) delay1 = 200;          // minimum delay 
        thread()->msleep(delay1);                // delay after screen update or status line update
    }

}

void Worker::waitForScreen() {
    // wait for curent graphic screen to be updated before changing islands etc.
    int n = 0;
    while (d->requestUpdate) {
        thread()->msleep(1);                     // wait for graphics screen update before changing
        if (++n > 20) break;                     // wait max 20 ms
    }
}

void statisticsInit0(AltruData * d) {
    // reset statistics counters before simulation start
    d->nSum = 0;                     // number of items in sums in fgsum
    d->nSteady = 0;                  // number of generations before steady state statistics
    d->sumMutations = 0;             // sum number of mutations
    d->sumExtinctions = 0;           // sum number of extinctions
    d->sumMigrants = 0;              // sum number of migrants
    statisticsInit1(d);
}

void statisticsInit1(AltruData * d) {
    // reset statistics counters before each generation
    for (int i = 0; i < d->nLoci; i++) {
        d->mutations[i][0] = d->mutations[i][1] = 0;    
        d->totalPhenotypes[i] = 0;
        d->genePool[i][0] = d->genePool[i][1] = 0;
    }
    d->demesDied = 0;
    d->migrantsTot = 0;
    d->totalPopulation = 0;
    d->inhabitedDemes = 0;
    d->altruistDemes = 0;
    //d->egoistDemes = 0;
    d->stopCause = stopCauseNone;
    d->sumExtinctions = 0;
}

void statisticsUpdate(AltruData * d) {
    // update statistics after each generation
    if (d->totalPopulation > 0) {
        d->fAltru = float(d->totalPhenotypes[0]) / float(d->totalPopulation);
    }
    else {
        d->fAltru = 0.f;
    }
    for (int i = 0; i < maxLoci; i++) {
        if (d->locusUsed[i] && d->totalPopulation > 0) {
            d->geneFraction[i] = float(d->genePool[i][1]) / float(d->totalPopulation * 2);
            d->sumMutations += d->mutations[i][0] + d->mutations[i][1];
        }
        else {
            d->geneFraction[i] = 0.f;
        }
    }
    d->sumExtinctions += d->demesDied;
    d->sumMigrants += d->migrantsTot;
}

void checkStopCriterion(AltruData * d) {
    // check if stop criterion reached

    int i;  // loop counter
    //int stopCause; // 1: aborted, 2: all died, 3: max number of generations reached, 4: stop criterion reached
    if (d->runState == state_stop && d->stopCause == stopCauseNone) {
        d->stopCause = stopCauseAbort;           // simulation stopped manually
        d->simulationResult = resultAborted;
    }
    if (d->totalPopulation == 0) {
        d->stopCause = stopCauseDied;            // all died
        d->simulationResult = resultDied;
    }
    if (d->generations >= d->maximumGenerations) {
        d->stopCause = stopCauseMaxGen;          // max number of generations reached
        d->simulationResult = resultPolymorphism;
    }
    if (d->generations >= d->minimumGenerations && d->totalPopulation > 0) { // stop criterion used only after minimum generations
        bool crit = true;                        // and-combination of criterion for all loci
        float mutantFraction = 0.f;              // fraction of mutant genes or phenotypes
        int locus = 0;
        if (d->stopCause == stopCauseNone) {     // check stop criterion
            switch (d->stopCriterion) {
            case stopCritNone:                   // no criterion
                break;
            case stopIfGenesUniform:             // fraction of all genes uniform
                for (i = 0; i < d->nLoci; i++) {
                    if (d->locusUsed[i]) {
                        mutantFraction = float(d->genePool[i][1]) / float(d->totalPopulation * 2);
                        if (mutantFraction < d->stopCriterionDegree && mutantFraction > 1.f - d->stopCriterionDegree) {
                            crit = false;
                        }
                    }
                }
                if (crit) {
                    mutantFraction = float(d->genePool[0][1]) / float(d->totalPopulation * 2);
                    d->simulationResult = mutantFraction > 0.5f ? resultAltruism : resultEgoism;
                    d->stopCause = stopCauseCriterion;  // stop criterion reached 
                }
                break;
            case stopIfPhenoUniform:  // fraction of all phenotypes >= StopCriterionDegree
                for (i = 0; i < d->nLoci; i++) {
                    if (d->locusUsed[i]) {
                        mutantFraction = d->totalPhenotypes[i] / d->totalPopulation;
                        if (mutantFraction < d->stopCriterionDegree && mutantFraction > 1.f - d->stopCriterionDegree) {
                            crit = false;
                        }
                    }
                }
                if (crit) {
                    mutantFraction = d->totalPhenotypes[0] / d->totalPopulation;
                    d->simulationResult = mutantFraction > 0.5f ? resultAltruism : resultEgoism;
                    d->stopCause = stopCauseCriterion;     // stop criterion reached 
                }
                break;
            default:  // any locus
                locus = (d->stopCriterion - stopIfGene0Uniform) / 2;  // find which locus criterion applies to
                if (locus > d->nLoci) break;               // locus not used
                if (d->genePool[locus][1] >= d->totalPopulation * 2 * d->stopCriterionDegree) {
                    d->stopCause = stopCauseCriterion;     // stop criterion reached: mutant fixated
                    d->simulationResult = resultAltruism;
                    break;
                }
                if ((d->stopCriterion - stopIfGene0Mutant & 1)
                && d->genePool[locus][1] <= d->totalPopulation * 2 * (1.f - d->stopCriterionDegree)) {
                    d->stopCause = stopCauseCriterion;     // stop criterion reached: wild type fixated
                    d->simulationResult = resultEgoism;
                }
                break;
            }
        }
    }
    if (d->stopCause != stopCauseNone) {
        // stop criterion reached or stopped for any reason
        d->runState = state_stop;   // stop simulation
        // calculate statistics
        if (d->totalPopulation > 0) {
            for (i = 0; i < d->nLoci; i++) {
                d->geneFraction[i] = float(d->genePool[i][1]) / float(d->totalPopulation * 2);
            }
            d->fAltru = float(d->totalPhenotypes[0]) / float(d->totalPopulation);
        }
        else {  // totalPopulation = 0
            for (i = 0; i < d->nLoci; i++) {
                d->geneFraction[i] = 0.f;
            }
            d->fAltru = 0.f;
        }
    }
}




void addNeighbor(AltruData * d, int32_t deme, int32_t neighbors[], int & n) {
    // used by function findNeighbors to add a neighbor island to a list, unless it is empty
    if (*(int32_t*)(d->demeData + deme * d->modelDemeStructureSize + d->modelPopOffset) != 0) { // check if not empty
        neighbors[n++] = deme;
    }
}

int findNeighbors(AltruData * d, int32_t deme, int32_t neighbors[]) {
    // find all (inhabited) neighbors to a deme, according to migration topology
    if (d->modelDemeStructureSize == 0) return 0;          // escape if not initialized
    switch (d->maxIslands) {        
    case 1: 
        return 0;                                          // there are no neighbors
    case 2:
        neighbors[0] = 1-deme;
        return 1;                                          // there is only one neighbor
    default:
        break;                                             // there is more than 1 neighbor
    }
    if (d->migrationTopology == topology3neighb) {         // each island has three neighbors
        d->rowLength = d->maxIslands / 2;                  // only two rows 
    }

    int32_t row = uint32_t(deme) / uint32_t(d->rowLength); // row number
    int32_t coloumn = deme - row * d->rowLength;           // column number
    int32_t lastRow = uint32_t(d->maxIslands - 1) / uint32_t(d->rowLength) - 1; // last row

    int n = 0;  // count number of neighbors
    int32_t j = 0;  // temporary index

    switch (d->migrationTopology) {

    case topologyLinear:  // linear organization
        if (deme > 0) addNeighbor(d, deme - 1, neighbors, n);                  // left neighbor
        if (deme < d->maxIslands-1) addNeighbor(d, deme + 1, neighbors, n);    // right neighbor
        break;

    case topology3neighb:  // two rows. 3 neighbors
        // continue in next case
    case topologyQuadratic:  // quadratic organization
        if (coloumn > 0) addNeighbor(d, deme - 1, neighbors, n);               // west neighbor
        if (row > 0) addNeighbor(d, deme - d->rowLength, neighbors, n);        // north neighbor
        if (row < lastRow) {
            addNeighbor(d, deme + d->rowLength, neighbors, n);                 // south neighbor
            if (coloumn < d->rowLength-1) addNeighbor(d, deme + 1, neighbors, n); // east neighbor
        }
        else {  // last row. no south neighbor
            if (deme < d->maxIslands-1) addNeighbor(d, deme + 1, neighbors, n);// east neighbor        
        }
        break;

    case topologyHoneycomb:  // honeycomb organization
        if (coloumn > 0) addNeighbor(d, deme - 1, neighbors, n);               // west neighbor
        if (row > 0) {
            addNeighbor(d, deme - d->rowLength, neighbors, n);                 // north neighbor
            if (row & 1) {  // odd row
                if (coloumn < d->rowLength-1) addNeighbor(d, deme - d->rowLength + 1, neighbors, n);  // north east neighbor
            }
            else {  // even row
                if (coloumn > 0) addNeighbor(d, deme - d->rowLength - 1, neighbors, n);  // north west neighbor
            }
        }
        if (row < lastRow) {      // not last row
            addNeighbor(d, deme + d->rowLength, neighbors, n);                 // south neighbor
            if (coloumn < d->rowLength - 1) {                                  // not last column
                addNeighbor(d, deme + 1, neighbors, n);                        // east neighbor
            }
            if (row & 1) {  // odd row
                if (coloumn < d->rowLength - 1) {
                    addNeighbor(d, deme + d->rowLength + 1, neighbors, n);     // south east neighbor
                }
            }
            else {  // even row
                if (coloumn > 0) {
                    addNeighbor(d, deme + d->rowLength - 1, neighbors, n);     // south west neighbor
                }
            }
        }
        else {  // last row. no south neighbor
            if (deme < d->maxIslands - 1) addNeighbor(d, deme + 1, neighbors, n); // east neighbor
        }
        break;

    case topologyOctal:  // octal organization
        if (coloumn > 0)  addNeighbor(d, deme - 1, neighbors, n);              // west neighbor
        if (row > 0) {
            addNeighbor(d, deme - d->rowLength, neighbors, n);                 // north neighbor
            if (coloumn < d->rowLength-1) addNeighbor(d, deme - d->rowLength + 1, neighbors, n); // north east neighbor
            if (coloumn > 0) addNeighbor(d, deme - d->rowLength - 1, neighbors, n); // north west neighbor            
        }
        if (row < lastRow) {      // not last row
            addNeighbor(d, deme + d->rowLength, neighbors, n);                 // south neighbor
            if (coloumn < d->rowLength - 1) {                                  // not last column
                addNeighbor(d, deme + 1, neighbors, n);                        // east neighbor
                addNeighbor(d, deme + d->rowLength + 1, neighbors, n);         // south east neighbor
            }
            if (coloumn > 0) {
                addNeighbor(d, deme + d->rowLength - 1, neighbors, n);         // south west neighbor
            }
        }
        else {  // last row. no south neighbor
            if (deme < d->maxIslands - 1) addNeighbor(d, deme + 1, neighbors, n); // east neighbor
        }
        break;

    case topologyRandom:                         // random migration. no organization
        j = d->ran->iRandom(0, d->maxIslands-2); // pick one of the neighbors
        if (j >= deme) j++;                      // skip own island
        neighbors[n++] = j;
        break;

    case topologyCommonPool:                     // one common migrant pool
        neighbors[n++] = d->maxIslands;          // last record is migrant pool
        break;

    case topologyFloating:                       // floating territories
        // implemented in terrain.cpp
        break;    
    }

    return n; // return number of islands in neighbors list
}


/***********************************************************************
  Stochastic functions applied to simulation
***********************************************************************/

void combineGenes(int32_t genePool[2], int32_t genotypes[3], RandomVariates * ran) {
    // stochastic combination of a gene pool into genotypes for a biallelic locus
    // genePool = total number of wild type and mutant genes at locus
    // ran = pointer to random variates generator
    // result:
    // genotypes[0] = number of (n,n) genotypes
    // genotypes[1] = number of (n,a) genotypes
    // genotypes[2] = number of (a,a) genotypes

    int32_t gn = genePool[0];                              // wild type genes
    int32_t ga = genePool[1];                              // mutant genes
    int32_t m  = uint32_t(gn + ga) / 2u;                   // total number of individuals
    if (gn > 2 * m) gn = 2 * m;                            // this can occur if gn is odd and ga is 0
    int32_t mn_ = ran->hypergeometric(m, gn, 2 * m);       // individuals with 'n' on first allele
    int32_t mnn = ran->hypergeometric(mn_, gn - mn_, m);   // individuals with 'n' on both alleles
    int32_t mna = gn - 2 * mnn;                            // heterozygote individuals
    int32_t maa = mnn + m - gn;                            // individuals with 'a' on both alleles
    // mnn + mna + maa = m
    genotypes[0] = mnn;
    genotypes[1] = mna;
    genotypes[2] = maa;
}

void differentialGrowth(int32_t genePool[2], double growthRate[3], int32_t genotypes[3], RandomVariates * ran) {
    // This simulates growth of a gene pool with fecundity selection where fitness or growth 
    // rate for a biallelic locus is specified for each genotype. 
    // Number of offspring follows Poisson distribution
    // genePool = total number of wild type (n) and mutant genes (a) at locus
    // growthRate is the average number of offspring of (n,n), (n,a) and (a,a) genotypes, respectively
    // ran = pointer to random variates generator
    // genotypes outputs the number of (n,n), (n,a) and (a,a) genotypes in the offspring generation
    int32_t gn = genePool[0];                              // wild type genes
    int32_t ga = genePool[1];                              // mutant genes
    int32_t m = uint32_t(gn + ga) / 2u;                    // total number of individuals
    int32_t mn_ = ran->hypergeometric(m, gn, 2 * m);       // individuals with 'n' on first allele
    int32_t mnn = ran->hypergeometric(mn_, gn - mn_, m);   // individuals with 'n' on both alleles
    int32_t mna = gn - 2 * mnn;                            // heterozygote individuals
    int32_t maa = mnn + m - gn;                            // individuals with 'a' on both alleles
    // growth of each genotype
    // (The sum of Poisson variates is a Poisson variate:)
    genotypes[0] = ran->poisson(mnn * growthRate[0]);
    genotypes[1] = ran->poisson(mna * growthRate[1]);
    genotypes[2] = ran->poisson(maa * growthRate[2]);
}

void growthAndSelection(int32_t genes[2], int32_t numChildren, int dominance, float relativeFitness, RandomVariates * ran) {
    // This simulates growth and selection at a secondary locus when the number of offspring is 
    // determined by other factors
    // genes is the number of neutral and mutant genes at a given locus. Both input and output
    // numChildren is the number of children determined by other factors
    // dominance is the dominance of the mutant gene
    // relativeFitness is the fitness of the mutant phenotype relative to the neutral wild type
    // ran is a pointer to the random variates generator
    int32_t genotypes[3];                                  // each genotype at the given locus
    double  distribution[3];                               // probability distribution for offspring of each genotype
    if (genes[0] + genes[1] == 0) return;                  // no parents

    if (relativeFitness == 1.f || numChildren == 0) {      // no selection at locus, only growth and drift
        genes[1] = ran->binomial(numChildren*2, double(genes[1]) / double(genes[0] + genes[1]));
        genes[0] = numChildren*2 - genes[1];
    }
    else {                                                 // growth and selection at locus
        combineGenes(genes, genotypes, ran);  // split into genotypes
        for (int i = 0; i < 3; i++) distribution[i] = genotypes[i];
        distribution[2] *= relativeFitness;                // fitnes ratio
        switch (dominance) {
        case recessive:
            break;
        case dominant:
            distribution[1] *= relativeFitness;
            break;
        case incompleteDominant:                           // half dominant
            distribution[1] *= sqrtf(relativeFitness);
        }
        ran->multinomial(genotypes, distribution, numChildren, 3);  // growth, selection, and drift
        genes[0] = genotypes[0] * 2 + genotypes[1];
        genes[1] = numChildren*2 - genes[0];
    }
}
