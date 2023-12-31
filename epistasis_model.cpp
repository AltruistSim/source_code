/***************************  epistasis_model.cpp   ***************************
* Author:        Agner Fog
* Date created:  1995-08-13
* Last modified: 2023-12-31
* Version:       3.00.00
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines a model of epistasis in a structured or viscous population.
* Epistasis can lead to evolution through punctuated equilibria
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"


void epistasisInitFunction(AltruData * d, int index);
void epistasisGenerationFunction(AltruData * d, int mode);


/******************************************************************************
*         structures and definitions
******************************************************************************/

// List of model-specific parameters
// These parameters will appear in parameter list and dialog box.
// Maximum number of model-specific parameters is maxModelSpecificParameters defined in altruist.h
// ParameterDef::type can be: 1: integer, 2: float, 3: boolean
// ParameterDef::num must be 1 if not an array
// Offsets and names must be unique
static const ParameterDef epistasisParameterDefinitions[] {
    {0, 1, 0, "Model-specific parameters"},
    {0, 0, 0, 0}
};

// Deme structure, describing each deme or territory
struct EpistasisDeme {
    int32_t nn;              // number of individuals in deme
    int32_t nmax;            // carrying capacity or max deme size
    int32_t gg[3][2];        // gene pool, locus A, B, C
    int32_t mutationAge[3];  // age of each mutant gene before combined with mutants at other loci
};

const int locusA = 0;        // A locus
const int locusB = 1;        // B locus
const int locusC = 2;        // C locus


/*
struct DemeFieldDescriptor {
    int32_t type;            // 0: end of list, 
                             // 1: size of Deme structure, 
                             // 2: carrying capacity (max individuals)
                             // 3: population (number of individuals)
                             // 4: gene count, 
                             // 5: genotype or phenotype count
                             // 8: group property
    int32_t varType;         // varInt16 or varInt32 or varFloat
    int32_t offset;          // offset into Deme structure
    int32_t statistic;       // 1: calculate sum, 2: calculate mean
                             // 4: possible stop criterion, stop when above certain value
                             // 8: possible stop criterion, stop when below certain value
    int32_t graphics;        // show in graphics display:
                             // 2:  max size
                             // 3:  population. divide by max size or nMaxPerDeme to get relative size
                             // 4:  area. divide by territorySizeMax to get relative size
                             // 10: gene count for mutant, primary locus
                             // 11: gene count for mutant, secondary locus
                             // 12: gene count for mutant, third locus
                             // 20: group property
    const char * name;       // name of variable
};
*/

// List of fields in Deme structure
static const DemeFieldDescriptor epistasisDemeDescriptors[] = {
    {1, varInt32, sizeof(EpistasisDeme), 0, 0, 0},
    {3, varInt32, demeFieldOffset(EpistasisDeme,nn),             0, 3,  "population"},
    {4, varInt32, demeFieldOffset(EpistasisDeme,gg[locusA][0]),  8, 0,  "nonA gene"},
    {4, varInt32, demeFieldOffset(EpistasisDeme,gg[locusA][1]),  5, 10, "A gene"},
    {4, varInt32, demeFieldOffset(EpistasisDeme,gg[locusB][0]),  2, 0,  "nonB gene"},
    {4, varInt32, demeFieldOffset(EpistasisDeme,gg[locusB][1]),  2, 11, "B gene"},
    {4, varInt32, demeFieldOffset(EpistasisDeme,gg[locusC][0]),  2, 0,  "nonC gene"},
    {4, varInt32, demeFieldOffset(EpistasisDeme,gg[locusC][1]),  2, 12, "C gene"},
    {0, 0, 0, 0, 0, 0}   // mark end of list
};


static const ModelDescriptor epistasisModel = {
    "Epistasis model",
    "Model where individual fitness depends on two or three mutually interacting loci "
    "in a structured or viscous population. "
    "This is useful for simulating Wright's shifting balance theory. "
    "Evolution through punctuated equilibria may happen when two genes have to be "
    "combined before a gain in fitness occurs, while each gene alone is reducing the fitness.",
    epistasisParameterDefinitions,
    epistasisDemeDescriptors,
    &epistasisInitFunction,
    &epistasisGenerationFunction
};

// This will register the epistasis model in the global model descriptor list
Construct epistasisConstructor(epistasisModel);

const char * userDataNames[] = {
    "A gene age",
    "B gene age",
    "C gene age",
    "total wild type",
    "total A type",
    "total B type",
    "total C type",
    "total AB type",
    "total ABC type",
    "total other types"
};


/******************************************************************************
*         initialization function
* state values:
* 1: model changed
* 2: parameters changed
******************************************************************************/

void epistasisInitFunction(AltruData * d, int state) {
    switch (state) {
    case model_initialize:  // initialization
        // model version
        d->modelVersionMajor = 3;
        d->modelVersionMinor = 0;
        d->modelDemeStructureSize = sizeof(EpistasisDeme);
        d->modelPopOffset = demeFieldOffset(EpistasisDeme,nn);

        // define names
        d->sLocusName[locusA] = "A";
        d->sLocusName[locusB] = "B";
        d->sLocusName[locusC] = "C";
        d->sGeneName[locusA][0] = "nonA";
        d->sGeneName[locusA][1] = "A";
        d->sGeneName[locusB][0] = "nonB";
        d->sGeneName[locusB][1] = "B";
        d->sGeneName[locusC][0] = "nonC";
        d->sGeneName[locusC][1] = "C";
        d->sPhenotypeName[locusA][0] = "nonA";
        d->sPhenotypeName[locusA][1] = "A";
        d->sPhenotypeName[locusB][0] = "nonB";
        d->sPhenotypeName[locusB][1] = "B";
        d->sPhenotypeName[locusC][0] = "nonC";
        d->sPhenotypeName[locusC][1] = "C";

        // enable and disable fields in parameters dialog boxes
        d->bIndividualPropertiesUsed = 0xFF;
        d->bGeographyParametersUsed = 0x241106;
        d->bGroupPropertiesUsed = 0;
        d->bSelModels = 1;
        d->bFitFunc = 0;
        d->bExtinctionPatterns = 0;
        d->bWarPatterns = 0;
        d->nLoci = 3; // number of loci
        if (d->locusUsed[locusC]) d->locusUsed[locusB] = true;
        if (d->locusUsed[locusB]) d->locusUsed[locusA] = true;
        d->graphicsTypeForModel = graphicsIslands;
        d->nIslands = d->maxIslands;
        d->selectionModel = selectionFecundity;
        d->bStopCriterionUsed = 1;
        for (int locus = 0; locus < d->nLoci; locus++) {
            if (d->locusUsed[locus]) {
                d->bStopCriterionUsed |= 3 << (locus*2 + 3);
            }
        }
        d->userDataNames = userDataNames;
        d->nUser = 10; // sizeof() doesn't work on userDataNames in Microsoft compiler
        d->bUserDataType = 0xFFF8;
        d->fitRowLocus = locusB;
        d->fitColLocus = locusA;
        d->fitSupColLocus = locusC;
        d->fitConditionName = 0;
        break;
    case model_parameters_changed:     // parameters changed
        d->nLoci = 3;
        d->nIslands = d->maxIslands;
        break;    
    }
}

static void errorCheck(AltruData * d, EpistasisDeme * deme) {
    // check deme data for consistency
    char text[64];
    int ideme = deme - (EpistasisDeme*)(d->demeData);

    if (deme->gg[locusC][0] + deme->gg[locusC][1] != 0 && !d->locusUsed[locusC]) {
        errors.reportError("unused locus != 0");
    }

    if (deme->gg[locusA][0] < 0 || deme->gg[locusA][1] < 0
        || deme->gg[locusA][0] + deme->gg[locusA][1] != deme->nn * 2
        || (deme->gg[locusA][0] + deme->gg[locusA][1] & 1)) {
        sprintf_s(text, "error generation %i, deme %i, nn %i, ggA %i %i",
            (int)d->generations, ideme,
            deme->nn, deme->gg[locusA][0], deme->gg[locusA][1]);
        errors.reportError(text);
    }
    if ((d->locusUsed[locusB] && deme->gg[locusB][0] + deme->gg[locusB][1] != deme->nn * 2) ||
        (d->locusUsed[locusC] && deme->gg[locusC][0] + deme->gg[locusC][1] != deme->nn * 2)) {
        sprintf_s(text, "error nn %i, ggB %i %i, ggC %i %i",
            deme->nn, deme->gg[locusB][0], deme->gg[locusB][1], deme->gg[locusC][0], deme->gg[locusC][1]);
        errors.reportError(text);
    }
}

// calculate fitness of all 27 possible genotypes from finess of 8 possible phenotypes and 
// dominance of each or the three loci
static void fitnessOfGenotypes(AltruData * d, double fitness[]) {
    // calculate fitness of each genotype from d->fit[] that contains fitness for each phenotype
    int ia, ib, ic;  // indexes for each locus: 0: homozygote wild type, 1: heterozygote, 2: homozygote mutant
    int j;           // index into fitness array

    // first, fill out fitness array for all homozygote types:
    for (ic = 0; ic < 2; ic++) {
        for (ib = 0; ib < 2; ib++) {
            for (ia = 0; ia < 2; ia++) {
                fitness[ia * 2 + ib * 6 + ic * 18] = d->fit[ia + ib * 2 + ic * 4];
            }
        }
        if (!d->locusUsed[locusC]) break;
    }

    // calculate A heterozygotes
    for (ic = 0; ic < 3; ic += 2) {
        for (ib = 0; ib < 3; ib += 2) {
            j = 1 + ib * 3 + ic * 9;   // heterozygote in A
            switch (d->dominance[locusA]) {
            case recessive:            // recessive
                fitness[j] = fitness[j-1];
                break;
            case dominant:             // dominant
                fitness[j] = fitness[j+1];
                break;
            case incompleteDominant:   // half dominant
                fitness[j] = (fitness[j-1] + fitness[j+1]) * 0.5;
                break;
            }
        }
        if (!d->locusUsed[locusC]) break;
    }

    // calculate B heterozygotes
    for (ic = 0; ic < 3; ic += 2) {
        for (ia = 0; ia < 3; ia++) {
            j = ia + 3 + ic * 9;       // heterozygote in B
            switch (d->dominance[locusB]) {
            case recessive:            // recessive
                fitness[j] = fitness[j-3];
                break;
            case dominant:             // dominant
                fitness[j] = fitness[j+3];
                break;
            case incompleteDominant:   // half dominant
                fitness[j] = (fitness[j-3] + fitness[j+3]) * 0.5;
                break;
            }
        }
        if (!d->locusUsed[locusC]) break;
    }

    // calculate C heterozygotes
    if (d->locusUsed[locusC]) {
        for (ib = 0; ib < 3; ib++) {
            for (ia = 0; ia < 3; ia++) {
                j = ia + ib * 3 + 9;   // heterozygote in C
                switch (d->dominance[locusC]) {
                case recessive:        // recessive
                    fitness[j] = fitness[j - 9];
                    break;
                case dominant:         // dominant
                    fitness[j] = fitness[j + 9];
                    break;
                case incompleteDominant:// half dominant
                    fitness[j] = (fitness[j - 9] + fitness[j + 9]) * 0.5;
                    break;
                }
            }
        }
    }
}

// stochastic division of a gene pool of three biallelic loci into all 27 possible genotypes
static void genePool2Genotypes(AltruData * d, EpistasisDeme * deme, int32_t genotypes[27]) {
    int32_t genoA[3];                  // genotypes of each combination at A locus
    int32_t genoB[3];                  // genotypes of each combination at B locus
    int32_t genoC[3] = {0,0,0};        // genotypes of each combination at C locus
    int32_t genoAB[9];                 // genotypes of each combination at A and B locus
    int i;

    for (i = 0; i < 27; i++) genotypes[i] = 0;   // reset all genotypes
    if (deme->gg[locusC][1] == 0) {
        // shortcut for no C mutants
        if (deme->gg[locusB][1] == 0) {
            // shortcut for no B mutants
            if (deme->gg[locusA][1] == 0) {
                // no mutants at all
                genotypes[0] = deme->gg[locusA][0] / 2u;
            }
            else {
                // only A mutants
                combineGenes(deme->gg[locusA], genoA, d->ran);
                for (i = 0; i < 3; i++) genotypes[i] = genoA[i];
            }
        }
        else if (deme->gg[locusA][1] == 0) {
            // only B mutants
            combineGenes(deme->gg[locusB], genoB, d->ran);
            for (i = 0; i < 3; i++) genotypes[i*3] = genoB[i];
        }
        else {
            // A and B mutants, no C mutants
            // get genotypes of each combination at A locus:
            combineGenes(deme->gg[locusA], genoA, d->ran);
            // get genotypes of each combination at B locus:
            combineGenes(deme->gg[locusB], genoB, d->ran);
            // distribute (0,0) at B locus into all combinations at A locus
            d->ran->multiHypergeometric(genotypes, genoA, genoB[0], 3);
            // calculate remaining A genotypes
            for (i = 0; i < 3; i++) genoA[i] -= genotypes[i];
            // distribute (0,B) at B locus into all combinations at A locus
            d->ran->multiHypergeometric(genotypes + 3, genoA, genoB[1], 3);
            // calculate remaining A genotypes
            for (i = 0; i < 3; i++) genotypes[i + 6] = genoA[i] - genotypes[i + 3];
        }
    }
    else {
        // both A, B, and C mutants
        // (We are making no shortcut for C and no A or no B, because this is rare)
        // get genotypes of each combination at A locus:
        combineGenes(deme->gg[locusA], genoA, d->ran);
        // get genotypes of each combination at B locus:
        combineGenes(deme->gg[locusB], genoB, d->ran);
        // get genotypes of each combination at C locus:
        combineGenes(deme->gg[locusC], genoC, d->ran);
        // distribute (0,0) at B locus into all combinations at A locus
        d->ran->multiHypergeometric(genoAB, genoA, genoB[0], 3);
        // calculate remaining A genotypes
        for (i = 0; i < 3; i++) genoA[i] -= genoAB[i];
        // distribute (0,B) at B locus into all combinations at A locus
        d->ran->multiHypergeometric(genoAB + 3, genoA, genoB[1], 3);
        // calculate remaining A genotypes
        for (i = 0; i < 3; i++) genoAB[i + 6] = genoA[i] - genoAB[i + 3];

        // distribute (0,0) at C locus into all combinations at A and B locus
        d->ran->multiHypergeometric(genotypes, genoAB, genoC[0], 9);
        // calculate remaining AB genotypes
        for (i = 0; i < 9; i++) genoAB[i] -= genotypes[i];
        // distribute (0,C) at C locus into all combinations at A and B locus
        d->ran->multiHypergeometric(genotypes + 9, genoAB, genoC[1], 9);
        // get remaining A and B genotypes
        for (i = 0; i < 9; i++) genotypes[i + 18] = genoAB[i] - genotypes[i + 9];
    }

    // count mutant phenotypes over all generations
    // totalPhenotypes is used as follows:
    // totalPhenotypes[0] : no A, no B, no C
    // totalPhenotypes[1] : A, no B, no C
    // totalPhenotypes[2] : B, no A, no C
    // totalPhenotypes[3] : C, no A, no B
    // totalPhenotypes[4] : A and B, no C
    // totalPhenotypes[5] : A and B and C
    // totalPhenotypes[6] : all other combinations
    if (deme->gg[locusC][1] == 0) {
        for (int i = 0; i < 9; i++) {
            genoAB[i] = genotypes[i];  // get genoAB if not given above
        }
    }

    int32_t phenoAB[4];                // phenotypes of each combination at A and B locus
    phenoAB[0] = genoAB[0];
    phenoAB[1] = genoAB[2];
    phenoAB[2] = genoAB[6];
    phenoAB[3] = genoAB[8];
    // partial dominance ignored. Treat as dominant
    if (d->dominance[locusB] == recessive) {
        if (d->dominance[locusA] == recessive) { // both recessive
            phenoAB[0] += genoAB[1] + genoAB[3] + genoAB[4];
            phenoAB[1] += genoAB[5];
            phenoAB[2] += genoAB[7];
        }
        else {  // A dominant, B recessive
            phenoAB[0] += genoAB[3];
            phenoAB[1] += genoAB[1] + genoAB[4] + genoAB[5];
            phenoAB[3] += genoAB[7];
        }
    }
    else {  // B dominant    
        if (d->dominance[locusA] == recessive) { // A recessive, B dominant
            phenoAB[0] += genoAB[1];
            phenoAB[2] += genoAB[3] + genoAB[4] + genoAB[7];
            phenoAB[3] += genoAB[5];
        }
        else {  // A and B dominant
            phenoAB[1] += genoAB[1];
            phenoAB[2] += genoAB[3];
            phenoAB[3] += genoAB[4] + genoAB[5] + genoAB[7];
        }    
    }
    if (d->locusUsed[locusC] && genoC[1] + genoC[2] && deme->nn) {
        // C present
        int c, ac, bc, abc;
        if (d->dominance[locusC] == recessive) {           // C recessive
            c   = phenoAB[0] * genoC[2] / deme->nn;        // C only
            ac  = phenoAB[1] * genoC[2] / deme->nn;        // C and A
            bc  = phenoAB[2] * genoC[2] / deme->nn;        // C and B
            abc = phenoAB[3] * genoC[2] / deme->nn;        // ABC
        }
        else {   // C dominant
            c   = phenoAB[0] * (genoC[1] + genoC[2]) / deme->nn;
            ac  = phenoAB[1] * (genoC[1] + genoC[2]) / deme->nn;
            bc  = phenoAB[2] * (genoC[1] + genoC[2]) / deme->nn;
            abc = phenoAB[3] * (genoC[1] + genoC[2]) / deme->nn;
        }
        d->totalPhenotypes[0] += phenoAB[0] - c;
        d->totalPhenotypes[1] += phenoAB[1] - ac;
        d->totalPhenotypes[2] += phenoAB[2] - bc;
        d->totalPhenotypes[3] += c;
        d->totalPhenotypes[4] += phenoAB[3] - abc;
        d->totalPhenotypes[5] += abc;
        int other = deme->nn - c - ac - bc - abc;
        if (other > 0) d->totalPhenotypes[6] += other;
    }
    else {
        // no C
        d->totalPhenotypes[0] += phenoAB[0];
        d->totalPhenotypes[1] += phenoAB[1];
        d->totalPhenotypes[2] += phenoAB[2];
        d->totalPhenotypes[4] += phenoAB[3];
    }
}

// reduce all 27 possible genotypes of three biallelic loci into a gene pool
static void genotypes2GenePool(AltruData * d, int32_t genotypes[27], EpistasisDeme * deme) {
    for (int i = 0; i < 3; i++) deme->gg[i][0] = deme->gg[i][1] = 0; // reset all
    for (int ia = 0; ia < 3; ia++) {                       // loop through locus A genotypes
        for (int ib = 0; ib < 3; ib++) {                   // loop through locus B genotypes
            for (int ic = 0; ic < 3; ic++) {               // loop through locus C genotypes
                int j = ia + ib * 3 + ic * 9;              // index into genotypes array
                int32_t g0, g1;                            // randomize genes for offspring of heterozygotes
                switch (ia) {
                case 0:
                    deme->gg[0][0] += genotypes[j] * 2;
                    break;
                case 1:
                    g0 = d->ran->binomial(genotypes[j] * 2, 0.5);
                    g1 = genotypes[j] * 2 - g0;
                    deme->gg[0][0] += g0;
                    deme->gg[0][1] += g1;
                    break;
                case 2:
                    deme->gg[0][1] += genotypes[j] * 2;
                    break;
                }
                switch (ib) {
                case 0:
                    deme->gg[1][0] += genotypes[j] * 2;
                    break;
                case 1:
                    g0 = d->ran->binomial(genotypes[j] * 2, 0.5);
                    g1 = genotypes[j] * 2 - g0;
                    deme->gg[1][0] += g0;
                    deme->gg[1][1] += g1;
                    break;
                case 2:
                    deme->gg[1][1] += genotypes[j] * 2;
                    break;
                }
                if (!d->locusUsed[2]) break;     // locus C not used
                switch (ic) {
                case 0:
                    deme->gg[2][0] += genotypes[j] * 2;
                    break;
                case 1:
                    g0 = d->ran->binomial(genotypes[j] * 2, 0.5);
                    g1 = genotypes[j] * 2 - g0;
                    deme->gg[2][0] += g0;
                    deme->gg[2][1] += g1;
                    break;
                case 2:
                    deme->gg[2][1] += genotypes[j] * 2;
                    break;
                }
            }
        }
    }
    deme->nn = (deme->gg[0][0] + deme->gg[0][1]) / 2u;
}


/******************************************************************************
            Generation function: simulate one generation

Discrete non-overlapping generations
mating, reproduction (with drift and possible fecundity selection), mutation, 
viability selection, population regulation, migration
******************************************************************************/

void epistasisGenerationFunction(AltruData * d, int state) {
    int32_t ideme;                                         // deme index
    EpistasisDeme * deme;                                  // point to deme
    int32_t genotypes[27];                                 // individuals of each genotype
    double fitness[27];                                    // fitness of each genotype
    int32_t neighbors[8];                                  // list of neighbor islands
    int numNeighbors;                                      // number of neighbor islands
    int32_t iNeighbor;                                     // index to neighbor islands
    EpistasisDeme * neighborDeme = 0;                      // pointer to neighbor deme

    // check allocated memory
    if (d->demeData == 0) return;

    if (state == state_start) {                            // initialization

        // initialize statistics variables
        statisticsInit0(d);

        // initialize demes
        d->nIslands = d->maxIslands;

        for (ideme = 0; ideme < d->nIslands; ideme++) {
            deme = (EpistasisDeme*)(d->demeData) + ideme;  // point to current deme
            // randomize island size
            int32_t populationSize;
            if (d->carryingCapacityStandardDeviation > 0.) {
                // randomize carrying capacity
                double upperLimit = d->nMaxPerDeme + d->carryingCapacityStandardDeviation * 7.;
                if (upperLimit > 16000.) upperLimit = d->carryingCapacity[0] + 10000.;
                populationSize = (int)lround(d->ran->normalTrunc(d->carryingCapacity[0], d->carryingCapacityStandardDeviation, 0., upperLimit));
                if (populationSize < d->minGroupSize) populationSize = d->minGroupSize;
            }
            else {
                // constant carrying capacity
                populationSize = d->nMaxPerDeme;
            }
            deme->nmax = populationSize;
            // populate island
            deme->nn = d->ran->poisson(deme->nmax / 2);    // random population size
            if (deme->nn > deme->nmax) deme->nn = deme->nmax; // limit population
            if (deme->nn < d->minGroupSize) deme->nn = d->minGroupSize; 
            // random gene pool of island
            for (int locus = 0; locus < d->nLoci; locus++) {
                if (d->locusUsed[locus]) {
                    deme->gg[locus][1] = d->ran->binomial(deme->nn * 2, d->fg0[locus]); // number of wild type genes
                    deme->gg[locus][0] = deme->nn * 2 - deme->gg[locus][1];             // number of mutant genes 
                }
                else {
                    deme->gg[locus][0] = deme->gg[locus][1] = 0;
                }
                deme->mutationAge[locus] = 0;
            }
            errorCheck(d, deme);
        }
        // reset migrant pool
        deme = (EpistasisDeme*)(d->demeData) + d->maxIslands;
        deme->nmax = deme->nn = 0;
        for (int locus = 0; locus < 3; locus++) {
            deme->gg[locus][0] = deme->gg[locus][1] = 0;
            deme->mutationAge[locus] = 0;
        }

        // start running
        d->runState = state_run;
    }

    // calculate fitness of each of the 27 possible genotypes
    fitnessOfGenotypes(d, fitness);

    // simulate one generation 
    if (state == state_run) {

        // reset statictics variables
        //statisticsInit1(d); // don't reset totalPhenotypes
        for (int i = 0; i < d->nLoci; i++) {
            d->mutations[i][0] = d->mutations[i][1] = 0;
            d->genePool[i][0] = d->genePool[i][1] = 0;
        }
        d->migrantsTot = 0;
        d->totalPopulation = 0;
        d->stopCause = stopCauseNone;

        int const numGenotypes = d->locusUsed[locusC] ? 27 : 9; // number of possible genotypees
        int const numUsedLoci  = d->locusUsed[locusC] ? 3 : 2;  // number of loci in use

        // deme loop 1: mutation, growth, selection
        for (ideme = 0; ideme < d->nIslands; ideme++) {
            deme = (EpistasisDeme*)(d->demeData) + ideme;       // point to current deme

            errorCheck(d, deme);

            // mutation
            for (int locus = 0; locus < numUsedLoci; locus++) {
                if (d->locusUsed[locus]) {
                    int32_t forwardMutations  = d->ran->binomial(deme->gg[locus][0], d->murate[locus][0]);
                    int32_t backwardMutations = d->ran->binomial(deme->gg[locus][1], d->murate[locus][1]);
                    deme->gg[locus][1] += forwardMutations - backwardMutations;
                    deme->gg[locus][0] -= forwardMutations - backwardMutations;
                    d->mutations[locus][0] += forwardMutations;
                    d->mutations[locus][1] += backwardMutations;
                    // set mutation age
                    if (forwardMutations && deme->gg[locus][1] && deme->mutationAge[locus] == 0) {
                        deme->mutationAge[locus] = 1; // new mutation
                    }
                }
            }

            // reproduction and growth
            if (d->selectionModel == selectionFecundity) {

                // split deme genes into 27 possible genotype combinations
                genePool2Genotypes(d, deme, genotypes);

                // growth of each genotype
                for (int g = 0; g < numGenotypes; g++) {
                    genotypes[g] = d->ran->poisson(genotypes[g] * fitness[g]);
                }

                // new deme gene pool
                genotypes2GenePool(d, genotypes, deme);
            }
            else {
                errors.reportError("other selection model not supported for epistasis");
            }
        }

        // deme loop 2: migration
        for (ideme = 0; ideme < d->nIslands; ideme++) {
            deme = (EpistasisDeme*)(d->demeData) + ideme;       // point to current deme

            errorCheck(d, deme);

            numNeighbors = findNeighbors(d, ideme, neighbors);  // find neighbor demes

            if (numNeighbors || d->immigrationPattern == immigrationCommonPool) {
                // immigrate
                float migrationRate = d->migrationRate[0];

                int32_t numMigrants = 0;                        // number of immigrants
                if (d->immigrationPattern == immigrationCommonPool) {
                    // migrants come from common pool
                    neighborDeme = (EpistasisDeme*)(d->demeData) + d->maxIslands; // point to common pool
                }
                else if (d->immigrationPattern == immigrationRandomGroup) {
                    int migrantSource = d->ran->iRandom(0, d->maxIslands - 2);
                    if (migrantSource >= ideme) migrantSource++;                  // don't migrate from same deme
                    neighborDeme = (EpistasisDeme*)(d->demeData) + migrantSource;
                }
                else {
                    // migrants come from random neighbor
                    iNeighbor = d->ran->iRandom(0, numNeighbors - 1);             // find random neighbor deme
                    neighborDeme = (EpistasisDeme*)(d->demeData) + neighbors[iNeighbor]; // point to neighbor deme
                }

                switch (d->immigrationPattern) {
                case immigrationCommonPool:      // from common pool
                case immigrationNeighbor:        // from random neighbor deme
                case immigrationRandomGroup:     // from random deme
                    numMigrants = d->ran->poisson(migrationRate * deme->nmax);
                    break;
                case immigrationProportional:    // prop. w. population, from neighbor
                    numMigrants = d->ran->poisson(migrationRate * deme->nn);
                    break;
                case immigrationVacant:          // prop. w. vacant capacity, from neighbor
                    numMigrants = d->ran->poisson(migrationRate * (deme->nmax - deme->nn + 1));
                    break;
                }

                if (numNeighbors > 0 && numMigrants > 0) {
                    if (numMigrants > neighborDeme->nn) numMigrants = neighborDeme->nn;
                    // select migrant genes
                    for (int locus = 0; locus < numUsedLoci; locus++) {
                        int32_t gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborDeme->gg[locus][0], neighborDeme->nn * 2);
                        int32_t gmig1 = numMigrants * 2 - gmig0;
                        deme->gg[locus][0] += gmig0;          // add migrant genes to current deme
                        deme->gg[locus][1] += gmig1;
                        neighborDeme->gg[locus][0] -= gmig0;  // subtract migrant genes from neighbor deme
                        neighborDeme->gg[locus][1] -= gmig1;
                        // transfer mutation age if mutant gene transferred
                        if (gmig1 && neighborDeme->mutationAge[locus] && !deme->mutationAge[locus]) {
                            deme->mutationAge[locus] = neighborDeme->mutationAge[locus];
                        }
                    }
                    deme->nn += numMigrants;
                    neighborDeme->nn -= numMigrants;
                    d->migrantsTot += numMigrants;  // count total number of migrants
                }

                // limit population to carrying capacity
                int32_t nn = deme->nn;
                if (nn > deme->nmax) {
                    for (int locus = 0; locus < numUsedLoci; locus++) {
                        int32_t gn = d->ran->hypergeometric(2 * deme->nmax, deme->gg[locus][0], 2 * nn);
                        int32_t ga = 2 * deme->nmax - gn;
                        deme->gg[locus][0] = gn;
                        deme->gg[locus][1] = ga;
                    }
                    deme->nn = deme->nmax;
                }
            }

            // update mutationAge
            for (int locus = 0; locus < numUsedLoci; locus++) {
                if (deme->gg[locus][1] == 0) deme->mutationAge[locus] = 0; // reset mutationAge if mutation died
            }
            // update mutationAge for each locus until it meets mutation at other locus
            if (deme->mutationAge[locusA] && !deme->mutationAge[locusB]) deme->mutationAge[locusA]++;
            if (deme->mutationAge[locusB] && !deme->mutationAge[locusA]) deme->mutationAge[locusB]++;
            if (deme->mutationAge[locusC] && !(deme->mutationAge[locusA] && deme->mutationAge[locusB])) deme->mutationAge[locusC]++;
        }

        // deme loop 3: update global gene pool and statistics
        for (ideme = 0; ideme < d->nIslands; ideme++) {
            deme = (EpistasisDeme*)(d->demeData) + ideme;     // point to current deme
            for (int locus = 0; locus < numUsedLoci; locus++) {
                d->genePool[locus][0] += deme->gg[locus][0];
                d->genePool[locus][1] += deme->gg[locus][1];
            }
            d->totalPopulation += deme->nn;
        }

        // make migrant pool
        if (d->migrationTopology == topologyCommonPool) {  // migrant pool used
            // migrant pool placed last in memory block
            EpistasisDeme * migrantPool = (EpistasisDeme*)(d->demeData) + d->maxIslands;
            double reduction = 1.;
            if (d->totalPopulation > 100000000) {  // avoid overflow
                reduction = 100000000. / d->totalPopulation;
            }
            for (int locus = 0; locus < numUsedLoci; locus++) {
                migrantPool->gg[locus][0] = int32_t(d->genePool[locus][0] * reduction);
                migrantPool->gg[locus][1] = int32_t(d->genePool[locus][1] * reduction);
            }
            migrantPool->nn = (migrantPool->gg[locusA][0] + migrantPool->gg[locusA][1]) / 2u;
            migrantPool->nmax = migrantPool->nn;
        }

        // update statistics counters
        //statisticsUpdate(d);
        d->generations++;  // count generations
        for (int locus = 0; locus < numUsedLoci; locus++) {
            if (d->locusUsed[locus] && d->totalPopulation > 0) {
                d->geneFraction[locus] = float(d->genePool[locus][1]) / float(d->totalPopulation * 2);
                d->sumMutations += d->mutations[locus][0] + d->mutations[locus][1];
            }
            else {
                d->geneFraction[locus] = 0.f;
            }
        }
        d->sumMigrants += d->migrantsTot;
        // user data for phenotypes
        for (int i = 0; i < 7; i++) d->userData[i+3].f = d->totalPhenotypes[i];
        // user data for mutation age
        deme = (EpistasisDeme*)(d->demeData) + d->nIslands/2;     // point to arbitrary deme
        for (int i = locusA; i <= locusC; i++) d->userData[i].i = deme->mutationAge[i];

        // check if finished
        checkStopCriterion(d);
    }
}
