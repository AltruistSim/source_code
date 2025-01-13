/***************************  epistasis_model.cpp   ***************************
* Author:        Agner Fog
* Date created:  1995-08-13
* Last modified: 2025-01-09
* Version:       3.003
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines a model of epistasis in a structured or viscous population.
* Epistasis can lead to evolution through punctuated equilibria
*
* (c) Copyright 2023-2025 Agner Fog.
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

// Group structure, describing each group or territory
struct EpistasisGroup {
    int32_t nn;              // number of individuals in group
    int32_t nmax;            // carrying capacity or max group size
    int32_t gg[3][2];        // gene pool, locus A, B, C
    int32_t mutationAge[3];  // age of each mutant gene before combined with mutants at other loci
};

const int locusA = 0;        // A locus
const int locusB = 1;        // B locus
const int locusC = 2;        // C locus


/*
struct GroupFieldDescriptor {
    int32_t type;            // 0: end of list, 
                             // 1: size of group structure, 
                             // 2: carrying capacity (max individuals)
                             // 3: population (number of individuals)
                             // 4: gene count, 
                             // 5: genotype or phenotype count
                             // 8: group property
    int32_t varType;         // varInt16 or varInt32 or varFloat
    int32_t offset;          // offset into group structure
    int32_t statistic;       // 1: calculate sum, 2: calculate mean
                             // 4: possible stop criterion, stop when above certain value
                             // 8: possible stop criterion, stop when below certain value
    int32_t graphics;        // show in graphics display:
                             // 2:  max size
                             // 3:  population. divide by max size or nMaxPerGroup to get relative size
                             // 4:  area. divide by territorySizeMax to get relative size
                             // 10: gene count for mutant, primary locus
                             // 11: gene count for mutant, secondary locus
                             // 12: gene count for mutant, third locus
                             // 20: group property
    const char * name;       // name of variable
};
*/

// List of fields in group structure
static const GroupFieldDescriptor epistasisGroupDescriptors[] = {
    {1, varInt32, sizeof(EpistasisGroup), 0, 0, 0},
    {3, varInt32, groupFieldOffset(EpistasisGroup,nn),             0, 3,  "population"},
    {4, varInt32, groupFieldOffset(EpistasisGroup,gg[locusA][0]),  8, 0,  "nonA gene"},
    {4, varInt32, groupFieldOffset(EpistasisGroup,gg[locusA][1]),  5, 10, "A gene"},
    {4, varInt32, groupFieldOffset(EpistasisGroup,gg[locusB][0]),  2, 0,  "nonB gene"},
    {4, varInt32, groupFieldOffset(EpistasisGroup,gg[locusB][1]),  2, 11, "B gene"},
    {4, varInt32, groupFieldOffset(EpistasisGroup,gg[locusC][0]),  2, 0,  "nonC gene"},
    {4, varInt32, groupFieldOffset(EpistasisGroup,gg[locusC][1]),  2, 12, "C gene"},
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
    epistasisGroupDescriptors,
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
        d->modelGroupStructureSize = sizeof(EpistasisGroup);
        d->modelPopOffset = groupFieldOffset(EpistasisGroup,nn);

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

static void errorCheck(AltruData * d, EpistasisGroup * group) {
    // check group data for consistency
    char text[64];
    int igroup = group - (EpistasisGroup*)(d->groupData);

    if (group->gg[locusC][0] + group->gg[locusC][1] != 0 && !d->locusUsed[locusC]) {
        errors.reportError("unused locus != 0");
    }

    if (group->gg[locusA][0] < 0 || group->gg[locusA][1] < 0
        || group->gg[locusA][0] + group->gg[locusA][1] != group->nn * 2
        || (group->gg[locusA][0] + group->gg[locusA][1] & 1)) {
        sprintf_s(text, "error generation %i, group %i, nn %i, ggA %i %i",
            (int)d->generations, igroup,
            group->nn, group->gg[locusA][0], group->gg[locusA][1]);
        errors.reportError(text);
    }
    if ((d->locusUsed[locusB] && group->gg[locusB][0] + group->gg[locusB][1] != group->nn * 2) ||
        (d->locusUsed[locusC] && group->gg[locusC][0] + group->gg[locusC][1] != group->nn * 2)) {
        sprintf_s(text, "error nn %i, ggB %i %i, ggC %i %i",
            group->nn, group->gg[locusB][0], group->gg[locusB][1], group->gg[locusC][0], group->gg[locusC][1]);
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
static void genePool2Genotypes(AltruData * d, EpistasisGroup * group, int32_t genotypes[27]) {
    int32_t genoA[3];                  // genotypes of each combination at A locus
    int32_t genoB[3];                  // genotypes of each combination at B locus
    int32_t genoC[3] = {0,0,0};        // genotypes of each combination at C locus
    int32_t genoAB[9];                 // genotypes of each combination at A and B locus
    int i;

    for (i = 0; i < 27; i++) genotypes[i] = 0;   // reset all genotypes
    if (group->gg[locusC][1] == 0) {
        // shortcut for no C mutants
        if (group->gg[locusB][1] == 0) {
            // shortcut for no B mutants
            if (group->gg[locusA][1] == 0) {
                // no mutants at all
                genotypes[0] = group->gg[locusA][0] / 2u;
            }
            else {
                // only A mutants
                combineGenes(group->gg[locusA], genoA, d->ran);
                for (i = 0; i < 3; i++) genotypes[i] = genoA[i];
            }
        }
        else if (group->gg[locusA][1] == 0) {
            // only B mutants
            combineGenes(group->gg[locusB], genoB, d->ran);
            for (i = 0; i < 3; i++) genotypes[i*3] = genoB[i];
        }
        else {
            // A and B mutants, no C mutants
            // get genotypes of each combination at A locus:
            combineGenes(group->gg[locusA], genoA, d->ran);
            // get genotypes of each combination at B locus:
            combineGenes(group->gg[locusB], genoB, d->ran);
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
        combineGenes(group->gg[locusA], genoA, d->ran);
        // get genotypes of each combination at B locus:
        combineGenes(group->gg[locusB], genoB, d->ran);
        // get genotypes of each combination at C locus:
        combineGenes(group->gg[locusC], genoC, d->ran);
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
    if (group->gg[locusC][1] == 0) {
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
    if (d->locusUsed[locusC] && genoC[1] + genoC[2] && group->nn) {
        // C present
        int c, ac, bc, abc;
        if (d->dominance[locusC] == recessive) {           // C recessive
            c   = phenoAB[0] * genoC[2] / group->nn;        // C only
            ac  = phenoAB[1] * genoC[2] / group->nn;        // C and A
            bc  = phenoAB[2] * genoC[2] / group->nn;        // C and B
            abc = phenoAB[3] * genoC[2] / group->nn;        // ABC
        }
        else {   // C dominant
            c   = phenoAB[0] * (genoC[1] + genoC[2]) / group->nn;
            ac  = phenoAB[1] * (genoC[1] + genoC[2]) / group->nn;
            bc  = phenoAB[2] * (genoC[1] + genoC[2]) / group->nn;
            abc = phenoAB[3] * (genoC[1] + genoC[2]) / group->nn;
        }
        d->totalPhenotypes[0] += phenoAB[0] - c;
        d->totalPhenotypes[1] += phenoAB[1] - ac;
        d->totalPhenotypes[2] += phenoAB[2] - bc;
        d->totalPhenotypes[3] += c;
        d->totalPhenotypes[4] += phenoAB[3] - abc;
        d->totalPhenotypes[5] += abc;
        int other = group->nn - c - ac - bc - abc;
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
static void genotypes2GenePool(AltruData * d, int32_t genotypes[27], EpistasisGroup * group) {
    for (int i = 0; i < 3; i++) group->gg[i][0] = group->gg[i][1] = 0; // reset all
    for (int ia = 0; ia < 3; ia++) {                       // loop through locus A genotypes
        for (int ib = 0; ib < 3; ib++) {                   // loop through locus B genotypes
            for (int ic = 0; ic < 3; ic++) {               // loop through locus C genotypes
                int j = ia + ib * 3 + ic * 9;              // index into genotypes array
                int32_t g0, g1;                            // randomize genes for offspring of heterozygotes
                switch (ia) {
                case 0:
                    group->gg[0][0] += genotypes[j] * 2;
                    break;
                case 1:
                    g0 = d->ran->binomial(genotypes[j] * 2, 0.5);
                    g1 = genotypes[j] * 2 - g0;
                    group->gg[0][0] += g0;
                    group->gg[0][1] += g1;
                    break;
                case 2:
                    group->gg[0][1] += genotypes[j] * 2;
                    break;
                }
                switch (ib) {
                case 0:
                    group->gg[1][0] += genotypes[j] * 2;
                    break;
                case 1:
                    g0 = d->ran->binomial(genotypes[j] * 2, 0.5);
                    g1 = genotypes[j] * 2 - g0;
                    group->gg[1][0] += g0;
                    group->gg[1][1] += g1;
                    break;
                case 2:
                    group->gg[1][1] += genotypes[j] * 2;
                    break;
                }
                if (!d->locusUsed[2]) break;     // locus C not used
                switch (ic) {
                case 0:
                    group->gg[2][0] += genotypes[j] * 2;
                    break;
                case 1:
                    g0 = d->ran->binomial(genotypes[j] * 2, 0.5);
                    g1 = genotypes[j] * 2 - g0;
                    group->gg[2][0] += g0;
                    group->gg[2][1] += g1;
                    break;
                case 2:
                    group->gg[2][1] += genotypes[j] * 2;
                    break;
                }
            }
        }
    }
    group->nn = (group->gg[0][0] + group->gg[0][1]) / 2u;
}


/******************************************************************************
            Generation function: simulate one generation

Discrete non-overlapping generations
mating, reproduction (with drift and possible fecundity selection), mutation, 
viability selection, population regulation, migration
******************************************************************************/

void epistasisGenerationFunction(AltruData * d, int state) {
    int32_t igroup;                                        // group index
    EpistasisGroup * group;                                // point to group
    int32_t genotypes[27];                                 // individuals of each genotype
    double fitness[27];                                    // fitness of each genotype
    int32_t neighbors[8];                                  // list of neighbor islands
    int numNeighbors;                                      // number of neighbor islands
    int32_t iNeighbor;                                     // index to neighbor islands
    EpistasisGroup * neighborGroup = 0;                    // pointer to neighbor group

    // check allocated memory
    if (d->groupData == 0) return;

    if (state == state_start) {                            // initialization

        // initialize statistics variables
        statisticsInit0(d);

        // initialize groups
        d->nIslands = d->maxIslands;

        for (igroup = 0; igroup < d->nIslands; igroup++) {
            group = (EpistasisGroup*)(d->groupData) + igroup;  // point to current group
            // randomize island size
            int32_t populationSize;
            if (d->carryingCapacityStandardDeviation > 0.) {
                // randomize carrying capacity
                double upperLimit = d->nMaxPerGroup + d->carryingCapacityStandardDeviation * 7.;
                if (upperLimit > 16000.) upperLimit = d->nMaxPerGroup + 10000.;
                populationSize = (int)lround(d->ran->normalTrunc(d->nMaxPerGroup, d->carryingCapacityStandardDeviation, 0., upperLimit));
                if (populationSize < d->minGroupSize) populationSize = d->minGroupSize;
            }
            else {
                // constant carrying capacity
                populationSize = d->nMaxPerGroup;
            }
            group->nmax = populationSize;
            // populate island
            group->nn = d->ran->poisson(group->nmax / 2);    // random population size
            if (group->nn > group->nmax) group->nn = group->nmax; // limit population
            if (group->nn < d->minGroupSize) group->nn = d->minGroupSize; 
            // random gene pool of island
            for (int locus = 0; locus < d->nLoci; locus++) {
                if (d->locusUsed[locus]) {
                    group->gg[locus][1] = d->ran->binomial(group->nn * 2, d->fg0[locus]); // number of wild type genes
                    group->gg[locus][0] = group->nn * 2 - group->gg[locus][1];             // number of mutant genes 
                }
                else {
                    group->gg[locus][0] = group->gg[locus][1] = 0;
                }
                group->mutationAge[locus] = 0;
            }
            errorCheck(d, group);
        }
        // reset migrant pool
        group = (EpistasisGroup*)(d->groupData) + d->maxIslands;
        group->nmax = group->nn = 0;
        for (int locus = 0; locus < 3; locus++) {
            group->gg[locus][0] = group->gg[locus][1] = 0;
            group->mutationAge[locus] = 0;
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

        // group loop 1: mutation, growth, selection
        for (igroup = 0; igroup < d->nIslands; igroup++) {
            group = (EpistasisGroup*)(d->groupData) + igroup;       // point to current group

            errorCheck(d, group);

            // mutation
            for (int locus = 0; locus < numUsedLoci; locus++) {
                if (d->locusUsed[locus]) {
                    int32_t forwardMutations  = d->ran->binomial(group->gg[locus][0], d->murate[locus][0]);
                    int32_t backwardMutations = d->ran->binomial(group->gg[locus][1], d->murate[locus][1]);
                    group->gg[locus][1] += forwardMutations - backwardMutations;
                    group->gg[locus][0] -= forwardMutations - backwardMutations;
                    d->mutations[locus][0] += forwardMutations;
                    d->mutations[locus][1] += backwardMutations;
                    // set mutation age
                    if (forwardMutations && group->gg[locus][1] && group->mutationAge[locus] == 0) {
                        group->mutationAge[locus] = 1; // new mutation
                    }
                }
            }

            // reproduction and growth
            if (d->selectionModel == selectionFecundity) {

                // split group genes into 27 possible genotype combinations
                genePool2Genotypes(d, group, genotypes);

                // growth of each genotype
                for (int g = 0; g < numGenotypes; g++) {
                    genotypes[g] = d->ran->poisson(genotypes[g] * fitness[g]);
                }

                // new group gene pool
                genotypes2GenePool(d, genotypes, group);
            }
            else {
                errors.reportError("other selection model not supported for epistasis");
            }
        }

        // group loop 2: migration and population control
        for (igroup = 0; igroup < d->nIslands; igroup++) {
            group = (EpistasisGroup*)(d->groupData) + igroup;       // point to current group

            errorCheck(d, group);

            numNeighbors = findNeighbors(d, igroup, neighbors);  // find neighbor groups

            if (numNeighbors || d->immigrationPattern == immigrationCommonPool) {
                // immigrate
                float migrationRate = d->migrationRate[0];

                int32_t numMigrants = 0;                        // number of immigrants
                if (d->immigrationPattern == immigrationCommonPool) {
                    // migrants come from common pool
                    neighborGroup = (EpistasisGroup*)(d->groupData) + d->maxIslands; // point to common pool
                }
                else if (d->immigrationPattern == immigrationRandomGroup) {
                    int migrantSource = d->ran->iRandom(0, d->maxIslands - 2);
                    if (migrantSource >= igroup) migrantSource++;                  // don't migrate from same group
                    neighborGroup = (EpistasisGroup*)(d->groupData) + migrantSource;
                }
                else {
                    // migrants come from random neighbor
                    iNeighbor = d->ran->iRandom(0, numNeighbors - 1);             // find random neighbor group
                    neighborGroup = (EpistasisGroup*)(d->groupData) + neighbors[iNeighbor]; // point to neighbor group
                }

                switch (d->immigrationPattern) {
                case immigrationCommonPool:      // from common pool
                case immigrationNeighbor:        // from random neighbor group
                case immigrationRandomGroup:     // from random group
                    //numMigrants = d->ran->poisson(migrationRate * group->nmax);
                    numMigrants = d->ran->binomial(group->nmax, migrationRate);
                    break;
                case immigrationProportional:    // prop. w. population, from neighbor
                    //numMigrants = d->ran->poisson(migrationRate * group->nn);
                    numMigrants = d->ran->binomial(group->nn, migrationRate);
                    break;
                case immigrationVacant:          // prop. w. vacant capacity, from neighbor
                    //numMigrants = d->ran->poisson(migrationRate * (group->nmax - group->nn + 1));
                    if (group->nn <= group->nmax) {
                        numMigrants = d->ran->binomial(group->nmax - group->nn + 1, migrationRate);
                    }
                    break;
                }

                if (numNeighbors > 0 && numMigrants > 0) {
                    if (numMigrants > neighborGroup->nn) numMigrants = neighborGroup->nn;
                    // select migrant genes
                    for (int locus = 0; locus < numUsedLoci; locus++) {
                        int32_t gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborGroup->gg[locus][0], neighborGroup->nn * 2);
                        int32_t gmig1 = numMigrants * 2 - gmig0;
                        group->gg[locus][0] += gmig0;          // add migrant genes to current group
                        group->gg[locus][1] += gmig1;
                        neighborGroup->gg[locus][0] -= gmig0;  // subtract migrant genes from neighbor group
                        neighborGroup->gg[locus][1] -= gmig1;
                        // transfer mutation age if mutant gene transferred
                        if (gmig1 && neighborGroup->mutationAge[locus] && !group->mutationAge[locus]) {
                            group->mutationAge[locus] = neighborGroup->mutationAge[locus];
                        }
                    }
                    group->nn += numMigrants;
                    neighborGroup->nn -= numMigrants;
                    d->migrantsTot += numMigrants;  // count total number of migrants
                }

                // limit population to carrying capacity
                int32_t nn = group->nn;
                if (nn > group->nmax) {
                    for (int locus = 0; locus < numUsedLoci; locus++) {
                        int32_t gn = d->ran->hypergeometric(2 * group->nmax, group->gg[locus][0], 2 * nn);
                        int32_t ga = 2 * group->nmax - gn;
                        group->gg[locus][0] = gn;
                        group->gg[locus][1] = ga;
                    }
                    group->nn = group->nmax;
                }
            }

            // update mutationAge
            for (int locus = 0; locus < numUsedLoci; locus++) {
                if (group->gg[locus][1] == 0) group->mutationAge[locus] = 0; // reset mutationAge if mutation died
            }
            // update mutationAge for each locus until it meets mutation at other locus
            if (group->mutationAge[locusA] && !group->mutationAge[locusB]) group->mutationAge[locusA]++;
            if (group->mutationAge[locusB] && !group->mutationAge[locusA]) group->mutationAge[locusB]++;
            if (group->mutationAge[locusC] && !(group->mutationAge[locusA] && group->mutationAge[locusB])) group->mutationAge[locusC]++;
        }

        // group loop 3: update global gene pool and statistics
        for (igroup = 0; igroup < d->nIslands; igroup++) {
            group = (EpistasisGroup*)(d->groupData) + igroup;     // point to current group
            for (int locus = 0; locus < numUsedLoci; locus++) {
                d->genePool[locus][0] += group->gg[locus][0];
                d->genePool[locus][1] += group->gg[locus][1];
            }
            d->totalPopulation += group->nn;
        }

        // make migrant pool
        if (d->migrationTopology == topologyCommonPool) {  // migrant pool used
            // migrant pool placed last in memory block
            EpistasisGroup * migrantPool = (EpistasisGroup*)(d->groupData) + d->maxIslands;
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
        group = (EpistasisGroup*)(d->groupData) + d->nIslands/2;     // point to arbitrary group
        for (int i = locusA; i <= locusC; i++) d->userData[i].i = group->mutationAge[i];

        // check if finished
        checkStopCriterion(d);
    }
}
