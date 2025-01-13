/**************************  regality_model.cpp   *****************************
* Author:        Agner Fog
* Date created:  2023-10-03
* Last modified: 2025-01-11
* Version:       3.003
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the regality model where groups can steal territory
* from each other, and the military strength depends on common support for
* a strong leader. 
* See: Fog, A: Warlike and Peaceful Societies. Open Book Publishers, 2017.
* https://doi.org/10.11647/obp.0128
*
* (c) Copyright 2023-2025 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"


// List of model-specific parameters
// These parameters will appear in parameter list and dialog box.
// Maximum number of model-specific parameters is maxModelSpecificParameters defined in altruist.h
// ParameterDef::type can be: 1: integer, 2: float, 3: boolean
// ParameterDef::num must be 1
// Offsets and name must be unique
static const ParameterDef regalityParameterDefinitions[] {
    {0, 1, 0, "Model-specific parameters"},
    {3, 1, altruDataOffset(modelspec_i[0]), "leader_regality_counts"},
    {0, 0, 0, 0}
};

// list of indexes into userData
static int const numUser = 0;                    // number of userData

const int locusRegality = locusAltruism;         // locus for regality gene

/*
struct GroupFieldDescriptor {
    int32_t type;            // 0: end of list, 
                             // 1: size of group structure, 
                             // 2: carrying capacity (max individuals)
                             // 3: population (number of individuals)
                             // 4: gene count, 
                             // 5: genotype or phenotype count
                             // 6: area
                             // 8: group property
    int32_t varType;         // varInt16 or varInt32 or varFloat
    int32_t offset;          // offset into group structure
    int32_t statistic;       // 1: calculate sum, 2: calculate mean
                             // 4: possible stop criterion, stop when above certain value
                             // 8: possible stop criterion, stop when below certain value
    int32_t graphics;        // show in graphics display:
                             // 2:  max size
                             // 3:  population. divide by max size or nMaxPerGroup to get relative size
                             // 6:  area. divide by territorySizeMax to get relative size
                             // 10: gene count for mutant, primary locus
                             // 11: gene count for mutant, secondary locus
                             // 12: gene count for mutant, third locus
                             // 20: group property
    const char * name;       // name of variable
}; */

// List of fields in group structure
static const GroupFieldDescriptor regalityGroupDescriptors[] = {
    {1, varInt32, sizeof(TerriGroup), 0, 0, 0},
    {3, varInt32, groupFieldOffset(TerriGroup,nn),           0, 3,  "population"},
    {2, varInt32, groupFieldOffset(TerriGroup,nmax),         0, 2,  "carrying capacity"},
    {6, varInt32, groupFieldOffset(TerriGroup,area),         2, 6,  "area"},
    {4, varInt32, groupFieldOffset(TerriGroup,gAltruism[0]), 8, 10, "neutral gene"},
    {4, varInt32, groupFieldOffset(TerriGroup,gAltruism[1]), 5, 10, "regality gene"},
    {8, varInt16, groupFieldOffset(TerriGroup,age),          0, 20, "age"},
    {8, varFloat, groupFieldOffset(TerriGroup,groupfit),     2, 20, "group fitness"},
    {8, varInt32, groupFieldOffset(TerriGroup,emigrationPotential), 0, 20, "emigration potential"},
    {0, 0, 0, 0, 0, 0}   // mark end of list
};

// function prototypes
void regalityInitFunction(AltruData * d, int index);
void regalityGenerationFunction(AltruData * d, int mode);


// define model
static const ModelDescriptor regalityModel = {
    "Regality model",
    "Groups have dynamic territories. A strong group can conquer territory from a weaker neighbor group. "
    "Group strength depends on common support for a strong leader.",
    regalityParameterDefinitions,
    regalityGroupDescriptors,
    &regalityInitFunction,
    &regalityGenerationFunction
};

// Register the regality model in the global model descriptor list:
Construct regalityConstructor(regalityModel);


// allocate memory buffers
static void allocateMemory(AltruData * d) {
    // allocate group list
    if (d->bufferSize < (d->maxIslands + 1) * sizeof(TerriGroup)) {
        if (d->bufferSize) delete[] d->groupData; // delete any previous smaller buffer
        d->groupData = new int8_t[(d->maxIslands + 1) * sizeof(TerriGroup)]();
        if (d->groupData == 0) errors.reportError("memory allocation failed");
        else d->bufferSize = (d->maxIslands + 1) * sizeof(TerriGroup);
    }
    // allocate point map
    if (d->extraBufferSize[ebMap] < d->totArea * sizeof(idType)) {
        if (d->extraBuffer[ebMap]) delete[] d->extraBuffer[ebMap]; // delete any previous smaller buffer
        d->extraBuffer[ebMap] = new int8_t[d->totArea * sizeof(idType)]();
        if (d->extraBuffer[ebMap] == 0) errors.reportError("memory allocation failed");
        else d->extraBufferSize[ebMap] = d->totArea * sizeof(idType);
    }
    // allocate for point list
    if (d->extraBufferSize[ebPointsSorted] < pointListSize * sizeof(TPoint)) {
        if (d->extraBuffer[ebPointsSorted]) delete[] d->extraBuffer[ebPointsSorted]; // delete any previous smaller buffer
        d->extraBuffer[ebPointsSorted] = new int8_t[pointListSize * sizeof(TPoint)]();
        d->extraBufferSize[ebPointsSorted] = pointListSize * sizeof(TPoint);
    }
    // allocate for point distance list
    if (d->extraBufferSize[ebPointDistSorted] < pointListSize * sizeof(TPointDist)) {
        if (d->extraBuffer[ebPointDistSorted]) delete[] d->extraBuffer[ebPointDistSorted]; // delete any previous smaller buffer
        d->extraBuffer[ebPointDistSorted] = new int8_t[pointListSize * sizeof(TPointDist)]();
        d->extraBufferSize[ebPointDistSorted] = pointListSize * sizeof(TPointDist);
    }
    // allocate order list
    if (d->extraBufferSize[ebOrder] < (d->maxIslands + 1) * sizeof(int32_t)) {
        if (d->extraBuffer[ebOrder]) delete[] d->extraBuffer[ebOrder]; // delete any previous smaller buffer
        d->extraBuffer[ebOrder] = new int8_t[(d->maxIslands + 1) * sizeof(int32_t)]();
        d->extraBufferSize[ebOrder] = (d->maxIslands + 1) * sizeof(int32_t);
    }
}


/******************************************************************************
*         initialization function
* state values:
* 1: model changed
* 2: parameters changed
******************************************************************************/

void regalityInitFunction(AltruData * d, int state) {
    switch (state) {
    case model_initialize:  // initialization
        // model version
        d->modelVersionMajor = 3;
        d->modelVersionMinor = 0;
        d->modelGroupStructureSize = sizeof(TerriGroup);
        d->modelPopOffset = groupFieldOffset(TerriGroup, nn); 
        // enable and disable fields in parameters dialog boxes
        d->bGeographyParametersUsed = 0x510D1;
        d->bGroupPropertiesUsed = 0x300D0;

        d->bSelModels = 0x17; // don't support independent viability selection                
        if (d->selectionModel == selectionIndependentViability) d->selectionModel = selectionFecundity;

        // names of genes and loci
        d->sLocusName[locusRegality]          = "regality";
        d->sGeneName[locusRegality][0]        = "neutral";
        d->sGeneName[locusRegality][1]        = "regality";
        d->sPhenotypeName[locusRegality][0]   = "neutral";
        d->sPhenotypeName[locusRegality][1]   = "regalist";
        d->sPhenotypeName[locusRegality+1][0] = "among neutral";
        d->sPhenotypeName[locusRegality+1][1] = "among regalist";

        d->bEmigrationPattern = 0b1111;                    // emigration patterns supported
        d->bImmigrationPattern = 0b111111;                 // immigration patterns supported
        d->bExtinctionPatterns = 0;
        d->bWarPatterns = 0b0111;
        d->bStopCriterionUsed = 0x1F;

        d->nLoci = 1;                                      // number of loci
        for (int i = 0; i < maxLoci; i++) {
            d->locusUsed[i] = i < d->nLoci;
        }

        d->graphicsTypeForModel = graphicsTerritories;
        d->nUser = numUser;
        for (int i = 0; i < d->nUser; i++) d->userData[i].i = 0;  // reset model-specific variables
        d->userDataNames = 0;
        d->numRows = 0;                                    // make sure area is recalculated
        d->fitRowLocus = locusRegality;
        d->fitColLocus = locusRegality;
        d->fitSupColLocus = locusConformity;
        d->fitConditionName = 1;
        if (d->territorySizeMax < d->territorySizeMin * 2 + 1) {
            // make sure territorySizeMax > 2*territorySizeMin so that big territories can be split in two
            if (d->territorySizeMax < 8) d->territorySizeMax = 8;
            d->territorySizeMin = d->territorySizeMax / 2u - 2;        
        }
        checkArea(d);                                      // calculate area
        allocateMemory(d);                                 // allocate memory
        break;
    case model_parameters_changed:  // parameters changed
        // check that area is not changed while running
        if (d->runState > state_start && d->totArea != d->lastArea) {
            errors.reportError("Error: Cannot change area while running");
            d->totArea = d->lastArea;
        }
        checkArea(d);                                      // calculate area
        if (d->maxIslands < d->totArea / d->territorySizeMin) d->maxIslands = d->totArea / d->territorySizeMin;
        allocateMemory(d);  // re-allocate memory in case parameters have been changed
        break;    
    }
}

// check group data for consistency
static void errorCheck(AltruData * d, TerriGroup * group) {
    char text[64];
    if (group->gAltruism[0] < 0 || group->gAltruism[1] < 0
        || group->gAltruism[0] + group->gAltruism[1] != group->nn * 2
        || (group->gAltruism[0] + group->gAltruism[1] & 1)) {
        sprintf_s(text, "error n %i, gAltru %i %i",
            group->nn, group->gAltruism[0], group->gAltruism[1]);
        errors.reportError(text);
    }
    if (group->nn > 0 && group->area == 0) {
        errors.reportError("Population on zero area");
    }
}

// find an empty TerriGroup record
static TerriGroup * findEmptyGroup(AltruData * d) {
    for (int i = 1; i <= d->maxIslands; i++) {
        TerriGroup * group = (TerriGroup*)(d->groupData) + i;
        if (group->area == 0 && group->nn == 0) {
            group->age = 0; 
            group->id = idType(i);
            if (i > d->nIslands) d->nIslands = i;          // update number of territories in use
            return group;
        }
    }
    return 0;  // not found
}



/******************************************************************************
            Generation function: simulate one generation

Discrete non-overlapping generations:
Loop 1:
mutation, immigration, 
select a leader,
calculate group fitness,
reproduction with growth, selection, and drift,
population regulation,
calculate group fitness again,
territorial war
Loop 2:
split big group
Loop 3:
statistics

******************************************************************************/

void regalityGenerationFunction(AltruData * d, int state) {
    int32_t iGroup;                                        // group index
    TerriGroup * group;                                    // point to group
    int32_t genotypes1[6];                                 // [0-2]: wild type, heterozygote, regalists before selection, [3-5]; same, for leader
    int32_t genotypes2[6];                                 // same, after selection
    double fitness[6];                                     // [0-2]: fitness of each genotype, [3-5]; same, for leader
    double survivalRate[3] = {0};                          // survival rate for each genotype under independent viability selection only
    TerriGroup * neighborGroup = 0;                        // pointer to neighbor group
    idType * map = (idType *)(d->extraBuffer[ebMap]);      // map of owner of each point
    Habitat terrain(d);                                    // object for manipulating points and territories
    float fractionOfRegality;                              // fracton of regality phenotype individuals


    if (state == state_start) {
        // initialization before first generation

        checkArea(d);                                      // calculate area of terrain
        allocateMemory(d);                                 // re-allocate memory in case parameters have been changed
        map = (idType *)(d->extraBuffer[ebMap]);           // map of owner of each point
        terrain.setData(d);                                // data have changed

        // initialize statistics variables
        statisticsInit0(d);

        // reset group array
        memset(d->groupData, 0, (d->maxIslands + 1) * sizeof(TerriGroup));

        // make a fence around habitat 
        // (this makes the code simpler because we don't have to check if we go outside the limits of x and y)
        int i, j;
        for (i = 0; i < d->totArea; i++) map[i] = 1;
        for (i = 0; i < d->rowLengthTerri; i++) map[i] = 0;
        for (j = 0; j < d->totArea; j += d->rowLengthTerri) {
            map[j] = map[j + d->rowLengthTerri - 1] = 0;
        }
        for (i = (d->numRows - 1) * d->rowLengthTerri; i <d->totArea; i++) map[i] = 0;

        // first make one big territory covering the whole habitat except the fence
        group = (TerriGroup*)(d->groupData) + 1;              // group 0 used for fence, first group is number 1
        group->area = (d->rowLengthTerri - 2) * (d->numRows - 2);
        group->sx = (d->rowLengthTerri-1)*group->area/2;
        group->sy = (d->numRows-1)*group->area/2;
        group->centerx = group->sx / group->area; 
        group->centery = group->sy / group->area; 
        group->id = 1;
        d->nIslands = 1;

        // make territories by splitting the big territory into many small ones
        int32_t size1 = d->territorySizeMax / 2;           // size of new territories
        j = 1;
        int k;
        while (j) {  // keep splitting all territories that are too big
            j = 0;
            int n = d->nIslands;
            for (k = 1; k <= n; k++) {
                group = (TerriGroup*)(d->groupData) + k;      // point to current group
                int h = group->area / (2 * size1);
                if (h > 1) {
                    // split this territory
                    TerriGroup * group2 = findEmptyGroup(d);  // find empty record and update nIslands
                    if (group2 == 0) {
                        errors.reportError("Not enough vacant territory records");
                        j = 0;  break;
                    }
                    terrain.splitArea(group, group2, h * size1);
                    j++;
                }
            }
        }

        // put population into all territories
        for (k = 1; k <= d->nIslands; k++) {
            group = (TerriGroup*)(d->groupData) + k;          // point to current group
            int32_t nn = (int32_t)lround(group->area * d->carryingCapacity[0]);
            if (nn < d->minGroupSize) nn = d->minGroupSize;
            group->nn = group->nmax = nn;
            group->gAltruism[1] = d->ran->binomial(nn*2, d->fg0[locusRegality]);
            group->gAltruism[0] = nn*2 - group->gAltruism[1];
            group->age = 0;
            group->emigrationPotential = 1;                 // emigrationPotential is unused
            group->groupfit = 0.f;                          // not calculated yet
        }
        // start running
        d->runState = state_run;
    }

    // initialize statistics counters
    d->migrantsTot = 0;
    d->groupsDied = 0;
    d->newColonies = 0;
    for (int i = 0; i < d->nLoci; i++) {
        d->mutations[i][0] = d->mutations[i][1] = 0;    
    }

    // war gain standard deviation
    double warGainStdDev = d->warIntensity * d->territorySizeMax * 0.1;

    // shuffle groups in random order
    int32_t * const orderList = (int32_t*)d->extraBuffer[ebOrder];      // list in random order
    d->ran->shuffle(orderList, 1, d->nIslands);                         // make list in random order. exclude number 0

    if (d->carryingCapacity[locusAltruism] <= 0.f) {
        d->totalPopulation = 0;
        d->stopCause = stopCauseDied;
        d->simulationResult = resultDied;
        d->runState = state_stop;
        return;
    }


    /**************************************************************************
    *                       First group loop:                                  *
    *  select groups in random order
    *  mutation
    *  find neighbors
    *  immigration from neighbors
    *  calculate group fitness
    *  select a leader
    *  differential growth, selection, drift
    *  limit population to carrying capacity of territory
    *  calculate group fitness again
    *  territorial war
    **************************************************************************/
    for (int32_t ter = 0; ter < d->nIslands; ter++) {
        // loop through groups in random order 
        iGroup = orderList[ter];                            // random order
        group = (TerriGroup*)(d->groupData) + iGroup;          // point to current group
        if (group->area == 0) continue;                     // skip unused records

        // mutation
        int32_t forwardMutations  = d->ran->binomial(group->gAltruism[0], d->murate[locusRegality][0]);
        int32_t backwardMutations = d->ran->binomial(group->gAltruism[1], d->murate[locusRegality][1]);
        group->gAltruism[1] += forwardMutations - backwardMutations;
        group->gAltruism[0] -= forwardMutations - backwardMutations;
        d->mutations[locusRegality][0] += forwardMutations;
        d->mutations[locusRegality][1] += backwardMutations;

        // find neighbors
        Neighbor neighborlist[maxNeighbors];
        int numNeighbors;
        // note: In the special case that a territory has an enclave, the findNeighbors function 
        // may find only the enclave or it may find all neighbors except the enclave.
        // It would be too time consuming to make a thorough search to find both enclaves and other neighbors.
        // The consequence is that a group may sometimes fail to fight war against an enclave
        // even if it would be able to conquer it completely.
        // If you want to fix this problem, you may identify territories with only one neighbor as enclaves.
        terrain.findNeighbors(group, neighborlist, &numNeighbors);
        // calculate total border length
        int32_t perimeter = 0;
        for (int i = 0; i < numNeighbors; i++) {
            perimeter += neighborlist[i].sharedBorder;
        }

        // immigration: pick immigrants from all neighbors
        int sNeighbor = 0;
        if (d->immigrationPattern != immigrationAllNeighbors) {
            // select a neighbor to migrate from. Probability proportional with shared border
            int n = d->ran->iRandom(0, perimeter-1);
            int s = 0;
            for (sNeighbor = 0; sNeighbor < numNeighbors; sNeighbor++) {
                s += neighborlist[sNeighbor].sharedBorder;
                if (s > n) break;
            }
        }
        TerriGroup globalGenePool;                         // used in case immigrationCommonPool
        TerriGroup * neighborGroup;                        // neighbor where immigrants come from
        float expectedMigrants = 0.f;                      // expected number of migrants before randomization
        float stdGroupPopulation = d->territorySizeMax * d->carryingCapacity[locusAltruism]; // population of a big group = d->nMaxPerGroup
        int numMigrants = 0;                               // number of immigrants
        int maxMigrants = 0;                               // max number of immigrants = population of migration source

        // loop for all neighbors
        int iNeighbor;
        for (iNeighbor = 0; iNeighbor < numNeighbors; iNeighbor++) {

            // select one or all neighbors
            if (d->immigrationPattern != immigrationAllNeighbors && iNeighbor != sNeighbor) continue;
            neighborGroup = (TerriGroup*)(d->groupData) + neighborlist[iNeighbor].id;

            if (neighborGroup->nn > 0) {
                // select immigration pattern
                switch (d->immigrationPattern) {
                case immigrationCommonPool:                // immigration from global gene pool
                    if (d->totalPopulation > INT32_MAX) errors.reportError("Global gene pool overflow");
                    globalGenePool.nn = d->totalPopulation;
                    globalGenePool.gAltruism[0] = d->genePool[locusRegality][0];
                    globalGenePool.gAltruism[1] = d->genePool[locusRegality][1];
                    globalGenePool.emigrationPotential = d->nMaxPerGroup;
                    neighborGroup = &globalGenePool;
                    //expectedMigrants = d->migrationRate[0] * group->nmax;
                    expectedMigrants = d->migrationRate[0];
                    maxMigrants = group->nmax;
                    break;
                case immigrationNeighbor:                  // immigration from one random neighbor
                    //expectedMigrants = d->migrationRate[0] * group->nmax * neighborGroup->nn / stdGroupPopulation; * neighborGroup->emigrationPotential / stdGroupPopulation;
                    expectedMigrants = d->migrationRate[0] * group->nmax / stdGroupPopulation;
                    maxMigrants = neighborGroup->emigrationPotential;
                    break;
                case immigrationAllNeighbors:              // immigration from all neighbors, proportional to shared border
                    //expectedMigrants = d->migrationRate[0] * group->nmax * neighborGroup->emigrationPotential * neighborlist[iNeighbor].sharedBorder / (stdGroupPopulation * perimeter);
                    expectedMigrants = d->migrationRate[0] * group->nmax
                        * neighborlist[iNeighbor].sharedBorder / (stdGroupPopulation * perimeter);
                    maxMigrants = neighborGroup->emigrationPotential;
                    break;
                case immigrationProportional:              // immigration from one random neighbor, proportional to population of target group
                    //expectedMigrants = d->migrationRate[0] * group->nn * neighborGroup->emigrationPotential / stdGroupPopulation;
                    expectedMigrants = d->migrationRate[0] * group->nn / stdGroupPopulation;
                    maxMigrants = neighborGroup->emigrationPotential;
                    break;
                case immigrationVacant:                    // immigration from one random neighbor, proportional to vacant area
                    //expectedMigrants = d->migrationRate[0] * (group->nmax - group->nn) * neighborGroup->emigrationPotential / stdGroupPopulation;
                    expectedMigrants = d->migrationRate[0] * (group->nmax - group->nn) / stdGroupPopulation;
                    if (expectedMigrants < 0) expectedMigrants = 0;
                    maxMigrants = neighborGroup->emigrationPotential;
                    break;
                case immigrationRandomGroup:               // immigration from a random group
                    {
                    int migrantSource = d->ran->iRandom(1, d->nIslands - 1);
                    if (migrantSource >= iGroup) migrantSource++;     // avoid migration from same group
                    neighborGroup = (TerriGroup*)(d->groupData) + migrantSource; 
                    //expectedMigrants = d->migrationRate[0] * group->nmax * neighborGroup->emigrationPotential / stdGroupPopulation;
                    expectedMigrants = d->migrationRate[0] * group->nmax / stdGroupPopulation;
                    maxMigrants = neighborGroup->emigrationPotential;
                    }
                }
                if (expectedMigrants > 0.f) { 
                    if (expectedMigrants > 1.f) expectedMigrants = 1.f;
                    numMigrants = d->ran->binomial(maxMigrants, expectedMigrants); // make number of migrants random
                }
                if (numMigrants > neighborGroup->nn) numMigrants = neighborGroup->nn; // cannot migrate more than there are
                // select migrant genes
                int32_t gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborGroup->gAltruism[0], neighborGroup->nn * 2);
                int32_t gmig1 = numMigrants * 2 - gmig0;
                group->gAltruism[0] += gmig0;               // add migrant genes to current group
                group->gAltruism[1] += gmig1;
                group->nn += numMigrants;
                neighborGroup->gAltruism[0] -= gmig0;       // subtract migrant genes from neighbor group
                neighborGroup->gAltruism[1] -= gmig1;
                neighborGroup->nn -= numMigrants;
                d->migrantsTot += numMigrants;             // count total number of migrants
            }
        }

        // split into genotypes
        combineGenes(group->gAltruism, genotypes1, d->ran);

        // select a leader
        for (int i = 3; i < 6; i++) genotypes1[i] = 0; // reset leader genotypes
        if (d->leaderSelection < 0.f) d->leaderSelection = 0.f;
        bool hasLeader = true;                                        // true if there is a leader
        if (group->nn == 0) {                                         // group is empty
            hasLeader = false;
        }
        else if (d->leaderSelection > 0.99E6f) {                     // leader must have regality
            if (genotypes1[2] == 0 && (d->dominance[locusRegality] != dominant || genotypes1[1] == 0)) {
                hasLeader = false;                                   // cannot make a leader because leader must be regalist and there are none
            }
        }
        else if (d->leaderSelection == 0.f) {                        // leader must be non-regal
            if (genotypes1[0] == 0 && (d->dominance[locusRegality] == dominant || genotypes1[1] == 0)) {
                hasLeader = false;                                   // there are non non-regal
            }
        }
        if (hasLeader) {
            if (d->leaderSelection == 1.f) {
                // unbiased leader selection among available genotypes
                int32_t n = d->ran->iRandom(0, group->nn - 1); // pick a random individual among all individuals
                int32_t s = 0;
                for (int i = 0; i < 3; i++) {
                    s += genotypes1[i];
                    if (s > n) {
                        genotypes1[i]--; genotypes1[i + 3]++;  // move selected individual to leader genotype array
                        break;
                    }
                }
            }
            else {
                // biased leader selection
                float weightedList[3];  // genotypes * (bias for being a leader)
                weightedList[0] = genotypes1[0];
                weightedList[2] = genotypes1[2] * d->leaderSelection;
                weightedList[1] = genotypes1[1];
                if (d->dominance[locusRegality] == dominant) weightedList[1] *= d->leaderSelection;
                else if (d->dominance[locusRegality] == incompleteDominant) weightedList[1] *= sqrtf(d->leaderSelection);
                float w = 0;  // sum of weightedList
                for (int i = 0; i < 3; i++) w += weightedList[i];
                float n = d->ran->randomf() * w;  // select random from weightedList
                float s = 0.f;
                int i;
                for (i = 0; i < 3; i++) {
                    s += weightedList[i];
                    if (s > n) break;
                }
                if (i < 3) {
                    genotypes1[i]--; genotypes1[i + 3]++;  // move selected individual to leader genotype array
                }
            }
        }

        bool leaderRegalityCounts = d->modelspec_i[0] != 0 && hasLeader;

        // calculate fraction of regality
        float regalPhenotypes = genotypes1[2];             // number of regality phenotypes
        switch (d->dominance[locusRegality]) {
        case recessive:                // regalism recessive
            if (leaderRegalityCounts) regalPhenotypes += genotypes1[5];
            break;
        case dominant:                 // regalism dominant
            regalPhenotypes += genotypes1[1];
            if (leaderRegalityCounts) regalPhenotypes += genotypes1[4];
            break;
        case incompleteDominant:       // half dominant
            regalPhenotypes += genotypes1[1] * 0.5f;
            if (leaderRegalityCounts) regalPhenotypes += genotypes1[4] * 0.5f;
        }
        if (!leaderRegalityCounts && hasLeader) regalPhenotypes++; // leader always supports himself

        fractionOfRegality = 0.f;
        if (group->nn > 0) {  // avoid division by 0          
            fractionOfRegality = regalPhenotypes / float(group->nn);
        }

        // calculate group fitness from fraction of regality and group fitness exponent
        switch (d->fitfunc) {
        case 0:              // one regalist is enough
            if (regalPhenotypes > 0.f) group->groupfit = 1.f;
            else group->groupfit = 0.f;
            break;
        case 1: case 3:      // convex/concave
            group->groupfit = powf(fractionOfRegality, d->groupFitCurvature);
            break;
        case 2:  default:    // linear
            group->groupfit = fractionOfRegality;
            break;
        case 4:              // all or nothing
            if (regalPhenotypes == group->nn) group->groupfit = 1.f;
            else group->groupfit = 0.f;
            break;
        }

        // calculate relative fitness for each phenotype
        fitness[0] = d->fit[0] * (1.f - group->groupfit) + d->fit[1] * group->groupfit; // fitness of neutral
        fitness[2] = d->fit[2] * (1.f - group->groupfit) + d->fit[3] * group->groupfit; // fitness of regality
        switch (d->dominance[locusRegality]) {   // find fitness of heterozygotes
        case recessive:                // regalism recessive
            fitness[1] = fitness[0];
            break;
        case dominant:                 // regalism dominant
            fitness[1] = fitness[2];
            break;
        case incompleteDominant:       // half dominant
            fitness[1] = 0.5 * (fitness[0] + fitness[2]);
            break;
        }
        // fitness of leader
        if (hasLeader) {
            for (int i = 0; i < 3; i++) {
                double fi = 1.f + d->leaderAdvantage * fractionOfRegality; // fitness factor of leader
                double fn = 1.;                                            // fitness factor of non-leaders
                if (fi < 0.) fi = 0.;
                if (group->nn > 1) {        // avoid division by 0                
                    fn = (group->nn - fi) / (group->nn - 1.);
                    if (fn < 0.) {
                        fn = 0.;  fi = group->nn;
                    }
                }
                else {
                    fi = 1.;
                }
                fitness[i+3] = fitness[i] * fi;
                fitness[i]   = fitness[i] * fn;
            }
        }
        else {
            for (int i = 3; i < 6; i++) {
                fitness[i] = 0.;
            }
        }

        // emigration
        switch(d->emigrationPattern) {
        case emigrationConstant:                           // emigration independent of population
            group->emigrationPotential = (d->nMaxPerGroup + d->minGroupSize) / 2;
            break;
        case emigrationProportional:                       // emigration proportional to population
            group->emigrationPotential = group->nn;
            break;
        case emigrationExcess:                             // emigration proportional to excess population
            group->emigrationPotential = (int32_t)lround(group->nn * d->growthRate - group->nmax);
            if (group->emigrationPotential < 1) group->emigrationPotential = 1;
            break;
        case emigrationGroupfit:                           // emigration proportional to group fitness * population
            group->emigrationPotential = group->nn * group->groupfit;
            break;
        }

        // growth and selection
        int32_t populationSize = 0;
        if (d->selectionModel == selectionFecundity) {
            // fecundity selection
            for (int i = 0; i < 3; i++) {
                // growth of each genotype, non-leader + leader
                genotypes2[i] = d->ran->poisson(genotypes1[i] * fitness[i] + genotypes1[i+3] * fitness[i+3]);
                populationSize += genotypes2[i];
                genotypes2[i+3] = 0;  // non-leaders and leader genotypes joined in genotypes2[0-2]
            }
            if (populationSize > group->nmax) {
                if (d->emigrationPattern == emigrationExcess) {
                    group->emigrationPotential = populationSize - group->nmax;
                }
                populationSize = group->nmax; // min(populationSize, nmax)
            }
            // reduce population to carrying capacity
            d->ran->multiHypergeometric(genotypes2, genotypes2, populationSize, 3);
        }
        else {
            // viability selection
            for (int i = 0; i < 6; i++) {
                // equal growth of all genotypes before selection
                populationSize += genotypes1[i] = d->ran->poisson(genotypes1[i] * d->growthRate);
            }
            if (populationSize > group->nmax) populationSize = group->nmax; // min(populationSize, nmax)

            // selection according to selection model
            switch (d->selectionModel) {

            case selectionFecundity:             // fecundity selection handled above
                break;

            case selectionPositiveViability:     // positive viability selection
                d->ran->multiWalleniusNCHyp(genotypes2, genotypes1, fitness, populationSize, 6);
                break;

            case selectionNegativeViability:     // negatitive viability selection
                d->ran->complMultiWalleniusNCHyp(genotypes2, genotypes1, fitness, populationSize, 6);
                break;

            case selectionMinimumViability:      // minimum viability selection                                
                d->ran->multiFishersNCHyp(genotypes2, genotypes1, fitness, populationSize, 6);
                break;
            }        
        }

        // get gene pool of group from genotypes
        int32_t nn = 0;
        for (int i = 0; i < 3; i++) {
            genotypes2[i] += genotypes2[i+3];    // add leader and non-leader genotypes
            genotypes2[i+3] = 0;
            nn += genotypes2[i];                 // count all
        }

#if true
        // drift in offspring of heterozygotes
        int32_t g0 = d->ran->binomial(genotypes2[1]*2, 0.5);
#else
        int32_t g0 = genotypes2[1];
#endif
        int32_t g1 = genotypes2[1]*2 - g0;

        // sum gene pool
        group->nn = nn;
        group->gAltruism[0] = genotypes2[0]*2 + g0;
        group->gAltruism[1] = genotypes2[2]*2 + g1;


        // the rest of this loop is only relevant for war
        if (d->warIntensity == 0.f) continue;

        // re-calculate fraction of regality
        regalPhenotypes = genotypes2[2] + genotypes2[5];         // number of regality phenotypes
        switch (d->dominance[locusRegality]) {
        case recessive:                // regalism recessive
            break;
        case dominant:                 // regalism dominant
            regalPhenotypes += genotypes2[1] + genotypes2[4];
            break;
        case incompleteDominant:       // half dominant
            regalPhenotypes += (genotypes2[1] + genotypes2[4]) * 0.5f;
        }
        fractionOfRegality = 0.f;
        if (group->nn) {                // avoid division by 0          
            fractionOfRegality = regalPhenotypes / float(group->nn);
        }

        // re-calculate group fitness from fraction of regalists and group fitness exponent
        switch (d->fitfunc) {
        case 0:                        // one regalist is enough
            if (regalPhenotypes > 0.f) group->groupfit = 1.f;
            else group->groupfit = 0.f;
            break;
        case 1: case 3:                // convex/concave
            group->groupfit = powf(fractionOfRegality, d->groupFitCurvature);
            break;
        case 2:  default:              // linear
            group->groupfit = fractionOfRegality;
            break;
        case 4:                        // all or nothing
            if (regalPhenotypes == group->nn) group->groupfit = 1.f;
            else group->groupfit = 0.f;
            break;
        }

        // territorial war
        switch (d->warPattern) {
        case warAgainstAll:
            break;
        case warAgainstBorder: {
            // select a neighbor to attack. Probability proportional with shared border
            int n = d->ran->iRandom(0, perimeter-1);
            int s = 0;
            for (sNeighbor = 0; sNeighbor < numNeighbors; sNeighbor++) {
                s += neighborlist[sNeighbor].sharedBorder;
                if (s > n) break;
            }
            break;}
        case warAgainstWeak: {
            // select a neighbor to attack. Probability inverse proportional with strength
            float strength, strength1, sumStrength = 0.f;
            int i;
            for (i = 0; i < numNeighbors; i++) {           // calculate strength of all neighbors
                neighborGroup = (TerriGroup*)(d->groupData) + neighborlist[i].id;
                strength = neighborGroup->groupfit * neighborGroup->nn;
                if (strength > 0.f) sumStrength += 1.f / strength;
            }
            strength1 = d->ran->randomf() * sumStrength;   // random point in sum of strengths
            sumStrength = 0.f;
            for (i = 0; i < numNeighbors; i++) {
                neighborGroup = (TerriGroup*)(d->groupData) + neighborlist[i].id;
                strength = neighborGroup->groupfit * neighborGroup->nn;
                if (strength > 0.f) sumStrength += 1.f / strength;
                else if (neighborGroup->area > 0) break;    // neighbor with zero strength and nonzero area
                if (sumStrength >= strength1) break;
            }
            sNeighbor = i;
            break;}
        }

        // loop for all neighbors:
        for (iNeighbor = 0; iNeighbor < numNeighbors; iNeighbor++) {

            // select one or all neighbors
            if (d->warPattern != warAgainstAll && iNeighbor != sNeighbor) continue;
            neighborGroup = (TerriGroup*)(d->groupData) + neighborlist[iNeighbor].id;

            // attack this neighbor. expected gain is warIntensity * (difference in strength)
            float expectedGain = d->warIntensity * d->territorySizeMax *
                (group->groupfit * group->nn - neighborGroup->groupfit * neighborGroup->nn) / d->nMaxPerGroup;
            if (expectedGain > 0.f) {
                float gain = d->ran->normal(expectedGain, warGainStdDev);
                int transferArea = (int)lround(gain);

                if (transferArea > 0) {  // ignore negative gains
                    // note: if this code is modified so that transferArea can be negative (i.e. we can lose 
                    // territory to neighbor) and warPattern == warAgainstAll, then we have to 
                    // check here if neighborGroup is still a neighbor to group (use isNeighbor function)

                    int32_t neighborAreaBefore = neighborGroup->area;

                    // take land from neighborGroup
                    terrain.conquer(neighborGroup, group, transferArea, neighborlist + iNeighbor);

                    // transfer survivors from loser to winner
                    transferArea = neighborAreaBefore - neighborGroup->area; // may occasionally transfer more than expected
                    // expected number of migrating survivors is proportional to the lost area
                    float meanSurvivors = d->surviv * transferArea / neighborAreaBefore;
                    if (meanSurvivors > 1.f) meanSurvivors = 1.f;  // can happen if surviv > 1
                    int32_t survivors = d->ran->binomial(neighborGroup->nn, meanSurvivors);
                    if (survivors > 0) {
                        // genes of survivors
                        int32_t tg0 = d->ran->hypergeometric(survivors * 2, neighborGroup->gAltruism[0], neighborGroup->nn * 2);
                        int32_t tg1 = survivors * 2 - tg0;
                        group->gAltruism[0] += tg0;             // add survivor genes to current group
                        group->gAltruism[1] += tg1;
                        group->nn += survivors;
                        neighborGroup->gAltruism[0] -= tg0;     // subtract survivor genes from neighbor group
                        neighborGroup->gAltruism[1] -= tg1;
                        neighborGroup->nn -= survivors;
                        d->migrantsTot += survivors;           // count total number of migrants
                        // transfer emigration potential
                        int32_t emiPot = neighborGroup->emigrationPotential * transferArea / neighborAreaBefore;
                        group->emigrationPotential += emiPot;
                        neighborGroup->emigrationPotential -= emiPot;
                    }
                    if (neighborGroup->area == 0) {
                        // total loss
                        d->groupsDied++;
                        neighborGroup->nn = 0;
                        neighborGroup->gAltruism[0] = 0;
                        neighborGroup->gAltruism[1] = 0;
                    }
                }
            }
        }
    }

    /**************************************************************************
    *              second group loop: split big groups                         *
    **************************************************************************/
    // loop through groups in linear order 
    for (iGroup = 1; iGroup <= d->nIslands; iGroup++) {
        group = (TerriGroup*)(d->groupData) + iGroup;          // point to current group
        if (group->area == 0) continue;                     // skip unused record

        // split groups with large area in two
        if (group->area > d->territorySizeMax && group->nn > 4) {
            int e;
            // find a vacant group record
            neighborGroup = findEmptyGroup(d);               // find empty record and update nIslands
            if (neighborGroup == 0) {
                errors.reportError("Not enough vacant territory records");
                break;
            }
            terrain.splitArea(group, neighborGroup, group->area / 2);
            d->newColonies++;
            neighborGroup->age = 0;

            // transfer half of population to new area
            if (group->nn > 1) {
                double mean = group->nn * 0.5;
                double limit1 = d->minGroupSize * 0.5;  if (limit1 >= mean) limit1 = 0.;
                double limit2 = group->nn - limit1;  if (limit2 <= mean) limit2 = group->nn;
                int32_t colonists = d->ran->normalTrunc(mean, group->nn * 0.1, limit1, limit2);
                // transfer half of gene pool
                int32_t tg0 = d->ran->hypergeometric(colonists * 2, group->gAltruism[0], group->nn * 2);
                int32_t tg1 = colonists * 2 - tg0;
                neighborGroup->gAltruism[0] = tg0;          // add survivor genes to current group
                neighborGroup->gAltruism[1] = tg1;
                neighborGroup->nn = colonists;
                group->gAltruism[0] -= tg0;                 // subtract survivor genes from neighbor group
                group->gAltruism[1] -= tg1;
                group->nn -= colonists;
            }
            // copy fitness. This approximate group fitness is only for display and statistics.
            // An accurate value will be calculated in next generation
            neighborGroup->groupfit = group->groupfit;
            // transfer half of emigration potential
            group->emigrationPotential = neighborGroup->emigrationPotential = group->emigrationPotential / 2u;
        }
    }

    /**************************************************************************
    *                 third group loop: statistics                             *
    **************************************************************************/

    // reset statistics counters before each generation
    for (int i = 0; i < d->nLoci; i++) {
        d->totalPhenotypes[i] = 0;
        d->genePool[i][0] = d->genePool[i][1] = 0;
    }
    d->inhabitedGroups = 0;
    d->totalPopulation = 0;
    d->altruistGroups = 0;
    d->stopCause = stopCauseNone;

    // loop through groups in linear order 
    for (iGroup = 1; iGroup <= d->nIslands; iGroup++) {
        group = (TerriGroup*)(d->groupData) + iGroup;          // point to current group
        if (group->area == 0) continue;                     // skip unused record
        if (group->nn == 0) continue;                       // skip empty territory
        group->age++;  if (group->age <= 0) group->age--;     // avoid overflow
        //if (group->nn == 0) emptyTerritories++;
        else d->inhabitedGroups++;
        d->totalPopulation += group->nn;
        d->genePool[locusRegality][0] += group->gAltruism[0];
        d->genePool[locusRegality][1] += group->gAltruism[1];
        d->totalPhenotypes[locusRegality] += int(group->groupfit * group->nn);
        if (group->gAltruism[1] > group->gAltruism[0]) d->altruistGroups++;
        
        errorCheck(d, group); // temporary error check. May be removed
    }
    // update statistics counters
    d->generations++;  // count generations
    d->sumMutations += d->mutations[locusRegality][0] + d->mutations[locusRegality][1];
    d->sumExtinctions += d->groupsDied;
    d->sumMigrants += d->migrantsTot;

    // check if finished
    checkStopCriterion(d);

    // temporary debug check
    //checkAllTerritories(d); // temporary check
}
