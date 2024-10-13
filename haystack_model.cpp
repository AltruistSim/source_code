/****************************  haystack_model.cpp   *****************************
* Author:        Agner Fog
* Date created:  2023-09-04
* Last modified: 2024-10-13
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the haystack model, also known as intrademic group selection
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

void haystackInitFunction(AltruData * d, int index);
void haystackGenerationFunction(AltruData * d, int mode);


/******************************************************************************
*         structures and definitions
******************************************************************************/

// List of model-specific parameters
// These parameters will appear in parameter list and dialog box.
// Maximum number of model-specific parameters is maxModelSpecificParameters defined in altruist.h
// ParameterDef::type can be: 1: integer, 2: float, 3: boolean
// ParameterDef::num must be 1 if not an array
// Offsets and names must be unique
static const ParameterDef haystackParameterDefinitions[] {
    {0, 1, 0, "Model-specific parameters"},
    {1, 1, altruDataOffset(modelspec_i[0]), "max_metapopulation"},
    {2, 1, altruDataOffset(modelspec_f[1]), "mix_period_growthrate"},
    {0, 0, 0, 0}
};

// Group structure, describing each group or territory
struct HaystackGroup {
    int32_t nmax;            // carrying capacity or max group size
    int32_t nn;              // number of individuals in group
    int32_t gAltruism[2];    // gene pool, egoism, altruism
    int32_t age;             // group age
    int32_t altruists;       // number of phenotypic altruists
    int32_t emigrationPotential; // excess available for emigration
    float   groupfit;        // group fitness
};

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
static const GroupFieldDescriptor haystackGroupDescriptors[] = {
    {1, varInt32, sizeof(HaystackGroup), 0, 0, 0},
    {2, varInt32, groupFieldOffset(HaystackGroup,nmax),           0, 2,  "carrying capacity"},
    {3, varInt32, groupFieldOffset(HaystackGroup,nn),             0, 3,  "population"},
    {4, varInt32, groupFieldOffset(HaystackGroup,gAltruism[0]),   8, 0,  "egoism gene"},
    {4, varInt32, groupFieldOffset(HaystackGroup,gAltruism[1]),   5, 10, "altruism gene"},
    {8, varInt32, groupFieldOffset(HaystackGroup,age),            0,  0, "age"},
    {8, varFloat, groupFieldOffset(HaystackGroup,groupfit),       2, 20, "fraction of altruists"},
    {8, varInt32, groupFieldOffset(TerriGroup,emigrationPotential), 0, 20, "emigration potential"},    
    {0, 0, 0, 0, 0, 0}   // mark end of list
};


static const ModelDescriptor haystackModel = {
    "Haystack model",
    "Intrademic group selection."
    " Haystacks are colonized periodically. Population is mixed when haystacks are removed. "
    "Atruist have lower individual fitness, but increase fitness of all members of the same haystack.",
    haystackParameterDefinitions,
    haystackGroupDescriptors,
    &haystackInitFunction,
    &haystackGenerationFunction
};

// This will register the Island model in the global model descriptor list
Construct haystackConstructor(haystackModel);

// values for haystackState
static const int stateMixing = 0;
static const int stateHaystack = 1;


/******************************************************************************
*         initialization function
* state values:
* 1: model changed
* 2: parameters changed
******************************************************************************/

void haystackInitFunction(AltruData * d, int state) {
    switch (state) {
    case model_initialize:  // initialization
        // model version
        d->modelVersionMajor = 3;
        d->modelVersionMinor = 0;
        d->modelGroupStructureSize = sizeof(HaystackGroup);
        d->modelPopOffset = groupFieldOffset(HaystackGroup,nn); // make nn available to findNeighbors function
        d->nLoci = 1;

        // names of genes and loci
        d->sLocusName[locusAltruism]          = "altruism";
        d->sGeneName[locusAltruism][0]        = "egoism";
        d->sGeneName[locusAltruism][1]        = "altruism";
        d->sPhenotypeName[locusAltruism][0]   = "egoist";
        d->sPhenotypeName[locusAltruism][1]   = "altruist";
        d->sPhenotypeName[locusAltruism+1][0] = "among egoists";
        d->sPhenotypeName[locusAltruism+1][1] = "among altruists";

        // enable and disable fields in parameters dialog boxes
        d->bGeographyParametersUsed = 0x3D1206;
        d->bIndividualPropertiesUsed = 0x0F;
        d->bGroupPropertiesUsed = 0x377;
        d->bEmigrationPattern = 0b1111;
        d->bImmigrationPattern = 0b1111111;
        d->bColonPat = 0b1110011;
        d->bExtinctionPatterns = 0;
        d->bWarPatterns = 0;
        d->bTopology = 0b01101100;
        d->bSelModels = 0b10111;       // disable independent viability selection
        d->bStopCriterionUsed = 0x1F;
        d->graphicsTypeForModel = graphicsIslands;
        d->nUser = 3;                  // model-specific variables
        for (int i = 0; i < d->nUser; i++) d->userData[i].i = 0;  // reset model-specific variables
        d->userDataNames = 0;
        d->fitRowLocus = locusAltruism;
        d->fitColLocus = locusAltruism;
        d->fitSupColLocus = locusConformity;
        d->fitConditionName = 1;
        d->locusUsed[1] = false;
        d->locusUsed[2] = false;
        break;
    case model_parameters_changed:     // parameters changed
        d->nIslands = d->maxIslands;
        break;    
    }
}

static void errorCheck(AltruData * d, HaystackGroup * group) {
    // check group data for consistency
        char text[64];
    if (group->gAltruism[0] < 0 || group->gAltruism[1] < 0 
        || group->gAltruism[0] + group->gAltruism[1] != group->nn * 2 
        || (group->gAltruism[0] + group->gAltruism[1] & 1)) {
        sprintf_s(text, "error n %i, gAltru %i %i",
            group->nn, group->gAltruism[0], group->gAltruism[1]);
        errors.reportError(text);
    }
}


/******************************************************************************
            Generation function: simulate one generation

Discrete non-overlapping generations.
mating, reproduction with drift and fecundity selection, mutation, 
population regulation.
Haystack period: growth in isolated groups. growth rate depends on fraction of 
altruists in group. Group extinction possible.
Mixing period: Groups mix up into metapopulation. Growth and population regulation
******************************************************************************/

void haystackGenerationFunction(AltruData * d, int state) {
    int32_t iGroup;                                        // group index
    HaystackGroup * group;                                 // point to group
    int32_t genotypes1[3];                                 // egoists, heterozygote, altruists before selection
    int32_t genotypes2[3];                                 // egoists, heterozygote, altruists after selection
    double fitness[3];                                     // fitness of each genotype
    int32_t neighbors[8];                                  // list of neighbor islands
    int numNeighbors;                                      // number of neighbor islands
    int32_t iNeighbor;                                     // index to neighbor islands
    int strongestNeighbor;                                 // index to strongest neighbor group
    HaystackGroup * neighborGroup = 0;                     // pointer to neighbor group
    double survivalRate[3] = {0};                          // survival rate for each genotype under independent viability selection only
    int32_t populationSize;                                // total population size


    const int inhabitedIslandsLength = 4096;
    int inhabitedIslandsNum;                               // number of inhabited islands in inhabitedIslands list
    // (may differ from d->inhabitedGroups if inhabitedIslands list is full)
    int32_t inhabitedIslands[inhabitedIslandsLength];      // used for list of inhabited islands

    if (d->groupData == 0) return;

    int32_t * const haystackState = &d->userData[0].i;     // 0: mixing period, 1: haystack period
    int32_t * const mixGenerations = &d->userData[1].i;    // number of generations currently spent in current mixing period 
    int32_t * const haystackGenerations = &d->userData[2].i;// number of generations currently spent in current haystack period 
    int32_t * const maxMetaPopulation = &d->modelspec_i[0];// carrying capacity in mixing state
    float   * const mixPeriodGrowthRate = &d->modelspec_f[1]; // growth rate during mixing period

    if (state == state_start) {
        // reset all groups
        memset(d->groupData, 0, (d->maxIslands+1) * sizeof(HaystackGroup));

        // populate group 0 as mixing population
        group = (HaystackGroup*)(d->groupData) + 0;           // point to group 0
        populationSize = d->colonySize * d->maxIslands;
        group->nmax = d->nMaxPerGroup;
        if (populationSize > 1000) populationSize = 1000;
        group->nn = populationSize;
        group->gAltruism[1] = d->ran->binomial(group->nn * 2, d->fg0[locusAltruism]); // number of altruism genes
        group->gAltruism[0] = group->nn * 2 - group->gAltruism[1]; // number of egoism genes 
        group->age = 1;                                     // don't set age = 0, that means extinction
        group->emigrationPotential = 1;
        group->groupfit = 0.f;                              // not calculated yet
        errorCheck(d, group);

        // initialize other variables
        int i;
        statisticsInit0(d);
        d->totalPopulation = populationSize;
        d->inhabitedGroups = 1;
        *haystackState = stateMixing;
        *mixGenerations = 1;
        d->nIslands = 1;

        // reset migrant pool
        group = (HaystackGroup*)(d->groupData) + d->maxIslands;
        group->nmax = group->nn = 0;
        group->gAltruism[0] = group->gAltruism[1] = 0;
        group->age = 0; 
        group->groupfit = 0.f;
        group->emigrationPotential = 0;

        // start running
        d->runState = state_run;
    }

    // simulate one generation 
    if (state == state_run) {
        if (*haystackState == stateMixing) {
            d->nIslands = 1;
            group = (HaystackGroup*)(d->groupData) + 0;       // point to group 0
            if (*mixGenerations == 0) {
                // first generation in mixing state
                // populate group 0 as mixing population
                float reduction = 1.f;
                if (d->totalPopulation > 1000000) {        // limit population size to 1000000
                    reduction = 1000000.f / d->totalPopulation;
                    group->gAltruism[0] = d->genePool[locusAltruism][0] * reduction;
                    group->gAltruism[1] = d->genePool[locusAltruism][1] * reduction;
                    group->nn = (group->gAltruism[0] + group->gAltruism[1]) / 2;
                }
                else {
                    group->gAltruism[0] = d->genePool[locusAltruism][0];
                    group->gAltruism[1] = d->genePool[locusAltruism][1];
                    group->nn = d->totalPopulation;
                }
            }

            // mutation
            int32_t forwardMutations = d->ran->binomial(group->gAltruism[0], d->murate[locusAltruism][0]);
            int32_t backwardMutations = d->ran->binomial(group->gAltruism[1], d->murate[locusAltruism][1]);
            group->gAltruism[1] += forwardMutations - backwardMutations;
            group->gAltruism[0] -= forwardMutations - backwardMutations;
            d->mutations[locusAltruism][0] += forwardMutations;
            d->mutations[locusAltruism][1] += backwardMutations;

            // growth in mixing state without selection, with drift
            if (group->nn > 1) {
                populationSize = d->ran->poisson(double(group->nn) * *mixPeriodGrowthRate);
                if (populationSize < d->minGroupSize) populationSize = d->minGroupSize;
                int32_t me = d->ran->binomial(populationSize * 2, group->gAltruism[0] / (group->nn * 2.));
                int32_t ma = populationSize * 2 - me;
                group->gAltruism[0] = me;
                group->gAltruism[1] = ma;
                group->nn = populationSize;

                // limit population size, with drift
                if (populationSize > *maxMetaPopulation) {
                    me = d->ran->hypergeometric(*maxMetaPopulation * 2, group->gAltruism[0], populationSize * 2);
                    ma = *maxMetaPopulation * 2 - me;
                    group->gAltruism[0] = me;
                    group->gAltruism[1] = ma;
                    group->nn = *maxMetaPopulation;
                }
            }
            group->age = *mixGenerations;

            // statistics
            statisticsInit1(d);
            d->inhabitedGroups = d->altruistGroups = 0;
            d->totalPopulation = group->nn;
            d->genePool[locusAltruism][0] = group->gAltruism[0];
            d->genePool[locusAltruism][1] = group->gAltruism[1];
            d->totalPhenotypes[locusAltruism] = group->altruists;
            group->groupfit = 0.f;
            if (group->nn > 0) {
                d->inhabitedGroups = 1;
                group->groupfit = group->gAltruism[1] / (group->nn * 2.0f);
                if (group->groupfit > 0.5f) d->altruistGroups = 1;
            }

            if (++*mixGenerations >= d->mixingPeriod) {
                // finished mixing state. switch to haystack state
                *haystackState = stateHaystack;
                *haystackGenerations = 0;
            }
        }
        else {
            // haystack state
            if (*haystackGenerations == 0) {
                // first generation in haystack state
                d->nIslands = d->maxIslands;

                int32_t ge, ga, nn;
                group = (HaystackGroup*)(d->groupData) + 0;   // mixing population
                ge = group->gAltruism[0];  ga = group->gAltruism[1];
                nn = group->nn;
                // group loop 0: populate groups
                for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
                    group = (HaystackGroup*)(d->groupData) + iGroup; // point to current group
                    group->nn = d->colonySize;
                    group->nmax = d->nMaxPerGroup;
                    int32_t me = d->ran->binomial(group->nn * 2, ge / (nn * 2.));
                    int32_t ma = group->nn * 2 - me;
                    group->gAltruism[0] = me;
                    group->gAltruism[1] = ma;
                    group->age = 1;
                    group->emigrationPotential = 0;
                    group->groupfit = 0.f;
                }
            }

            // reset statictics variables
            inhabitedIslandsNum = 0;  d->groupsDied = 0;
            d->mutations[0][0] = d->mutations[0][1] = 0;

            // group loop 1: mutation, growth, selection, extinction
            for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
                group = (HaystackGroup*)(d->groupData) + iGroup; // point to current group

                errorCheck(d, group);

                // mutation
                int32_t forwardMutations = d->ran->binomial(group->gAltruism[0], d->murate[locusAltruism][0]);
                int32_t backwardMutations = d->ran->binomial(group->gAltruism[1], d->murate[locusAltruism][1]);
                group->gAltruism[1] += forwardMutations - backwardMutations;
                group->gAltruism[0] -= forwardMutations - backwardMutations;
                d->mutations[locusAltruism][0] += forwardMutations;
                d->mutations[locusAltruism][1] += backwardMutations;

                // reproduction and growth
                if (d->selectionModel != selectionFecundity) {  // same growth rate for all except under fecundity selection
                    int32_t nn2 = d->ran->poisson(double(group->nn) * d->growthRate);    // new population size
                    int32_t me = d->ran->binomial(nn2 * 2, group->gAltruism[0] / (group->nn * 2.));
                    int32_t ma = nn2 * 2 - me;
                    group->gAltruism[0] = me;
                    group->gAltruism[1] = ma;
                }
                group->nn = (group->gAltruism[0] + group->gAltruism[1]) / 2u;

                // calculate fraction of altruists
                combineGenes(group->gAltruism, genotypes1, d->ran);
                group->altruists = genotypes1[2];      // number of phenotypic altruists
                switch (d->dominance[locusAltruism]) {
                case recessive:                       // altruism recessive
                    break;
                case dominant:                        // altruism dominant
                    group->altruists += genotypes1[1];
                    break;
                case incompleteDominant:               // half dominant
                    group->altruists += genotypes1[1] / 2;
                }
                float fractionOfAltruists = 0.f;
                if (group->nn) {                       // avoid division by 0          
                    fractionOfAltruists = float(group->altruists) / float(group->nn);
                }

                // calculate group fitness from fraction of altruists and group fitness exponent
                switch (d->fitfunc) {
                case 0:                               // one altruist is enough
                    if (group->altruists > 0) group->groupfit = 1.f;
                    else group->groupfit = 0.f;
                    break;
                case 1: case 3:                       // convex/concave
                    group->groupfit = powf(fractionOfAltruists, d->groupFitCurvature);
                    break;
                case 2:  default:                     // linear
                    group->groupfit = fractionOfAltruists;
                    break;
                case 4:                               // all or nothing
                    if (group->altruists == group->nn) group->groupfit = 1.f;
                    else group->groupfit = 0.f;
                    break;
                }

                // calculate relative fitness for egoists and altruists
                fitness[0] = d->fit[0] * (1.f - group->groupfit) + d->fit[1] * group->groupfit; // fitness of egoists
                fitness[2] = d->fit[2] * (1.f - group->groupfit) + d->fit[3] * group->groupfit; // fitness of altruists
                switch (d->dominance[0]) {  // find fitness of heterozygotes
                case 0: // altruism recessive
                    fitness[1] = fitness[0];  break;
                case 1: // altruism dominant
                    fitness[1] = fitness[2];  break;
                case 2: // half dominant
                    fitness[1] = 0.5 * (fitness[0] + fitness[2]);  break;
                }

                int32_t populationSize = group->nn;         // new population size
                if (populationSize > group->nmax) populationSize = group->nmax;
                int32_t populationExcess = 0;              // excess population

                // selection according to selection model
                switch (d->selectionModel) {
                case selectionFecundity: default:   // fecundity selection
                    differentialGrowth(group->gAltruism, fitness, genotypes1, d->ran);
                    for (int i = 0; i < 3; i++) {
                        genotypes2[i] = genotypes1[i];
                    }
                    populationSize = genotypes1[0] + genotypes1[1] + genotypes1[2];
                    if (populationSize > group->nmax) {
                        populationExcess = populationSize - group->nmax;
                        if (populationExcess < 0) populationExcess = 0;
                        //populationSize = group->nmax;
                    }
                    // wait with reducing the population size until after migration
                    group->emigrationPotential = populationExcess;
                    break;

                case selectionPositiveViability: // positive viability selection
                    d->ran->multiWalleniusNCHyp(genotypes2, genotypes1, fitness, populationSize, 3);
                    break;

                case selectionNegativeViability: // negatitive viability selection
                    d->ran->complMultiWalleniusNCHyp(genotypes2, genotypes1, fitness, populationSize, 3);
                    break;

                case selectionMinimumViability:  // minimum viability selection                                
                    d->ran->multiFishersNCHyp(genotypes2, genotypes1, fitness, populationSize, 3);
                    break;
                }

                // emigration
                switch (d->emigrationPattern) {
                case emigrationConstant:                           // emigration independent of population
                    group->emigrationPotential = d->nMaxPerGroup / 2u;
                    break;
                case emigrationProportional:                       // emigration proportional to population
                    group->emigrationPotential = group->nn;
                    break;
                case emigrationExcess:                             // emigration proportional to excess population
                    if (d->selectionModel != selectionFecundity) {
                        // (emigrationPotential has been calculated above if selectionModel == selectionFecundity)
                        for (int i = 0; i < 3; i++) {
                            populationExcess += (int32_t)lround(genotypes1[i] * fitness[i]);
                        }
                        populationExcess -= group->nmax;
                        if (populationExcess < 0) populationExcess = 0;
                        group->emigrationPotential = populationExcess;
                    }
                    break;
                case emigrationGroupfit:
                    group->emigrationPotential = (int32_t)lround(group->nn * group->groupfit);
                    break;
                }

                // drift in offspring of heterozygotes
                int32_t g0 = d->ran->binomial(genotypes2[1]*2, 0.5);
                int32_t g1 = genotypes2[1]*2 - g0;

                // genes after selection
                group->gAltruism[0] = genotypes2[0] * 2 + g0;
                group->gAltruism[1] = genotypes2[2] * 2 + g1;
                group->nn = (group->gAltruism[0] + group->gAltruism[1]) / 2;  // number of individuals
                // calculate fraction of altruists after selection
                group->altruists = genotypes2[2]; // number of phenotypic altruists
                switch (d->dominance[locusAltruism]) {
                case 0: // altruism recessive
                    break;
                case 1: // altruism dominant
                    group->altruists += genotypes2[1];
                    break;
                case 2: // half dominant
                    group->altruists += genotypes2[1] / 2;
                }

                //errorCheck(d, group);

                // extinction
                float relativeSize = float(group->nn) / float(group->nmax); // relative size of group
                if (relativeSize > 1.f) relativeSize = 1.0f;

                float extinctionRate =
                    (1.f - group->groupfit) * ((1.f - relativeSize) * d->extinctionRate[0] + relativeSize * d->extinctionRate[1])
                    + group->groupfit * ((1.f - relativeSize) * d->extinctionRate[2] + relativeSize * d->extinctionRate[3]);

                // kill with probability extinctionRate
                bool kill = d->ran->bernoulli(extinctionRate);

                if (kill) {  // this group goes extinct
                    d->groupsDied++;
                    group->emigrationPotential = 0;
                    group->age = 0;
                    if (d->surviv > 0.f) {  // a fraction survives
                        int32_t numSurvivors = int32_t(d->surviv * group->nn);
                        group->gAltruism[0] = d->ran->hypergeometric(numSurvivors * 2, group->gAltruism[0], group->nn * 2);
                        group->gAltruism[1] = numSurvivors * 2 - group->gAltruism[0];
                        group->nn = numSurvivors;
                        group->groupfit = 0.f;
                    }
                    else {   // all die
                        group->gAltruism[0] = group->gAltruism[1] = 0;
                        group->nn = 0;
                        group->groupfit = 0.f;
                    }
                }

                // make list of inhabited islands/haystacks 
                if (group->nn > 0 && inhabitedIslandsNum < inhabitedIslandsLength) {
                    inhabitedIslands[inhabitedIslandsNum++] = iGroup;
                }
            }

            // group loop 2: migration or colonization
            for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
                group = (HaystackGroup*)(d->groupData) + iGroup;     // point to current group

                //errorCheck(d, group);

                // make list of neighbors
                numNeighbors = findNeighbors(d, iGroup, neighbors);
                neighborGroup = 0;
                float migrationRate = d->migrationRate[0];
                HaystackGroup combinedNeighbors;
                int migrationOrColonizationPattern;  // immigration and colonization use similar patterns
                if (group->age == 0 || group->nn == 0) {  // group is extinct or empty. recolonize
                    group->age = 0;                      // reset age if group is empty but not extinct
                    migrationOrColonizationPattern = d->colonizationPattern;
                }
                else {   // group is not extinct. calculate immigration
                    migrationOrColonizationPattern = d->immigrationPattern;
                }

                neighborGroup = 0;
                if (numNeighbors) {
                    // immigration or colonization
                    switch (migrationOrColonizationPattern) {
                    case immigrationCommonPool:
                        // migrants come from common pool
                        neighborGroup = (HaystackGroup*)(d->groupData) + d->maxIslands; // point to common pool
                        break;
                    case immigrationNeighbor:
                    case immigrationProportional:                        
                    case immigrationVacant: 
                        // migrants come from random neighbor
                        iNeighbor = d->ran->iRandom(0, numNeighbors - 1);   // find random neighbor group
                        neighborGroup = (HaystackGroup*)(d->groupData) + neighbors[iNeighbor];  // point to neighbor group
                        break;
                    case immigrationAllNeighbors: {// sum of all neighbors
                        memset(&combinedNeighbors, 0, sizeof(combinedNeighbors));  // set to all 0
                        double sumGene0 = 0, sumGene1 = 0, sumEmigr = 0;
                        for (int i = 0; i < numNeighbors; i++) {
                            neighborGroup = (HaystackGroup*)(d->groupData) + neighbors[i];  // point to neighbor group
                            sumGene0 += double(neighborGroup->gAltruism[0]) * neighborGroup->emigrationPotential;
                            sumGene1 += double(neighborGroup->gAltruism[1]) * neighborGroup->emigrationPotential;
                            sumEmigr += neighborGroup->emigrationPotential;
                        }
                        combinedNeighbors.gAltruism[0] = lround(sumGene0 / sumEmigr);
                        combinedNeighbors.gAltruism[1] = lround(sumGene1 / sumEmigr);
                        uint32_t numg = combinedNeighbors.gAltruism[0] + combinedNeighbors.gAltruism[1];
                        if (numg & 1) { // number of genes cannot be odd. make even
                            combinedNeighbors.gAltruism[0] &= -2;
                            combinedNeighbors.gAltruism[1] &= -2;
                        }
                        combinedNeighbors.nn = numg / 2u;
                        combinedNeighbors.emigrationPotential = sumEmigr / numNeighbors;
                        neighborGroup = &combinedNeighbors;
                        break;}
                    case immigrationRandomGroup:
                        if (inhabitedIslandsNum > 0) {  // immigration from random inhabited island
                            int i = d->ran->iRandom(0, inhabitedIslandsNum - 1);         // find random inhabited island
                            iNeighbor = inhabitedIslands[i];
                            neighborGroup = (HaystackGroup*)(d->groupData) + iNeighbor;     // point to random group
                        }
                        break;
                    case immigrationStrongNeighbor: {
                        // choose among all neighbors with probability proportional to emigrationPotential
                        int32_t sumEmig = 0;  // sum of emigration potential of all neighbors
                        for (int i = 0; i < numNeighbors; i++) {
                            neighborGroup = (HaystackGroup*)(d->groupData) + neighbors[i];  // point to neighbor group
                            sumEmig += neighborGroup->emigrationPotential;
                        }
                        int32_t choose = d->ran->iRandom(0, sumEmig-1);                  // random number for choosing neighbor
                        sumEmig = 0;
                        if (sumEmig == 0) {
                            neighborGroup = 0; break;
                        }
                        for (int i = 0; i < numNeighbors; i++) {
                            neighborGroup = (HaystackGroup*)(d->groupData) + neighbors[i];  // point to neighbor group
                            sumEmig += neighborGroup->emigrationPotential;
                            if (sumEmig >= choose) {
                                iNeighbor = i;
                                break;
                            }
                        }
                        break;}
                    }

                    int32_t numMigrants = 0;         // number of immigrants or colonizers

                    if (group->age == 0) {  // group is extinct. recolonize
                        //numMigrants = d->colonySize;
                        numMigrants = d->ran->poisson(d->colonySize);
                    }
                    else {  // group is not extinct. calculate number of immigrants
                        switch (d->immigrationPattern) {
                        case immigrationCommonPool:      // from common pool
                            numMigrants = d->ran->poisson(migrationRate * group->nmax);
                            break;
                        case immigrationNeighbor:        // from random neighbor group
                        case immigrationRandomGroup:     // from random group
                        case immigrationAllNeighbors:    // from all neighbors
                            if (neighborGroup) {
                                numMigrants = d->ran->poisson(migrationRate * neighborGroup->emigrationPotential);
                            }
                            break;
                        case immigrationStrongNeighbor:  // from neighbor selected according to emigrationPotential
                            if (neighborGroup) {
                                numMigrants = d->ran->poisson(migrationRate * neighborGroup->nn);
                            }
                            break;
                        case immigrationProportional:    // prop. w. own population, from neighbor
                            numMigrants = d->ran->poisson(migrationRate * group->nn);
                            break;
                        case immigrationVacant:          // prop. w. vacant capacity, from neighbor
                            numMigrants = group->nmax - group->nn + 1;
                            if (numMigrants < 0) numMigrants = 0;
                            else numMigrants = d->ran->poisson(migrationRate * numMigrants);
                            break;
                        }
                    }

                    // get migrant genes
                    if (neighborGroup && numMigrants > 0) {
                        if (numMigrants > neighborGroup->nn) numMigrants = neighborGroup->nn;
                        // select migrant genes
                        int32_t gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborGroup->gAltruism[0], neighborGroup->nn * 2);
                        int32_t gmig1 = numMigrants * 2 - gmig0;
                        group->gAltruism[0] += gmig0;           // add migrant genes to current group
                        group->gAltruism[1] += gmig1;
                        group->nn += numMigrants;
                        neighborGroup->gAltruism[0] -= gmig0;   // subtract migrant genes from neighbor group
                        neighborGroup->gAltruism[1] -= gmig1;
                        neighborGroup->nn -= numMigrants;
                        d->migrantsTot += numMigrants;         // count total number of migrants
                    }
                }

                // limit population
                if (group->nn > group->nmax) {
                    int32_t nn2 = group->nmax * 2;          // number of genes after population limitation
                    int32_t g0 = d->ran->hypergeometric(nn2, group->gAltruism[0], group->nn * 2);
                    int32_t g1 = nn2 - g0;
                    group->gAltruism[0] = g0;
                    group->gAltruism[1] = g1;
                    group->nn = group->nmax;
                }
                group->age++;                               // update age

                errorCheck(d, group);
            }

            // group loop 3: update global gene pool and statistics
            statisticsInit1(d);

            for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
                group = (HaystackGroup*)(d->groupData) + iGroup;// point to current group
                d->genePool[locusAltruism][0] += group->gAltruism[0];
                d->genePool[locusAltruism][1] += group->gAltruism[1];
                d->totalPopulation += group->nn;
                d->totalPhenotypes[locusAltruism] += group->altruists;
                if (group->groupfit > 0.5f) d->altruistGroups++;
                if (group->nn >= d->minGroupSize) d->inhabitedGroups++;
            }

            // make migrant pool
            if (d->immigrationPattern == immigrationCommonPool || d->colonizationPattern == colonizationCommonPool) {
                double poolGenes[2] = { 0, 0 };
                for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {  // group loop
                    // make contribution to migrant pool proportional to emigration potential
                    group = (HaystackGroup*)(d->groupData) + iGroup;// point to current group
                    poolGenes[0] += double(group->gAltruism[0]) * group->emigrationPotential;
                    poolGenes[1] += double(group->gAltruism[1]) * group->emigrationPotential;
                }
                double poolNn = (poolGenes[0] + poolGenes[1]) * 0.5;
                double reduction = 1.;
                if (poolNn > 1000000.) reduction = 1000000. / poolNn;

                // migrant pool placed last in memory block
                HaystackGroup * migrantPool = (HaystackGroup*)(d->groupData) + d->maxIslands;
                migrantPool->gAltruism[0] = (int)lround(poolGenes[0] * reduction);
                migrantPool->gAltruism[1] = (int)lround(poolGenes[1] * reduction);
                migrantPool->nn = (migrantPool->gAltruism[0] + migrantPool->gAltruism[1]) / 2u;
                migrantPool->emigrationPotential = d->nMaxPerGroup;
                migrantPool->nmax = migrantPool->nn;
                migrantPool->groupfit = 0.5f;
            }

            if (++*haystackGenerations >= d->haystackPeriod) {
                // finished haystack state. switch to mixing state
                *haystackState = stateMixing;
                *mixGenerations = 0;
            }
        }
        // update statistics counters
        statisticsUpdate(d);
        d->generations++;  // count generations

        // check if finished
        checkStopCriterion(d);
    }
}

