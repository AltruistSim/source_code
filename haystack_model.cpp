/****************************  haystack_model.cpp   *****************************
* Author:        Agner Fog
* Date created:  2023-09-04
* Last modified: 2023-12-31
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the island model with fixed geographic boundaries between groups
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

// Deme structure, describing each deme or territory
struct HaystackDeme {
    int32_t nmax;            // carrying capacity or max deme size
    int32_t nn;              // number of individuals in deme
    int32_t gAltruism[2];    // gene pool, egoism, altruism
    int32_t age;             // group age
    int32_t altruists;       // number of phenotypic altruists
    int32_t emigrationPotential; // excess available for emigration
    float   groupfit;        // group fitness
};

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
static const DemeFieldDescriptor haystackDemeDescriptors[] = {
    {1, varInt32, sizeof(HaystackDeme), 0, 0, 0},
    {2, varInt32, demeFieldOffset(HaystackDeme,nmax),           0, 2,  "carrying capacity"},
    {3, varInt32, demeFieldOffset(HaystackDeme,nn),             0, 3,  "population"},
    {4, varInt32, demeFieldOffset(HaystackDeme,gAltruism[0]),   8, 0,  "egoism gene"},
    {4, varInt32, demeFieldOffset(HaystackDeme,gAltruism[1]),   5, 10, "altruism gene"},
    {8, varInt32, demeFieldOffset(HaystackDeme,age),            0,  0, "age"},
    {8, varFloat, demeFieldOffset(HaystackDeme,groupfit),       2, 20, "fraction of altruists"},
    {8, varInt32, demeFieldOffset(TerriDeme,emigrationPotential), 0, 20, "emigration potential"},    
    {0, 0, 0, 0, 0, 0}   // mark end of list
};


static const ModelDescriptor haystackModel = {
    "Haystack model",
    "Intrademic group selection."
    " Haystacks are colonized periodically. Population is mixed when haystacks are removed. "
    "Atruist have lower individual fitness, but increase fitness of all members of the same haystack.",
    haystackParameterDefinitions,
    haystackDemeDescriptors,
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
        d->modelDemeStructureSize = sizeof(HaystackDeme);
        d->modelPopOffset = demeFieldOffset(HaystackDeme,nn); // make nn available to findNeighbors function
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

static void errorCheck(AltruData * d, HaystackDeme * deme) {
    // check deme data for consistency
        char text[64];
    if (deme->gAltruism[0] < 0 || deme->gAltruism[1] < 0 
        || deme->gAltruism[0] + deme->gAltruism[1] != deme->nn * 2 
        || (deme->gAltruism[0] + deme->gAltruism[1] & 1)) {
        sprintf_s(text, "error n %i, gAltru %i %i",
            deme->nn, deme->gAltruism[0], deme->gAltruism[1]);
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
    int32_t ideme;                                         // deme index
    HaystackDeme * deme;                                   // point to deme
    int32_t genotypes1[3];                                 // egoists, heterozygote, altruists before selection
    int32_t genotypes2[3];                                 // egoists, heterozygote, altruists after selection
    double fitness[3];                                     // fitness of each genotype
    int32_t neighbors[8];                                  // list of neighbor islands
    int numNeighbors;                                      // number of neighbor islands
    int32_t iNeighbor;                                     // index to neighbor islands
    int strongestNeighbor;                                 // index to strongest neighbor deme
    HaystackDeme * neighborDeme = 0;                       // pointer to neighbor deme
    double survivalRate[3] = {0};                          // survival rate for each genotype under independent viability selection only
    int32_t populationSize;                                // total population size


    const int inhabitedIslandsLength = 4096;
    int inhabitedIslandsNum;                               // number of inhabited islands in inhabitedIslands list
    // (may differ from d->inhabitedDemes if inhabitedIslands list is full)
    int32_t inhabitedIslands[inhabitedIslandsLength];      // used for list of inhabited islands

    if (d->demeData == 0) return;

    int32_t * const haystackState = &d->userData[0].i;     // 0: mixing period, 1: haystack period
    int32_t * const mixGenerations = &d->userData[1].i;    // number of generations currently spent in current mixing period 
    int32_t * const haystackGenerations = &d->userData[2].i;// number of generations currently spent in current haystack period 
    int32_t * const maxMetaPopulation = &d->modelspec_i[0];// carrying capacity in mixing state
    float   * const mixPeriodGrowthRate = &d->modelspec_f[1]; // growth rate during mixing period

    if (state == state_start) {
        // reset all demes
        memset(d->demeData, 0, (d->maxIslands+1) * sizeof(HaystackDeme));

        // populate deme 0 as mixing population
        deme = (HaystackDeme*)(d->demeData) + 0;           // point to deme 0
        populationSize = d->colonySize * d->maxIslands;
        deme->nmax = d->nMaxPerDeme;
        if (populationSize > 1000) populationSize = 1000;
        deme->nn = populationSize;
        deme->gAltruism[1] = d->ran->binomial(deme->nn * 2, d->fg0[locusAltruism]); // number of altruism genes
        deme->gAltruism[0] = deme->nn * 2 - deme->gAltruism[1]; // number of egoism genes 
        deme->age = 1;                                     // don't set age = 0, that means extinction
        deme->emigrationPotential = 1;
        deme->groupfit = 0.f;                              // not calculated yet
        errorCheck(d, deme);

        // initialize other variables
        int i;
        statisticsInit0(d);
        d->totalPopulation = populationSize;
        d->inhabitedDemes = 1;
        *haystackState = stateMixing;
        *mixGenerations = 1;
        d->nIslands = 1;

        // reset migrant pool
        deme = (HaystackDeme*)(d->demeData) + d->maxIslands;
        deme->nmax = deme->nn = 0;
        deme->gAltruism[0] = deme->gAltruism[1] = 0;
        deme->age = 0; 
        deme->groupfit = 0.f;
        deme->emigrationPotential = 0;

        // start running
        d->runState = state_run;
    }

    // simulate one generation 
    if (state == state_run) {
        if (*haystackState == stateMixing) {
            d->nIslands = 1;
            deme = (HaystackDeme*)(d->demeData) + 0;       // point to deme 0
            if (*mixGenerations == 0) {
                // first generation in mixing state
                // populate deme 0 as mixing population
                float reduction = 1.f;
                if (d->totalPopulation > 1000000) {        // limit population size to 1000000
                    reduction = 1000000.f / d->totalPopulation;
                    deme->gAltruism[0] = d->genePool[locusAltruism][0] * reduction;
                    deme->gAltruism[1] = d->genePool[locusAltruism][1] * reduction;
                    deme->nn = (deme->gAltruism[0] + deme->gAltruism[1]) / 2;
                }
                else {
                    deme->gAltruism[0] = d->genePool[locusAltruism][0];
                    deme->gAltruism[1] = d->genePool[locusAltruism][1];
                    deme->nn = d->totalPopulation;
                }
            }

            // mutation
            int32_t forwardMutations = d->ran->binomial(deme->gAltruism[0], d->murate[locusAltruism][0]);
            int32_t backwardMutations = d->ran->binomial(deme->gAltruism[1], d->murate[locusAltruism][1]);
            deme->gAltruism[1] += forwardMutations - backwardMutations;
            deme->gAltruism[0] -= forwardMutations - backwardMutations;
            d->mutations[locusAltruism][0] += forwardMutations;
            d->mutations[locusAltruism][1] += backwardMutations;

            // growth in mixing state without selection, with drift
            if (deme->nn > 1) {
                populationSize = d->ran->poisson(double(deme->nn) * *mixPeriodGrowthRate);
                if (populationSize < d->minGroupSize) populationSize = d->minGroupSize;
                int32_t me = d->ran->binomial(populationSize * 2, deme->gAltruism[0] / (deme->nn * 2.));
                int32_t ma = populationSize * 2 - me;
                deme->gAltruism[0] = me;
                deme->gAltruism[1] = ma;
                deme->nn = populationSize;

                // limit population size, with drift
                if (populationSize > *maxMetaPopulation) {
                    me = d->ran->hypergeometric(*maxMetaPopulation * 2, deme->gAltruism[0], populationSize * 2);
                    ma = *maxMetaPopulation * 2 - me;
                    deme->gAltruism[0] = me;
                    deme->gAltruism[1] = ma;
                    deme->nn = *maxMetaPopulation;
                }
            }
            deme->age = *mixGenerations;

            // statistics
            statisticsInit1(d);
            d->inhabitedDemes = d->altruistDemes = 0;
            d->totalPopulation = deme->nn;
            d->genePool[locusAltruism][0] = deme->gAltruism[0];
            d->genePool[locusAltruism][1] = deme->gAltruism[1];
            d->totalPhenotypes[locusAltruism] = deme->altruists;
            deme->groupfit = 0.f;
            if (deme->nn > 0) {
                d->inhabitedDemes = 1;
                deme->groupfit = deme->gAltruism[1] / (deme->nn * 2.0f);
                if (deme->groupfit > 0.5f) d->altruistDemes = 1; //else d->egoistDemes = 1;
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
                deme = (HaystackDeme*)(d->demeData) + 0;   // mixing population
                ge = deme->gAltruism[0];  ga = deme->gAltruism[1];
                nn = deme->nn;
                // deme loop 0: populate demes
                for (ideme = 0; ideme < d->maxIslands; ideme++) {
                    deme = (HaystackDeme*)(d->demeData) + ideme; // point to current deme
                    deme->nn = d->colonySize;
                    deme->nmax = d->nMaxPerDeme;
                    int32_t me = d->ran->binomial(deme->nn * 2, ge / (nn * 2.));
                    int32_t ma = deme->nn * 2 - me;
                    deme->gAltruism[0] = me;
                    deme->gAltruism[1] = ma;
                    deme->age = 1;
                    deme->emigrationPotential = 0;
                    deme->groupfit = 0.f;
                }
            }

            // reset statictics variables
            inhabitedIslandsNum = 0;  d->demesDied = 0;
            d->mutations[0][0] = d->mutations[0][1] = 0;

            // deme loop 1: mutation, growth, selection, extinction
            for (ideme = 0; ideme < d->maxIslands; ideme++) {
                deme = (HaystackDeme*)(d->demeData) + ideme; // point to current deme

                errorCheck(d, deme);

                // mutation
                int32_t forwardMutations = d->ran->binomial(deme->gAltruism[0], d->murate[locusAltruism][0]);
                int32_t backwardMutations = d->ran->binomial(deme->gAltruism[1], d->murate[locusAltruism][1]);
                deme->gAltruism[1] += forwardMutations - backwardMutations;
                deme->gAltruism[0] -= forwardMutations - backwardMutations;
                d->mutations[locusAltruism][0] += forwardMutations;
                d->mutations[locusAltruism][1] += backwardMutations;

                // reproduction and growth
                if (d->selectionModel != selectionFecundity) {  // same growth rate for all except under fecundity selection
                    int32_t nn2 = d->ran->poisson(double(deme->nn) * d->growthRate);    // new population size
                    int32_t me = d->ran->binomial(nn2 * 2, deme->gAltruism[0] / (deme->nn * 2.));
                    int32_t ma = nn2 * 2 - me;
                    deme->gAltruism[0] = me;
                    deme->gAltruism[1] = ma;
                }
                deme->nn = (deme->gAltruism[0] + deme->gAltruism[1]) / 2u;

                // calculate fraction of altruists
                combineGenes(deme->gAltruism, genotypes1, d->ran);
                deme->altruists = genotypes1[2];      // number of phenotypic altruists
                switch (d->dominance[locusAltruism]) {
                case recessive:                       // altruism recessive
                    break;
                case dominant:                        // altruism dominant
                    deme->altruists += genotypes1[1];
                    break;
                case incompleteDominant:               // half dominant
                    deme->altruists += genotypes1[1] / 2;
                }
                float fractionOfAltruists = 0.f;
                if (deme->nn) {                       // avoid division by 0          
                    fractionOfAltruists = float(deme->altruists) / float(deme->nn);
                }

                // calculate group fitness from fraction of altruists and group fitness exponent
                switch (d->fitfunc) {
                case 0:                               // one altruist is enough
                    if (deme->altruists > 0) deme->groupfit = 1.f;
                    else deme->groupfit = 0.f;
                    break;
                case 1: case 3:                       // convex/concave
                    deme->groupfit = powf(fractionOfAltruists, d->fitExpo);
                    break;
                case 2:  default:                     // linear
                    deme->groupfit = fractionOfAltruists;
                    break;
                case 4:                               // all or nothing
                    if (deme->altruists == deme->nn) deme->groupfit = 1.f;
                    else deme->groupfit = 0.f;
                    break;
                }

                // calculate relative fitness for egoists and altruists
                fitness[0] = d->fit[0] * (1.f - deme->groupfit) + d->fit[1] * deme->groupfit; // fitness of egoists
                fitness[2] = d->fit[2] * (1.f - deme->groupfit) + d->fit[3] * deme->groupfit; // fitness of altruists
                switch (d->dominance[0]) {  // find fitness of heterozygotes
                case 0: // altruism recessive
                    fitness[1] = fitness[0];  break;
                case 1: // altruism dominant
                    fitness[1] = fitness[2];  break;
                case 2: // half dominant
                    fitness[1] = 0.5 * (fitness[0] + fitness[2]);  break;
                }

                int32_t populationSize = deme->nn;         // new population size
                if (populationSize > deme->nmax) populationSize = deme->nmax;
                int32_t populationExcess = 0;              // excess population

                // selection according to selection model
                switch (d->selectionModel) {
                case selectionFecundity: default:   // fecundity selection
                    differentialGrowth(deme->gAltruism, fitness, genotypes1, d->ran);
                    for (int i = 0; i < 3; i++) {
                        genotypes2[i] = genotypes1[i];
                    }
                    populationSize = genotypes1[0] + genotypes1[1] + genotypes1[2];
                    if (populationSize > deme->nmax) {
                        populationExcess = populationSize - deme->nmax;
                        if (populationExcess < 0) populationExcess = 0;
                        //populationSize = deme->nmax;
                    }
                    // wait with reducing the population size until after migration
                    deme->emigrationPotential = populationExcess;
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
                    deme->emigrationPotential = d->nMaxPerDeme / 2u;
                    break;
                case emigrationProportional:                       // emigration proportional to population
                    deme->emigrationPotential = deme->nn;
                    break;
                case emigrationExcess:                             // emigration proportional to excess population
                    if (d->selectionModel != selectionFecundity) {
                        // (emigrationPotential has been calculated above if selectionModel == selectionFecundity)
                        for (int i = 0; i < 3; i++) {
                            populationExcess += (int32_t)lround(genotypes1[i] * fitness[i]);
                        }
                        populationExcess -= deme->nmax;
                        if (populationExcess < 0) populationExcess = 0;
                        deme->emigrationPotential = populationExcess;
                    }
                    break;
                case emigrationGroupfit:
                    deme->emigrationPotential = (int32_t)lround(deme->nn * deme->groupfit);
                    break;
                }

                // drift in offspring of heterozygotes
                int32_t g0 = d->ran->binomial(genotypes2[1]*2, 0.5);
                int32_t g1 = genotypes2[1]*2 - g0;

                // genes after selection
                deme->gAltruism[0] = genotypes2[0] * 2 + g0;
                deme->gAltruism[1] = genotypes2[2] * 2 + g1;
                deme->nn = (deme->gAltruism[0] + deme->gAltruism[1]) / 2;  // number of individuals
                // calculate fraction of altruists after selection
                deme->altruists = genotypes2[2]; // number of phenotypic altruists
                switch (d->dominance[locusAltruism]) {
                case 0: // altruism recessive
                    break;
                case 1: // altruism dominant
                    deme->altruists += genotypes2[1];
                    break;
                case 2: // half dominant
                    deme->altruists += genotypes2[1] / 2;
                }

                //errorCheck(d, deme);

                // extinction
                float relativeSize = float(deme->nn) / float(deme->nmax); // relative size of deme
                if (relativeSize > 1.f) relativeSize = 1.0f;

                float extinctionRate =
                    (1.f - deme->groupfit) * ((1.f - relativeSize) * d->extinctionRate[0] + relativeSize * d->extinctionRate[1])
                    + deme->groupfit * ((1.f - relativeSize) * d->extinctionRate[2] + relativeSize * d->extinctionRate[3]);

                // kill with probability extinctionRate
                bool kill = d->ran->bernoulli(extinctionRate);

                if (kill) {  // this deme goes extinct
                    d->demesDied++;
                    deme->emigrationPotential = 0;
                    deme->age = 0;
                    if (d->surviv > 0.f) {  // a fraction survives
                        int32_t numSurvivors = int32_t(d->surviv * deme->nn);
                        deme->gAltruism[0] = d->ran->hypergeometric(numSurvivors * 2, deme->gAltruism[0], deme->nn * 2);
                        deme->gAltruism[1] = numSurvivors * 2 - deme->gAltruism[0];
                        deme->nn = numSurvivors;
                        deme->groupfit = 0.f;
                    }
                    else {   // all die
                        deme->gAltruism[0] = deme->gAltruism[1] = 0;
                        deme->nn = 0;
                        deme->groupfit = 0.f;
                    }
                }

                // make list of inhabited islands/haystacks 
                if (deme->nn > 0 && inhabitedIslandsNum < inhabitedIslandsLength) {
                    inhabitedIslands[inhabitedIslandsNum++] = ideme;
                }
            }

            // deme loop 2: migration or colonization
            for (ideme = 0; ideme < d->maxIslands; ideme++) {
                deme = (HaystackDeme*)(d->demeData) + ideme;     // point to current deme

                //errorCheck(d, deme);

                // make list of neighbors
                numNeighbors = findNeighbors(d, ideme, neighbors);
                neighborDeme = 0;
                float migrationRate = d->migrationRate[0];
                HaystackDeme combinedNeighbors;
                int migrationOrColonizationPattern;  // immigration and colonization use similar patterns
                if (deme->age == 0 || deme->nn == 0) {  // group is extinct or empty. recolonize
                    deme->age = 0;                      // reset age if group is empty but not extinct
                    migrationOrColonizationPattern = d->colonizationPattern;
                }
                else {   // group is not extinct. calculate immigration
                    migrationOrColonizationPattern = d->immigrationPattern;
                }

                neighborDeme = 0;
                if (numNeighbors) {
                    // immigration or colonization
                    switch (migrationOrColonizationPattern) {
                    case immigrationCommonPool:
                        // migrants come from common pool
                        neighborDeme = (HaystackDeme*)(d->demeData) + d->maxIslands; // point to common pool
                        break;
                    case immigrationNeighbor:
                    case immigrationProportional:                        
                    case immigrationVacant: 
                        // migrants come from random neighbor
                        iNeighbor = d->ran->iRandom(0, numNeighbors - 1);   // find random neighbor deme
                        neighborDeme = (HaystackDeme*)(d->demeData) + neighbors[iNeighbor];  // point to neighbor deme
                        break;
                    case immigrationAllNeighbors: {// sum of all neighbors
                        memset(&combinedNeighbors, 0, sizeof(combinedNeighbors));  // set to all 0
                        double sumGene0 = 0, sumGene1 = 0, sumEmigr = 0;
                        for (int i = 0; i < numNeighbors; i++) {
                            neighborDeme = (HaystackDeme*)(d->demeData) + neighbors[i];  // point to neighbor deme
                            sumGene0 += double(neighborDeme->gAltruism[0]) * neighborDeme->emigrationPotential;
                            sumGene1 += double(neighborDeme->gAltruism[1]) * neighborDeme->emigrationPotential;
                            sumEmigr += neighborDeme->emigrationPotential;
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
                        neighborDeme = &combinedNeighbors;
                        break;}
                    case immigrationRandomGroup:
                        if (inhabitedIslandsNum > 0) {  // immigration from random inhabited island
                            int i = d->ran->iRandom(0, inhabitedIslandsNum - 1);         // find random inhabited island
                            iNeighbor = inhabitedIslands[i];
                            neighborDeme = (HaystackDeme*)(d->demeData) + iNeighbor;     // point to random deme
                        }
                        break;
                    case immigrationStrongNeighbor: {
                        // choose among all neighbors with probability proportional to emigrationPotential
                        int32_t sumEmig = 0;  // sum of emigration potential of all neighbors
                        for (int i = 0; i < numNeighbors; i++) {
                            neighborDeme = (HaystackDeme*)(d->demeData) + neighbors[i];  // point to neighbor deme
                            sumEmig += neighborDeme->emigrationPotential;
                        }
                        int32_t choose = d->ran->iRandom(0, sumEmig-1);                  // random number for choosing neighbor
                        sumEmig = 0;
                        if (sumEmig == 0) {
                            neighborDeme = 0; break;
                        }
                        for (int i = 0; i < numNeighbors; i++) {
                            neighborDeme = (HaystackDeme*)(d->demeData) + neighbors[i];  // point to neighbor deme
                            sumEmig += neighborDeme->emigrationPotential;
                            if (sumEmig >= choose) {
                                iNeighbor = i;
                                break;
                            }
                        }
                        break;}
                    }

                    int32_t numMigrants = 0;         // number of immigrants or colonizers

                    if (deme->age == 0) {  // group is extinct. recolonize
                        //numMigrants = d->colonySize;
                        numMigrants = d->ran->poisson(d->colonySize);
                    }
                    else {  // group is not extinct. calculate number of immigrants
                        switch (d->immigrationPattern) {
                        case immigrationCommonPool:      // from common pool
                            numMigrants = d->ran->poisson(migrationRate * deme->nmax);
                            break;
                        case immigrationNeighbor:        // from random neighbor deme
                        case immigrationRandomGroup:     // from random deme
                        case immigrationAllNeighbors:    // from all neighbors
                        case immigrationStrongNeighbor:  // from neighbor selected according to emigrationPotential
                            if (neighborDeme) {
                                numMigrants = d->ran->poisson(migrationRate * neighborDeme->emigrationPotential);
                            }
                            break;
                        case immigrationProportional:    // prop. w. population, from neighbor
                            numMigrants = d->ran->poisson(migrationRate * deme->nn);
                            break;
                        case immigrationVacant:          // prop. w. vacant capacity, from neighbor
                            numMigrants = deme->nmax - deme->nn + 1;
                            if (numMigrants < 0) numMigrants = 0;
                            else numMigrants = d->ran->poisson(migrationRate * numMigrants);
                            break;
                        }
                    }

                    // get migrant genes
                    if (neighborDeme && numMigrants > 0) {
                        if (numMigrants > neighborDeme->nn) numMigrants = neighborDeme->nn;
                        // select migrant genes
                        int32_t gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborDeme->gAltruism[0], neighborDeme->nn * 2);
                        int32_t gmig1 = numMigrants * 2 - gmig0;
                        deme->gAltruism[0] += gmig0;           // add migrant genes to current deme
                        deme->gAltruism[1] += gmig1;
                        deme->nn += numMigrants;
                        neighborDeme->gAltruism[0] -= gmig0;   // subtract migrant genes from neighbor deme
                        neighborDeme->gAltruism[1] -= gmig1;
                        neighborDeme->nn -= numMigrants;
                        d->migrantsTot += numMigrants;         // count total number of migrants
                    }
                }

                // limit population
                if (deme->nn > deme->nmax) {
                    int32_t nn2 = deme->nmax * 2;          // number of genes after population limitation
                    int32_t g0 = d->ran->hypergeometric(nn2, deme->gAltruism[0], deme->nn * 2);
                    int32_t g1 = nn2 - g0;
                    deme->gAltruism[0] = g0;
                    deme->gAltruism[1] = g1;
                    deme->nn = deme->nmax;
                }
                deme->age++;                               // update age

                errorCheck(d, deme);
            }

            // deme loop 3: update global gene pool and statistics
            statisticsInit1(d);

            for (ideme = 0; ideme < d->maxIslands; ideme++) {
                deme = (HaystackDeme*)(d->demeData) + ideme;// point to current deme
                d->genePool[locusAltruism][0] += deme->gAltruism[0];
                d->genePool[locusAltruism][1] += deme->gAltruism[1];
                d->totalPopulation += deme->nn;
                d->totalPhenotypes[locusAltruism] += deme->altruists;
                if (deme->groupfit > 0.5f) d->altruistDemes++; //else d->egoistDemes++;
                if (deme->nn >= d->minGroupSize) d->inhabitedDemes++;
            }

            // make migrant pool
            if (d->immigrationPattern == immigrationCommonPool || d->colonizationPattern == colonizationCommonPool) {
                double poolGenes[2] = { 0, 0 };
                for (ideme = 0; ideme < d->maxIslands; ideme++) {  // deme loop
                    // make contribution to migrant pool proportional to emigration potential
                    deme = (HaystackDeme*)(d->demeData) + ideme;// point to current deme
                    poolGenes[0] += double(deme->gAltruism[0]) * deme->emigrationPotential;
                    poolGenes[1] += double(deme->gAltruism[1]) * deme->emigrationPotential;
                }
                double poolNn = (poolGenes[0] + poolGenes[1]) * 0.5;
                double reduction = 1.;
                if (poolNn > 1000000.) reduction = 1000000. / poolNn;

                // migrant pool placed last in memory block
                HaystackDeme * migrantPool = (HaystackDeme*)(d->demeData) + d->maxIslands;
                migrantPool->gAltruism[0] = (int)lround(poolGenes[0] * reduction);
                migrantPool->gAltruism[1] = (int)lround(poolGenes[1] * reduction);
                migrantPool->nn = (migrantPool->gAltruism[0] + migrantPool->gAltruism[1]) / 2u;
                migrantPool->emigrationPotential = d->nMaxPerDeme;
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

