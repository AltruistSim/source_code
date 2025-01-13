/****************************  island_model.cpp   *****************************
* Author:        Agner Fog
* Date created:  1994-09-25
* Last modified: 2025-01-09
* Version:       3.003
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the island model with fixed geographic boundaries between groups
*
* (c) Copyright 2023-2025 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

// Comment out this to save cache space if conformity and endogamy loci not used.
// This improves speed by approximately 5%
#define MORELOCI

void islandInitFunction(AltruData * d, int index);
void islandGenerationFunction(AltruData * d, int mode);


/******************************************************************************
*         structures and definitions
******************************************************************************/

// List of model-specific parameters
// These parameters will appear in parameter list and dialog box.
// Maximum number of model-specific parameters is maxModelSpecificParameters defined in altruist.h
// ParameterDef::type can be: 1: integer, 2: float, 3: boolean
// ParameterDef::num must be 1 if not an array
// Offsets and names must be unique
static const ParameterDef islandParameterDefinitions[] {
    {0, 1, 0, "Model-specific parameters"},
    {0, 0, 0, 0}
};

// group structure, describing each group or territory
struct IslandGroup {
    int32_t nn;              // number of individuals in group
    int32_t nmax;            // carrying capacity or max group size
    int32_t gAltruism[2];    // gene pool, egoism, altruism
#ifdef MORELOCI
    int32_t gEndogamy[2];    // gene pool, exogamy, endogamy
    int32_t gConformity[2];  // gene pool, nonconformity, conformity
#endif
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
static const GroupFieldDescriptor islandGroupDescriptors[] = {
    {1, varInt32, sizeof(IslandGroup), 0, 0, 0},
    {3, varInt32, groupFieldOffset(IslandGroup,nn),             0, 3,  "population"},
    {2, varInt32, groupFieldOffset(IslandGroup,nmax),           0, 2,  "carrying capacity"},
    {4, varInt32, groupFieldOffset(IslandGroup,gAltruism[0]),   8, 0,  "egoism gene"},
    {4, varInt32, groupFieldOffset(IslandGroup,gAltruism[1]),   5, 10, "altruism gene"},
#ifdef MORELOCI
    {4, varInt32, groupFieldOffset(IslandGroup,gEndogamy[0]),   2, 0,  "exogamy gene"},
    {4, varInt32, groupFieldOffset(IslandGroup,gEndogamy[1]),   2, 11, "endogamy gene"},
    {4, varInt32, groupFieldOffset(IslandGroup,gConformity[0]), 2, 0,  "nonconformity gene"},
    {4, varInt32, groupFieldOffset(IslandGroup,gConformity[1]), 2, 12, "conformity gene"},
#endif
    {8, varInt32, groupFieldOffset(IslandGroup,age),            0, 20, "age"},
    {8, varInt32, groupFieldOffset(IslandGroup,altruists),      0, 20, "phenotypic altruists"},
    {8, varInt32, groupFieldOffset(IslandGroup,emigrationPotential), 0, 20, "emigration potential"},    
    {8, varFloat, groupFieldOffset(IslandGroup,groupfit),       2, 20, "group fitness"},
    {0, 0, 0, 0, 0, 0}       // mark end of list
};


static const ModelDescriptor islandModel = {
    "Island model",
    "Groups on isolated islands. Probability of extinction depends on group fitness."
    " Vacant islands are recolonized. "
    "Atruism will decrease individual fitness, but increase group fitness."
    "\nThis model can be extended with endogamy and conformity.",
    islandParameterDefinitions,
    islandGroupDescriptors,
    &islandInitFunction,
    &islandGenerationFunction
};

// This will register the Island model in the global model descriptor list
Construct islandConstructor(islandModel);


/******************************************************************************
*         initialization function
* state values:
* 1: model changed
* 2: parameters changed
******************************************************************************/

void islandInitFunction(AltruData * d, int state) {
    switch (state) {
    case model_initialize:  // initialization

        // model version
        d->modelVersionMajor = 3;
        d->modelVersionMinor = 0;
        d->modelGroupStructureSize = sizeof(IslandGroup);
        d->modelPopOffset = groupFieldOffset(IslandGroup,nn);// make nn available to findNeighbors function
        // enable and disable fields in parameters dialog boxes
        d->bGeographyParametersUsed = 0x351312;
        d->bIndividualPropertiesUsed = 0xF;
        if (d->locusUsed[locusEndogamy])   d->bIndividualPropertiesUsed |= 0x00100;
        if (d->locusUsed[locusConformity]) d->bIndividualPropertiesUsed |= 0x104FF;
        d->bGroupPropertiesUsed = 0x1076;
        d->bImmigrationPattern = 0b1111111;                     // immigration patterns supported
        d->bColonPat = 0b1110011;                     // colonization patterns supported
        d->bEmigrationPattern = 0b1111;
        d->bExtinctionPatterns = 0b0111;
        d->bTopology = 0b01111111;
        d->bWarPatterns = 0;
        d->bStopCriterionUsed = 7;
        d->nLoci = 3;                                      // number of loci
        if (d->locusUsed[locusEndogamy]) d->bGeographyParametersUsed |= 0x2000;  // migration under endogamy
        d->graphicsTypeForModel = graphicsIslands;
        for (int locus = 0; locus < d->nLoci; locus++) {
            if (d->locusUsed[locus]) {
                d->bStopCriterionUsed |= 3 << (locus*2 + 3);
            }
        }

        d->sLocusName[locusAltruism]          = "altruism";
        d->sGeneName[locusAltruism][0]        = "egoism";
        d->sGeneName[locusAltruism][1]        = "altruism";
        d->sPhenotypeName[locusAltruism][0]   = "egoist";
        d->sPhenotypeName[locusAltruism][1]   = "altruist";
        d->sLocusName[locusEndogamy]          = "endogamy";
        d->sGeneName[locusEndogamy][0]        = "exogamy";
        d->sGeneName[locusEndogamy][1]        = "endogamy";
        d->sPhenotypeName[locusEndogamy][0]   = "exogamist";
        d->sPhenotypeName[locusEndogamy][1]   = "endogamist";
        d->sLocusName[locusConformity]        = "conformity";
        d->sGeneName[locusConformity][0]      = "nonconformity";
        d->sGeneName[locusConformity][1]      = "conformity";
        d->sPhenotypeName[locusConformity][0] = "nonconformist";
        d->sPhenotypeName[locusConformity][1] = "conformist";

        d->fitRowLocus = locusAltruism;
        d->fitColLocus = locusAltruism;
        d->fitSupColLocus = locusConformity;
        d->fitConditionName = 1;
        break;

    case model_parameters_changed:  // parameters changed
        d->nIslands = d->maxIslands;
        if (d->locusUsed[locusEndogamy]) d->bGeographyParametersUsed |= 0x2000;  // migration under endogamy
        break;
    }
}

static void errorCheck(AltruData * d, IslandGroup * group) {
    // check group data for consistency
    char text[64];
    if (group->gAltruism[0] < 0 || group->gAltruism[1] < 0
        || group->gAltruism[0] + group->gAltruism[1] != group->nn * 2
        || (group->gAltruism[0] + group->gAltruism[1] & 1)) {
        sprintf_s(text, "error n %i, gAltru %i %i",
            group->nn, group->gAltruism[0], group->gAltruism[1]);
        errors.reportError(text);
    }
#ifdef MORELOCI
    if ((d->locusUsed[locusEndogamy] && group->gEndogamy[0] + group->gEndogamy[1] != group->nn * 2) ||
        (d->locusUsed[locusConformity] && group->gConformity[0] + group->gConformity[1] != group->nn * 2)) {
        sprintf_s(text, "error n %i, gEndo %i %i, gConf %i %i",
            group->nn, group->gEndogamy[0], group->gEndogamy[1], group->gConformity[0], group->gConformity[1]);
        errors.reportError(text);
    }
#endif
}


/******************************************************************************
            Generation function: simulate one generation

Discrete non-overlapping generations (Wright-Fisher model, He et al 2017):
mating, reproduction (with drift and possible fecundity selection), mutation, 
viability selection, population regulation, extinction, migration
******************************************************************************/

void islandGenerationFunction(AltruData * d, int state) {
    int32_t iGroup;                                        // group index
    IslandGroup * group;                                   // point to group
    int32_t genotypes1[3];                                 // egoists, heterozygote, altruists before selection
    int32_t genotypes2[3];                                 // egoists, heterozygote, altruists after selection
    double fitness[3];                                     // fitness of each genotype
    double survivalRate[3] = {0};                          // survival rate for each genotype under independent viability selection only

    int32_t neighbors[8];                                  // list of neighbor islands
    int numNeighbors;                                      // number of neighbor islands
    int32_t iNeighbor;                                     // index to neighbor islands
    IslandGroup * neighborGroup = 0;                       // pointer to neighbor group
    int inhabitedIslandsNum;                               // number of inhabited islands in inhabitedIslands list

    // check allocated memory
    if (d->groupData == 0) return;
    if (d->extraBufferSize[0] < d->maxIslands * sizeof(int32_t)) {
        // allocate memory for list of inhabited islands
        if (d->extraBuffer[0]) delete[] d->extraBuffer[0]; // delete any previous smaller buffer
        d->extraBuffer[0] = new int8_t[d->maxIslands * sizeof(int32_t)]();
        if (d->extraBuffer[0] == 0) errors.reportError("memory allocation failed");
        d->extraBufferSize[0] = d->maxIslands * sizeof(int32_t);
    }
    int32_t * inhabitedIslands = (int32_t*)d->extraBuffer[0];


    if (state == state_start) {

        // initialize statistics variables
        statisticsInit0(d);

        // initialize groups
        d->nIslands = d->maxIslands;

        for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
            group = (IslandGroup*)(d->groupData) + iGroup;     // point to current group
            //group->nmax = d->nMaxPerGroup;                 // max population
            // randomize island size
            int32_t nmax;
            if (d->carryingCapacityStandardDeviation > 0.) {
                // randomize carrying capacity
                double upperLimit = d->carryingCapacity[0] + d->carryingCapacityStandardDeviation * 5.;
                if (upperLimit > 16000.) upperLimit = d->carryingCapacity[0] + 10000.;
                nmax = (int)lround(d->ran->normalTrunc(d->carryingCapacity[0], d->carryingCapacityStandardDeviation, 0., upperLimit));
                if (nmax < d->minGroupSize) nmax = d->minGroupSize;
            }
            else {
                // constant carrying capacity
                nmax = (int)lround(d->carryingCapacity[0]);
            }
            group->nmax = nmax;
            d->nMaxPerGroup = nmax;

            group->nn = d->ran->poisson(group->nmax / 2u);        // random population size
            if (group->nn > group->nmax) group->nn = group->nmax;   // limit population
            if (group->nn < d->minGroupSize) group->nn = 0;       // extinction if below minimum
            group->gAltruism[1] = d->ran->binomial(group->nn * 2, d->fg0[locusAltruism]); // number of altruism genes
            group->gAltruism[0] = group->nn * 2 - group->gAltruism[1];                     // number of egoism genes 
#ifdef MORELOCI
            if (d->locusUsed[locusEndogamy]) {
                group->gEndogamy[1] = d->ran->binomial(group->nn * 2, d->fg0[locusEndogamy]); // number of endogamy genes
                group->gEndogamy[0] = group->nn * 2 - group->gEndogamy[1];                     // number of exogamy genes 
            }
            if (d->locusUsed[locusConformity]) {
                group->gConformity[1] = d->ran->binomial(group->nn * 2, d->fg0[locusConformity]); // number of conformity genes
                group->gConformity[0] = group->nn * 2 - group->gConformity[1];                     // number of nonconformity genes 
            }
#endif
            group->age = 1;                                 // don't set age = 0, that means extinction
            group->emigrationPotential = 1; // not calculated yet
            group->groupfit = 0.f;      // not calculated yet
            errorCheck(d, group);
        }
        // reset migrant pool
        group = (IslandGroup*)(d->groupData) + d->maxIslands;
        group->nmax = group->nn = 0;
        group->gAltruism[0] = group->gAltruism[1] = 0;
#ifdef MORELOCI
        group->gConformity[0] = group->gConformity[1] = 0;
        group->gEndogamy[0]   = group->gEndogamy[1]   = 0;
#endif
        group->age = 1; 
        group->groupfit = 0.f;
        group->emigrationPotential = 0;

        // start running
        d->runState = state_run;
    }

    // simulate one generation 
    if (state == state_run) {

        // reset statictics variables
        statisticsInit1(d);

        inhabitedIslandsNum = 0;

        // group loop 1: mutation, growth, selection, group selection
        for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
            group = (IslandGroup*)(d->groupData) + iGroup;     // point to current group

            errorCheck(d, group);

            // mutation
            int32_t forwardMutations = d->ran->binomial(group->gAltruism[0], d->murate[locusAltruism][0]);
            int32_t backwardMutations = d->ran->binomial(group->gAltruism[1], d->murate[locusAltruism][1]);
            group->gAltruism[1] += forwardMutations - backwardMutations;
            group->gAltruism[0] -= forwardMutations - backwardMutations;
            d->mutations[locusAltruism][0] += forwardMutations;
            d->mutations[locusAltruism][1] += backwardMutations;
#ifdef MORELOCI
            if (d->locusUsed[locusEndogamy]) {  // endogamy locus
                forwardMutations = d->ran->binomial(group->gEndogamy[0], d->murate[locusEndogamy][0]);
                backwardMutations = d->ran->binomial(group->gEndogamy[1], d->murate[locusEndogamy][1]);
                group->gEndogamy[1] += forwardMutations - backwardMutations;
                group->gEndogamy[0] -= forwardMutations - backwardMutations;
                d->mutations[locusEndogamy][0] += forwardMutations;
                d->mutations[locusEndogamy][1] += backwardMutations;
            }
            if (d->locusUsed[locusConformity]) {  // conformity locus
                forwardMutations = d->ran->binomial(group->gConformity[0], d->murate[locusConformity][0]);
                backwardMutations = d->ran->binomial(group->gConformity[1], d->murate[locusConformity][1]);
                group->gConformity[1] += forwardMutations - backwardMutations;
                group->gConformity[0] -= forwardMutations - backwardMutations;
                d->mutations[locusConformity][0] += forwardMutations;
                d->mutations[locusConformity][1] += backwardMutations;
            }
#endif

            // calculate fraction of altruists
            combineGenes(group->gAltruism, genotypes1, d->ran);
            group->altruists = genotypes1[2];     // number of phenotypic altruists
            switch (d->dominance[locusAltruism]) {
            case recessive:            // altruism recessive
                break;
            case dominant:             // altruism dominant
                group->altruists += genotypes1[1];
                break;
            case incompleteDominant:    // half dominant
                group->altruists += genotypes1[1] / 2;
            }
            float fractionOfAltruists = 0.f;
            if (group->nn) {            // avoid division by 0          
                fractionOfAltruists = float(group->altruists) / float(group->nn);
            }

            // calculate group fitness from fraction of altruists and group fitness exponent
            switch (d->fitfunc) {
            case 0:                    // one altruist is enough
                if (group->altruists > 0) group->groupfit = 1.f;
                else group->groupfit = 0.f;
                break;
            case 1: case 3:            // convex/concave
                group->groupfit = powf(fractionOfAltruists, d->groupFitCurvature);
                break;
            case 2:  default:          // linear
                group->groupfit = fractionOfAltruists;
                break;
            case 4:                    // all or nothing
                if (group->altruists == group->nn) group->groupfit = 1.f;
                else group->groupfit = 0.f;
                break;
            }

            // reproduction and growth
            if (d->selectionModel != selectionFecundity) {  // same growth rate for all except under fecundity selection
                int32_t nn2 = d->ran->poisson(double(group->nn) * d->growthRate);    // new population size
                int32_t me = d->ran->binomial(nn2 * 2, group->gAltruism[0] / (group->nn * 2.));
                int32_t ma = nn2 * 2 - me;
                group->gAltruism[0] = me;
                group->gAltruism[1] = ma;

                // split into genotypes
                combineGenes(group->gAltruism, genotypes1, d->ran);
#ifdef MORELOCI
                int32_t nGenes = group->gAltruism[0] + group->gAltruism[1];  // total number of genes
                if (d->locusUsed[locusEndogamy]) {
                    group->gEndogamy[0] = d->ran->binomial(nGenes, double(group->gEndogamy[0]) / double(group->gEndogamy[0] + group->gEndogamy[1]));
                    group->gEndogamy[1] = nGenes - group->gEndogamy[0];
                }
                if (d->locusUsed[locusConformity]) {
                    group->gConformity[0] = d->ran->binomial(nGenes, double(group->gConformity[0]) / double(group->gConformity[0] + group->gConformity[1]));
                    group->gConformity[1] = nGenes - group->gConformity[0];
                }
#endif
            }
            group->nn = (group->gAltruism[0] + group->gAltruism[1]) / 2u;
            group->emigrationPotential = group->nn - group->nmax;               // excess population can emigrate
            if (group->emigrationPotential < 0) group->emigrationPotential = 0;

            // find fitness of egoists and altruists
            float fit[4]; // fitness of egoists among egoists, egoists among altruists, altruists among egoists, altruists among altruists
            for (int i = 0; i < 4; i++) fit[i] = d->fit[i];
#ifdef MORELOCI
            if (d->locusUsed[locusConformity]) {
                // calculate level of conformity
                int32_t conformityGenotypes[3];
                // split gene pool into genotypes
                combineGenes(group->gConformity, conformityGenotypes, d->ran);
                // find number of phenotypic conformists
                int32_t conformists = conformityGenotypes[2];
                switch (d->dominance[locusConformity]) {
                case recessive:          // conformity recessive
                    break;
                case dominant:           // conformity dominant
                    conformists += conformityGenotypes[1];
                    break;
                case incompleteDominant: // half dominant
                    conformists += conformityGenotypes[1] / 2;
                }
                float conformity = 0.f;  // fraction of conformists
                if (group->nn > 0) {
                    conformity = float(conformists) / float(group->nn);
                }
                // adjust fitness coefficients to level of conformity
                for (int i = 0; i < 4; i++) {
                    fit[i] = d->fit[i] * (1.f - conformity) + d->fit[i+4] * conformity;
                }
            }
#endif
            // calculate relative fitness for egoists and altruists
            fitness[0] = fit[0] * (1.f - group->groupfit) + fit[1] * group->groupfit; // fitness of egoists
            fitness[2] = fit[2] * (1.f - group->groupfit) + fit[3] * group->groupfit; // fitness of altruists
            switch (d->dominance[locusAltruism]) {  // find fitness of heterozygotes
            case recessive: // altruism recessive
                fitness[1] = fitness[0];  break;
            case dominant: // altruism dominant
                fitness[1] = fitness[2];  break;
            case incompleteDominant: // half dominant
                fitness[1] = 0.5 * (fitness[0] + fitness[2]);  break;
            }

            int32_t populationSize = group->nn;
            if (populationSize > group->nmax) populationSize = group->nmax;
            
            // selection according to selection model
            switch (d->selectionModel) {

            case selectionFecundity: default:   // fecundity selection
                differentialGrowth(group->gAltruism, fitness, genotypes2, d->ran);
                populationSize = genotypes2[0] + genotypes2[1] + genotypes2[2];
                group->emigrationPotential = populationSize - group->nmax;
                if (group->emigrationPotential < 0) group->emigrationPotential = 0;
                if (populationSize > group->nmax) populationSize = group->nmax;
                // don't reduce the population size until after migration
                //d->ran->multiHypergeometric(genotypes2, genotypes1, populationSize, 3);
                break;

            case selectionPositiveViability:     // positive viability selection
                d->ran->multiWalleniusNCHyp(genotypes2, genotypes1, fitness, populationSize, 3);
                break;

            case selectionNegativeViability:     // negatitive viability selection
                d->ran->complMultiWalleniusNCHyp(genotypes2, genotypes1, fitness, populationSize, 3);
                break;

            case selectionMinimumViability:      // minimum viability selection                                
                d->ran->multiFishersNCHyp(genotypes2, genotypes1, fitness, populationSize, 3);
                break;

            case selectionIndependentViability: { // independent viability selection
                // Calculate survival rates for each phenotype under independent viability selection.
                // The population is kept constant by making survivalRate = 1 / growthRate
                // Population reduction takes place in two stages of approximately sqrt(survivalRate) each:
                // 1. independent viability selection, 
                // 2. proportional population reduction with no selection (except for genetic drift)
                double fitnessRatio = fitness[2] / fitness[0];
                survivalRate[0] = 1. / sqrt(d->growthRate * fitnessRatio);
                if (survivalRate[0] > 1.) survivalRate[0] = 1.;
                survivalRate[2] = survivalRate[0] * fitnessRatio;
                if (survivalRate[2] > 1.) {
                    survivalRate[0] /= survivalRate[2];
                    survivalRate[2] = 1.;
                }
                switch (d->dominance[locusAltruism]) {  // find survivalRate of heterozygotes
                case recessive: // altruism recessive
                    survivalRate[1] = survivalRate[0];  break;
                case dominant: // altruism dominant
                    survivalRate[1] = survivalRate[2];  break;
                case incompleteDominant: // half dominant
                    survivalRate[1] = 0.5 * (survivalRate[0] + survivalRate[2]);  break;
                }
                // select according to survival rate
                for (int i = 0; i < 3; i++) {
                    genotypes2[i] = d->ran->binomial(genotypes1[i], survivalRate[i]);
                }
                // wait with reducing population till after migration
                /*int32_t size2 = genotypes1[0] + genotypes1[1] + genotypes1[2];
                if (populationSize > size2) populationSize = size2;                
                d->ran->multiHypergeometric(genotypes2, genotypes1, populationSize, 3);*/
                break;}
            }

            // drift in offspring of heterozygotes
            int32_t g0 = d->ran->binomial(genotypes2[1]*2, 0.5);
            int32_t g1 = genotypes2[1]*2 - g0;

            // genes after selection
            group->gAltruism[0] = genotypes2[0] * 2 + g0;
            group->gAltruism[1] = genotypes2[2] * 2 + g1;
            group->nn = (group->gAltruism[0]+group->gAltruism[1]) / 2u; // number of individuals
#ifdef MORELOCI
            if (d->locusUsed[locusEndogamy]) {   // growth, selection, and drift at endogamy locus         
                growthAndSelection(group->gEndogamy, group->nn, d->dominance[locusEndogamy], d->fit2[0], d->ran);
            }
            if (d->locusUsed[locusConformity]) { // growth, selection, and drift at conformity locus         
                growthAndSelection(group->gConformity, group->nn, d->dominance[locusConformity], d->fit2[2], d->ran);
            }
#endif

            // emigration
            switch (d->emigrationPattern) {
            case emigrationConstant:                           // emigration independent of population
                group->emigrationPotential = d->nMaxPerGroup / 2u;
                break;
            case emigrationProportional:                       // emigration proportional to population
                group->emigrationPotential = group->nn;
                break;
            case emigrationExcess:                             // emigration proportional to excess population
                // emigrationPotential has been calculated above
                /*
                int32_t populationExcess = 0;
                for (int i = 0; i < 3; i++) {
                    populationExcess += (int32_t)lround(genotypes1[i] * fitness[i]);
                }
                populationExcess -= group->nmax;
                if (populationExcess < 0) populationExcess = 0;
                group->emigrationPotential = populationExcess;*/
                break;
            case emigrationGroupfit:
                group->emigrationPotential = (int32_t)lround(group->nn * group->groupfit);
                break;
            }

            // find neighbor groups if needed
            numNeighbors = 0;
            iNeighbor = 0;
            float fitness2 = 1.f;      // fitness relative to neighbor

            if (d->extinctionPattern == extinctionNeighbor || d->extinctionPattern == extinctionStrongest) {
                // make list of neighbors
                numNeighbors = findNeighbors(d, iGroup, neighbors);
            }

            switch (d->extinctionPattern) {
            case extinctionOwn:        // extinction probability depends on own group fitness only
                fitness2 = group->groupfit;
                break;

            case extinctionNeighbor:   // extinction probability depends on fitness relative to random neighbor
                if (numNeighbors) {
                    iNeighbor = d->ran->iRandom(0, numNeighbors-1);   // find random neighbor group
                    neighborGroup = (IslandGroup*)(d->groupData) + neighbors[iNeighbor];  // point to neighbor group
                    fitness2 = (1.f + group->groupfit - neighborGroup->groupfit) * 0.5;
                }
                break;

            case extinctionStrongest:  // extinction probability depends on fitness relative to strongest neighbor
                if (numNeighbors) {
                    // find strongest neighbor
                    float strongest = 0.f;  // fitness*population of strongest neighbor
                    IslandGroup * strongestNeighbor = 0;
                    for (int j = 0; j < numNeighbors; j++) {
                        neighborGroup = (IslandGroup*)(d->groupData) + neighbors[j];  // point to neighbor group
                        if (neighborGroup->groupfit * neighborGroup->nn > strongest) {
                            strongest = neighborGroup->groupfit * neighborGroup->nn;
                            strongestNeighbor = neighborGroup;
                        }
                    }
                    if (strongestNeighbor) {
                        fitness2 = (1.f + group->groupfit - strongestNeighbor->groupfit) * 0.5;
                    }
                }
            }
            
            // extinction
            float relativeSize = float(group->nn) / float(group->nmax); // relative size of group
            if (relativeSize > 1.f) relativeSize = 1.0f;

            float extinctionRate = 
                (1.f-fitness2) * ((1.f-relativeSize) * d->extinctionRate[0] + relativeSize * d->extinctionRate[1])
                +    fitness2  * ((1.f-relativeSize) * d->extinctionRate[2] + relativeSize * d->extinctionRate[3]);

            // kill with probability extinctionRate
            bool kill = d->ran->bernoulli(extinctionRate);

            if (kill) {                // this group goes extinct
                group->age = 0;         // age = 0 means extinct
                d->groupsDied++;
                if (d->surviv > 0.f) { // a fraction survives
                    int32_t numSurvivors = d->ran->binomial(group->nn, d->surviv);
                    group->gAltruism[0] = d->ran->hypergeometric(numSurvivors * 2, group->gAltruism[0], group->nn * 2);
                    group->gAltruism[1] = numSurvivors * 2 - group->gAltruism[0];
                    group->nn = numSurvivors;
                    group->groupfit = 0.f;
#ifdef MORELOCI
                    if (d->locusUsed[locusEndogamy]) {
                        group->gEndogamy[0] = d->ran->hypergeometric(numSurvivors * 2, group->gEndogamy[0], group->nn * 2);
                        group->gEndogamy[1] = numSurvivors * 2 - group->gEndogamy[0];
                    }
                    if (d->locusUsed[locusConformity]) {
                        group->gConformity[0] = d->ran->hypergeometric(numSurvivors * 2, group->gConformity[0], group->nn * 2);
                        group->gConformity[1] = numSurvivors * 2 - group->gConformity[0];
                    }
#endif
                }
                else {   // all die
                    group->gAltruism[0] = group->gAltruism[1] = 0;
#ifdef MORELOCI
                    group->gConformity[0] = group->gConformity[1] = 0;
                    group->gEndogamy[0] = group->gEndogamy[1] = 0;
#endif
                    group->nn = 0;
                    group->groupfit = 0.f;
                }
            }

            // make list of inhabited islands 
            if (group->nn > 0 && inhabitedIslandsNum < d->maxIslands) {
                inhabitedIslands[inhabitedIslandsNum++] = iGroup;
            }
        }

        // group loop 2: migration, recolonization
        for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
            group = (IslandGroup*)(d->groupData) + iGroup;     // point to current group

            errorCheck(d, group);

            // make list of neighbors
            numNeighbors = findNeighbors(d, iGroup, neighbors); // find neighbor groups
            neighborGroup = 0;
            IslandGroup neighborPool;         // all neighbors pooled

            if (numNeighbors) {
                // immigration and colonization use similar patterns
                int migrationOrColonizationPattern;
                if (group->age == 0 || group->nn == 0) {  // group is extinct or empty. recolonize
                    group->age = 0;                      // reset age if group is empty but not extinct
                    migrationOrColonizationPattern = d->colonizationPattern;
                }
                else {   // group is not extinct. calculate immigration
                    migrationOrColonizationPattern = d->immigrationPattern;
                }

                // immigration or colonization
                // find which group to immigrate or colonize from
                switch (migrationOrColonizationPattern) {
                case immigrationCommonPool:
                    // migrants come from common pool
                    neighborGroup = (IslandGroup*)(d->groupData) + d->maxIslands; // point to common pool
                    break;
                case immigrationNeighbor:
                case immigrationProportional:
                case immigrationVacant:
                    // migrants come from random neighbor
                    iNeighbor = d->ran->iRandom(0, numNeighbors - 1);   // find random neighbor group
                    neighborGroup = (IslandGroup*)(d->groupData) + neighbors[iNeighbor];  // point to neighbor group
                    break;
                case immigrationAllNeighbors: {// sum of all neighbors
                    // make weighted sum of all neighbor gene pools
                    memset(&neighborPool, 0, sizeof(neighborPool));  // set to zero
                    double sumGene[6] = {0,0,0,0,0,0};               // sum of neighbor genes * emigrationPotential
                    double sumEmigr = 0;                             // sum of emigrationPotential
                    for (int i = 0; i < numNeighbors; i++) {
                        neighborGroup = (IslandGroup*)(d->groupData) + neighbors[i];  // point to neighbor group
                        sumGene[0] += double(neighborGroup->gAltruism[0]) * neighborGroup->emigrationPotential;
                        sumGene[1] += double(neighborGroup->gAltruism[1]) * neighborGroup->emigrationPotential;
                        sumEmigr += neighborGroup->emigrationPotential;
#ifdef MORELOCI
                        sumGene[2] += neighborGroup->gEndogamy[0] * neighborGroup->emigrationPotential;
                        //sumGene[3] += neighborGroup->gEndogamy[1] * neighborGroup->emigrationPotential;#endif
                        sumGene[4] += neighborGroup->gConformity[0] * neighborGroup->emigrationPotential;
                        //sumGene[5] += neighborGroup->gConformity[1] * neighborGroup->emigrationPotential;
#endif
                    }
                    neighborPool.gAltruism[0] = lround(sumGene[0] / sumEmigr);
                    neighborPool.gAltruism[1] = lround(sumGene[1] / sumEmigr);
                    uint32_t numg = neighborPool.gAltruism[0] + neighborPool.gAltruism[1];
                    if (numg & 1) { // number of genes cannot be odd. make even
                        numg--;
                        neighborPool.gAltruism[0] &= -2;
                        neighborPool.gAltruism[1] &= -2;
                    }
                    neighborPool.nn = numg / 2u;
                    neighborPool.emigrationPotential = sumEmigr / numNeighbors;
#ifdef MORELOCI
                    neighborPool.gEndogamy[0] = lround(sumGene[2] / sumEmigr);
                    if (neighborPool.gEndogamy[0] > numg) neighborPool.gEndogamy[0] = numg;
                    neighborPool.gEndogamy[1] = numg - neighborPool.gEndogamy[0];

                    neighborPool.gConformity[0] = lround(sumGene[4] / sumEmigr);
                    if (neighborPool.gConformity[0] > numg) neighborPool.gConformity[0] = numg;
                    neighborPool.gConformity[1] = numg - neighborPool.gConformity[0];
#endif
                    neighborGroup = &neighborPool;
                    break; }
                case immigrationRandomGroup:
                    if (inhabitedIslandsNum > 0) {  // immigration from random inhabited island
                        int i = d->ran->iRandom(0, inhabitedIslandsNum - 1);         // find random inhabited island
                        iNeighbor = inhabitedIslands[i];
                        neighborGroup = (IslandGroup*)(d->groupData) + iNeighbor;     // point to random group
                    }
                    break;
                case immigrationStrongNeighbor: {
                    // choose among all neighbors with probability proportional to emigrationPotential
                    int32_t sumEmig = 0;  // sum of emigration potential of all neighbors
                    for (int i = 0; i < numNeighbors; i++) {
                        neighborGroup = (IslandGroup*)(d->groupData) + neighbors[i];  // point to neighbor group
                        sumEmig += neighborGroup->emigrationPotential;
                    }
                    if (sumEmig == 0) {
                        neighborGroup = 0; break;
                    }
                    int32_t choose = d->ran->iRandom(0, sumEmig - 1);              // random number for choosing neighbor
                    sumEmig = 0;
                    for (int i = 0; i < numNeighbors; i++) {
                        neighborGroup = (IslandGroup*)(d->groupData) + neighbors[i];  // point to neighbor group
                        sumEmig += neighborGroup->emigrationPotential;
                        if (sumEmig >= choose) {
                            iNeighbor = i;
                            break;
                        }
                    }
                    break;}
                }
                // determine number of immigrants or colonists
                int32_t numMigrants = 0;         // number of immigrants or colonists

                if (group->age == 0) {  // group is extinct. recolonize
                    numMigrants = d->ran->poisson(d->colonySize);
                }
                else {  // group is not extinct. calculate number of immigrants
                    float migrationRate = d->migrationRate[0];
#ifdef MORELOCI     
                    if (d->locusUsed[locusEndogamy]) {
                        // adjust migration rate to fraction of endogamists.
                        // find number of phenotypic endogamists
                        int32_t endogamyGenotypes[3];
                        // split gene pool into genotypes
                        combineGenes(group->gEndogamy, endogamyGenotypes, d->ran);
                        // find number of phenotypic endogamists
                        int32_t endogamists = endogamyGenotypes[2];
                        switch (d->dominance[locusEndogamy]) {
                        case recessive:              // endogamy recessive
                            break;
                        case dominant:               // endogamy dominant
                            endogamists += endogamyGenotypes[1];
                            break;
                        case incompleteDominant:      // half dominant
                            endogamists += endogamyGenotypes[1] / 2;
                        }
                        float endogamy = 0.f;        // fraction of endogamists
                        if (group->nn > 0) {
                            endogamy = float(endogamists) / float(group->nn);
                        }
                        // migrationRate as a function of endogamy
                        migrationRate = endogamy * d->migrationRate[1] + (1.f - endogamy) * d->migrationRate[0];
                    }
#endif
                    // find number of immigrants according to immigrationPattern
                    switch (d->immigrationPattern) {
                    case immigrationCommonPool:      // from common pool
                        numMigrants = d->ran->poisson(migrationRate * group->nmax);
                        break;
                    case immigrationNeighbor:        // from random neighbor group
                    case immigrationRandomGroup:     // from random group
                    case immigrationAllNeighbors:    // from all neighbors
                        if (neighborGroup) {
                            //numMigrants = d->ran->poisson(migrationRate * neighborGroup->emigrationPotential);
                            numMigrants = d->ran->binomial(neighborGroup->emigrationPotential, migrationRate);
                        }
                        break;
                    case immigrationStrongNeighbor:  // from neighbor selected according to emigrationPotential
                        if (neighborGroup) {
                            //numMigrants = d->ran->poisson(migrationRate * neighborGroup->nn);
                            numMigrants = d->ran->binomial(neighborGroup->nn, migrationRate);
                        }
                        break;
                    case immigrationProportional:    // prop. w. population, from neighbor
                        //numMigrants = d->ran->poisson(migrationRate * group->nn);
                        numMigrants = d->ran->binomial(neighborGroup->nn, migrationRate);
                        break;
                    case immigrationVacant:          // prop. w. vacant capacity, from neighbor
                        numMigrants = group->nmax - group->nn + 1;
                        if (numMigrants < 0) numMigrants = 0;
                        else {
                            numMigrants = d->ran->binomial(numMigrants, migrationRate);
                            //numMigrants = d->ran->poisson(migrationRate * numMigrants);
                        }
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
#ifdef MORELOCI
                    if (d->locusUsed[locusEndogamy]) {
                        gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborGroup->gEndogamy[0], neighborGroup->nn * 2);
                        gmig1 = numMigrants * 2 - gmig0;
                        group->gEndogamy[0] += gmig0;       // add colonist genes to current group
                        group->gEndogamy[1] += gmig1;
                        neighborGroup->gEndogamy[0] -= gmig0; // subtract colonist genes from neighbor group
                        neighborGroup->gEndogamy[1] -= gmig1;
                    }
                    if (d->locusUsed[locusConformity]) {
                        gmig0 = d->ran->hypergeometric(numMigrants * 2, neighborGroup->gConformity[0], neighborGroup->nn * 2);
                        gmig1 = numMigrants * 2 - gmig0;
                        group->gConformity[0] += gmig0;     // add colonist genes to current group
                        group->gConformity[1] += gmig1;
                        neighborGroup->gConformity[0] -= gmig0; // subtract colonist genes from neighbor group
                        neighborGroup->gConformity[1] -= gmig1;
                    }
#endif
                    neighborGroup->nn -= numMigrants;
                    d->migrantsTot += numMigrants;         // count total number of migrants
                }
            }

            // limit population
            if (group->nn > group->nmax) {
                int32_t nn2 = group->nmax * 2; // number of genes after population limitation
                int32_t g0 = d->ran->hypergeometric(nn2, group->gAltruism[0], group->nn * 2);
                int32_t g1 = nn2 - g0;
                group->gAltruism[0] = g0;
                group->gAltruism[1] = g1;
#ifdef MORELOCI
                if (d->locusUsed[locusEndogamy]) {
                    g0 = d->ran->hypergeometric(nn2, group->gEndogamy[0], group->nn * 2);
                    g1 = nn2 - g0;
                    group->gEndogamy[0] = g0;  // add migrant genes to current group
                    group->gEndogamy[1] = g1;
                }
                if (d->locusUsed[locusConformity]) {
                    g0 = d->ran->hypergeometric(nn2, group->gConformity[0], group->nn * 2);
                    g1 = nn2 - g0;
                    group->gConformity[0] = g0; // add migrant genes to current group
                    group->gConformity[1] = g1;
                }
#endif
                group->nn = group->nmax;
            }

            errorCheck(d, group);    // may be removed
        }

        // group loop 3: update global gene pool and statistics
        for (iGroup = 0; iGroup < d->maxIslands; iGroup++) {
            group = (IslandGroup*)(d->groupData) + iGroup;     // point to current group
            d->genePool[locusAltruism][0] += group->gAltruism[0];
            d->genePool[locusAltruism][1] += group->gAltruism[1];
#ifdef MORELOCI
            d->genePool[locusEndogamy][0] += group->gEndogamy[0];
            d->genePool[locusEndogamy][1] += group->gEndogamy[1];
            d->genePool[locusConformity][0] += group->gConformity[0];
            d->genePool[locusConformity][1] += group->gConformity[1];
#endif
            d->totalPopulation += group->nn;
            d->totalPhenotypes[locusAltruism] += group->altruists;
            if (group->groupfit > 0.5f) d->altruistGroups++;
            if (group->nn >= d->minGroupSize) d->inhabitedGroups++;
            group->age++;
        }

        // make migrant pool
        if (d->migrationTopology == topologyCommonPool || d->colonizationPattern == colonizationCommonPool) {  // migrant pool used
            // migrant pool placed last in memory block
            IslandGroup * migrantPool = (IslandGroup*)(d->groupData) + d->maxIslands;
            double reduction = 1.;
            if (d->totalPopulation > 100000000) {  // avoid overflow
                reduction = 100000000. / d->totalPopulation;
            }
            migrantPool->gAltruism[0] = int32_t(d->genePool[locusAltruism][0] * reduction);
            migrantPool->gAltruism[1] = int32_t(d->genePool[locusAltruism][1] * reduction);
#ifdef MORELOCI
            migrantPool->gEndogamy[0] = int32_t(d->genePool[locusEndogamy][0] * reduction);
            migrantPool->gEndogamy[1] = int32_t(d->genePool[locusEndogamy][1] * reduction);
            migrantPool->gConformity[0] = int32_t(d->genePool[locusConformity][0] * reduction);
            migrantPool->gConformity[1] = int32_t(d->genePool[locusConformity][1] * reduction);
#endif
            migrantPool->nn = (migrantPool->gAltruism[0] + migrantPool->gAltruism[1]) / 2;
            migrantPool->nmax = migrantPool->nn;
            migrantPool->groupfit = 0.5f;
        }
        // update statistics counters
        statisticsUpdate(d);
        d->generations++;  // count generations

        // check if finished
        checkStopCriterion(d);
    }
}
