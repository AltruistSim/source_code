/***********************  parameter_file_io.cpp   *****************************
* Author:        Agner Fog
* Date created:  2023-04-15
* Last modified: 2023-12-31
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the procedures and lists for reading and writing 
* parameter files.
*
* (c) Copyright 20232024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

/*
struct ParameterDef {
    uint32_t type;           // 0: header, 1: uint32_t, 2: float, 4: enumerated names
    uint32_t num;            // array size
    uint32_t offset;         // offset into AltruData (this works as a member pointer)
    const char * name;       // name of parameter
};
*/

// list of parameter definitions
ParameterDef const parameterDefinitions[] {
    {0, 1, 0, "Model"},
    {4, 1, altruDataOffset(iModel), "model"},

    {0, 1, 0, "Geography"},
    {1, 1, altruDataOffset(totArea), "totalArea"},
    {1, 1, altruDataOffset(maxIslands), "maxNumGroups"},
    {1, 1, altruDataOffset(nMaxPerDeme), "maxPerGroup"},
    {2, 2, altruDataOffset(carryingCapacity), "carryingCapacity"},
    {2, 1, altruDataOffset(carryingCapacityStandardDeviation), "carryingCapacityStandardDeviation"},
    {1, 1, altruDataOffset(territorySizeMax), "maxTerritoryArea"},
    {1, 1, altruDataOffset(territorySizeMin), "minTerritoryArea"},
    {1, 1, altruDataOffset(minGroupSize), "minGroupSize"},
    {1, 1, altruDataOffset(colonySize), "recolonizationGroupSize"},

    {0, 1, 0, "Migration"},
    {2, 2, altruDataOffset(migrationRate), "migrationRate"},
    {1, 1, altruDataOffset(emigrationPattern), "emigrationPattern"},
    {1, 1, altruDataOffset(immigrationPattern), "immigrationPattern"},
    {1, 1, altruDataOffset(colonizationPattern), "recolonizationPattern"},
    {1, 1, altruDataOffset(migrationTopology), "migrationTopology"},

    {0, 1, 0, "Loci"},
    {1, maxLoci,   altruDataOffset(locusUsed), "locusUsed"},
    {1, maxLoci,   altruDataOffset(dominance), "dominance"},
    {2, maxLoci,   altruDataOffset(fg0), "initialFractionOfMutant"},
    {2, maxLoci*2, altruDataOffset(murate), "mutationRates"},

    {0, 1, 0, "Individual fitness"},
    {2, 8, altruDataOffset(fit), "individualFitness"},
    {2, 4, altruDataOffset(fit2), "extraLociFitness"},
    {2, 1, altruDataOffset(growthRate), "growthRate"},
    {1, 1, altruDataOffset(selectionModel), "selectionModel"},

    {0, 1, 0, "Group properties"},
    {2, 1, altruDataOffset(surviv), "survivalRate"},
    {2, 1, altruDataOffset(warIntensity), "warIntensity"},
    {2, 4, altruDataOffset(extinctionRate), "demeExtinctionRates"},
    {1, 1, altruDataOffset(extinctionPattern), "extinctionPattern"},
    {1, 1, altruDataOffset(warPattern), "warPattern"},
    {1, 1, altruDataOffset(fitfunc), "groupFitnessFunction"},
    {2, 1, altruDataOffset(fitExpo), "groupFitnessExponent"},
    {2, 1, altruDataOffset(leaderAdvantage), "leaderAdvantage"},
    {2, 1, altruDataOffset(leaderSelection), "leaderSelection"},
    {1, 1, altruDataOffset(haystackPeriod), "haystackPeriod"},
    {1, 1, altruDataOffset(mixingPeriod), "mixingPeriod"},

    {0, 1, 0, "Run control"},
    {1, 1, altruDataOffset(seed), "randomSeed"},
    {1, 1, altruDataOffset(minimumGenerations), "minimumGenerations"},
    {1, 1, altruDataOffset(maximumGenerations), "maximumGenerations"},
    {1, 1, altruDataOffset(stopCriterion), "stopCriterion"},
    {2, 1, altruDataOffset(stopCriterionDegree), "stopCriterionDegree"},
    {1, 1, altruDataOffset(delayms), "delayScreen"},
    {1, 1, altruDataOffset(sweepsUsed), "sweepsUsed"},
    {1, maxSweeps, altruDataOffset(sweepType), "sweepType"},
    {1, maxSweeps, altruDataOffset(sweepParameter), "sweepParameter"},
    {2, maxSweeps, altruDataOffset(sweepStartValue), "sweepStartValue"},
    {2, maxSweeps, altruDataOffset(sweepEndValue), "sweepEndValue"},
    {2, maxSweeps, altruDataOffset(sweepStep), "sweepStep"},

    {0, 1, 0, "Output control"},
    {1, 1, altruDataOffset(bOutOptions), "fileOutputOptions"},
    {0,0,0,0}
};

// size of parameterDefinitions
uint32_t const parameterDefinitionsSize = ARR_LEN(parameterDefinitions);


// functions for reading and writing files
void Altruist::readParameterFile(QString filename) {
    QFile file(filename);
    bool success = file.open(QIODeviceBase::ReadOnly);
    if (!success) {
        errorMessage((QString("cannot read file ") + filename));
        return;
    }
    const int textlen = 256;
    char text[textlen];                // temporary text buffer
    char * sep;                        // point to token separator
    char * a, *b;                      // point to begin and end of token
    int par, i, j;                     // loop counters
    int32_t ival;                      // integer value
    float fval;                        // float value
    ParameterDef const * pardef;       // parameter definition record
    int8_t * field;                    // point to field in d
    bool modelKnown = false;           // model name has been read
    bool found = false;                // parameter name found in list

    // loop through lines in file
    while (!file.atEnd()) {
        file.readLine(text, textlen);
        if (text[0] == '[') continue;  // skip header
        sep = strchr(text, '=');
        if (sep == nullptr) continue;
        *sep = 0;                      // end first token
        // trim name. a = begin of name, b = end of name
        a = text; b = sep;
        while (a < b && *a <= ' ') a++;
        while (b > a && b[-1] <= ' ') b--;
        *b = 0;
        found = false;

        // loop through two parameter lists, the global and the model-specific list, to search for name a
        for (int parlist = 0; parlist < 2; parlist++) {
            if (parlist == 0) {
                pardef = parameterDefinitions;
            }
            else {
                if (modelKnown && d.iModel < models.getNum()) {
                    ModelDescriptor * m = models.getModel(d.iModel);
                    pardef = m->extraParameters;
                }
                else pardef--;
            }
            // loop through parameter list
            for (; pardef->num != 0 || pardef->name != nullptr; pardef++) {
                if (pardef->offset >= sizeof(AltruData)) continue;   // check if out of range
                if (pardef->name == nullptr) continue;
                if (strcasecmp(a, pardef->name) == 0 && pardef->type != 0) { // matching name found
                    found = true;
                    field = (int8_t*)&d + pardef->offset;
                    for (i = 0; i < pardef->num; i++) {    // loop through list of values
                        // trim parameter value
                        a = sep + 1;
                        while (a < b && *a <= ' ') a++;
                        switch (pardef->type) {
                        case 1:  // int32
                            ival = atoi(a);                // read integer value
                            *(int32_t*)field = ival;       // store value in d field
                            field += sizeof(int32_t);
                            break;
                        case 2:  // float
                            fval = atof(a);                // read float value
                            *(float*)field = fval;         // store value in d field
                            field += sizeof(float);
                            break;
                        case 3:  // int32 as bool
                            ival = atoi(a);  // read integer value
                            *(int32_t*)field = ival != 0;  // store value in d field
                            field += sizeof(int32_t);
                            break;
                        case 4:  // enumerated names
                            if (!sep) break;
                            b = a + strlen(a);
                            while (b > a && (b[-1] <= ' ' || b[-1] == ',')) b--;
                            *b = 0;
                            if (strcmp(pardef->name, "model") == 0) {
                                for (j = 0; j < models.getNum(); j++) {
                                    ModelDescriptor * m = models.getModel(j);
                                    if (strcasecmp(a, m->name) == 0) {
                                        d.iModel = j;      // matching model name found
                                        modelKnown = true;
                                        break;
                                    }
                                }
                                if (j == models.getNum()){ // model name not found
                                    errorMessage(QString("Unknown model ") + a);
                                }
                            }
                            else {
                                errorMessage("Unknown enumeration in parameterDefinitions");
                            }
                            break;
                        }
                        sep = strchr(a, ',');              // next value on line
                        if (sep == nullptr) break;
                    }
                    break;                                 // stop searching                
                }
            }
            if (found) break;
        }
        if (!found) {
            errorMessage(QString("Unknown parameter ") + a);
        }
    }
    file.close();
    (*(models.getModel(d.iModel)->initFunction))(&d, model_initialize); // initialize new model
}


void Altruist::writeParameterFile(QString filename) {
    QFile file(filename);
    bool success = file.open(QIODeviceBase::WriteOnly, (QFileDevice::Permission)0x6666);
    if (!success) {
        errorMessage((QString("cannot write file ") + filename));
        return;
    }
    const int textlen = 256;
    char text[textlen];                // temporary text buffer
    ParameterDef const * pardef;
    int i;                             // loop counter
    int n;                             // character counter
    int8_t * field;                    // point to field in d
    // loop through two parameter lists, the global and the model-specific list
    for (int parlist = 0; parlist < 2; parlist++) {
        if (parlist == 0) {
            pardef = parameterDefinitions;
        }
        else {
            if (d.iModel < models.getNum()) {
                ModelDescriptor * m = models.getModel(d.iModel);
                pardef = m->extraParameters;
            }
        }
        // loop through parameter list
        for (; pardef->num != 0 || pardef->name != nullptr; pardef++) {
            if (pardef->offset >= sizeof(AltruData)) continue;  // check if out of range
            field = (int8_t*)&d + pardef->offset; // pointer to field in d
            // write line according to field type
            switch (pardef->type) {
            case 0:  // header
                snprintf(text, textlen, "\n[%s]\n", pardef->name);
                break;
            case 1:  // int
                n = snprintf(text, textlen, "%s=%i", pardef->name, *(int32_t*)field);
                for (i = 1; i < pardef->num; i++) {   // if more than one value
                    field += sizeof(int32_t);
                    n += snprintf(text + n, textlen - n, ", %i", *(int32_t*)field);
                }
                snprintf(text + n, textlen - n, "\n");
                break;
            case 2:  // float
                n = snprintf(text, textlen, "%s=%G", pardef->name, *(float*)field);
                for (i = 1; i < pardef->num; i++) {   // if more than one value
                    field += sizeof(float);
                    n += snprintf(text + n, textlen - n, ", %G", *(float*)field);
                }
                snprintf(text + n, textlen - n, "\n");
                break;
            case 3:  // int32 as bool
                n = snprintf(text, textlen, "%s=%i", pardef->name, *(int32_t*)field != 0);
                for (i = 1; i < pardef->num; i++) {   // if more than one value
                    field += sizeof(int32_t);
                    n += snprintf(text + n, textlen - n, ", %i", *(int32_t*)field != 0);
                }
                snprintf(text + n, textlen - n, "\n");
                break;
            case 4:  // enumerated names
                if (strcmp(pardef->name, "model") == 0) {
                    ModelDescriptor * m = models.getModel(d.iModel);
                    snprintf(text, textlen, "%s=%s\n", pardef->name, m->name);
                }
                else {
                    errorMessage("Unknown enumeration in parameterDefinitions");
                }
                break;
            default:
                continue;
            }
            file.write(text);
        }
    }
    file.close();
}
