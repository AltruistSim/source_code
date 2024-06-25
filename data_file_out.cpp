/**************************  data_file_out.cpp   ******************************
* Author:        Agner Fog
* Date created:  2023-08-31
* Last modified: 2023-12-31
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines functions for writing result data files
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

void Worker::fileOutStart() {
    // start writing data output file
    char text[1024];                                       // text buffer
    if (!d->makeOutputFile) return;
    if (outputFile.isOpen()) outputFile.close();           // close any previously opened file
    outputFile.setFileName(d->outputFilePath + "/" + d->outputFileName);
    // delete previous file with same name only if this is the first result set since file name change:
    bool ok = outputFile.open(resultSets ? QIODeviceBase::Append : QIODeviceBase::WriteOnly);
    if (!ok) {
        errors.reportError((QString("Failed to open data output file ") + outputFile.fileName()).toUtf8());
        return;    
    }
    const char stars[] = "\n#*******************************************************************************\n";
    if (outputFile.size() > 1) {
        // appending to an existing file
        outputFile.write(stars);    
    }
    // write title, result set number, model name, time
    QDateTime date = QDateTime::currentDateTime();
    QString datestring = date.toString("yyyy-MM-dd hh:mm:ss");
    snprintf(text, sizeof(text), "\n# %s number %i\n# Model: %s.  %s\n", d->outputTitle.toUtf8().data(), 
        ++resultSets, models.getModel(d->iModel)->name, datestring.toUtf8().data());
    outputFile.write(text);

    if (d->bOutOptions & 1) {   // write initial parameters
        // geography
        snprintf(text, sizeof(text), "\n# Geography: maxNumGroups = %i, maxPerGroup = %i"
            "\n# carryingCapacity = %.1f, %.1f, carryingCapacityStandardDeviation = %.2f"
            "\n# minGroupSize = %i, recolonizationGroupSize = %i",
            d->maxIslands, d->nMaxPerDeme, d->carryingCapacity[0], d->carryingCapacity[1],
            d->carryingCapacityStandardDeviation, d->minGroupSize, d->colonySize);
        outputFile.write(text);
        // migration
        snprintf(text, sizeof(text), "\n# Migration: MigrationRate = %.3G, %.3G\n# immigrationPattern = %i"
            "recolonizationPattern = %i, migrationTopology = %i",
            d->migrationRate[0], d->migrationRate[1], 
            d->immigrationPattern, d->colonizationPattern, d->migrationTopology);
        outputFile.write(text);
        // loci
        for (int i = 0; i < maxLoci; i++) {
            if (d->locusUsed[i]) {
                snprintf(text, sizeof(text), "\n# Locus %s: dominance = %i"
                    "\n# initial fraction %.3f, mutation rates %.2G %.2G",
                    d->sLocusName[i], d->dominance[i], d->fg0[i],
                    d->murate[i][0], d->murate[i][1]);
                outputFile.write(text);
            }
        }
        // individual fitness
        outputFile.write("\n# Fitness: ");
        int numFit = d->locusUsed[locusConformity] ? 8 : 4;// number of fitness measures to print
        text[0] = 0;
        for (int i = 0; i < numFit; i++) {
            snprintf(text + strlen(text), sizeof(text), "%.3f%c ", d->fit[i], (i < numFit - 1 ? ',' : ' '));        
        }
        outputFile.write(text);
        snprintf(text, sizeof(text), "\n# growthRate = %.3f, selectionModel = %i",
            d->growthRate, d->selectionModel);
        outputFile.write(text);
        // group fitness
        snprintf(text, sizeof(text), "\n# Group properties: groupFitnessFunction = %i, groupFitnessExponent = %.2G"
            "\n# extinctionPattern = %i, warPattern = %i, survivalDegree = %.3G, warIntensity = %.3G"
            "\n# demeExtinctionRates = %.3f, %.3f, %.3f, %.3f",
            d->fitfunc, d->fitExpo, d->extinctionPattern, d->warPattern, d->surviv, d->warIntensity,
            d->extinctionRate[0], d->extinctionRate[1], d->extinctionRate[2], d->extinctionRate[3]);
        outputFile.write(text);
        snprintf(text, sizeof(text), "\n# Run control: randomSeed = %i"
            ", minimumGenerations = %i, maximumGenerations = %i"
            "\n# stopCriterion = %i, stopCriterionDegree = %.3f",
            d->seed, d->minimumGenerations, d->maximumGenerations,
            d->stopCriterion, d->stopCriterionDegree);
        outputFile.write(text);
    }
    // write any parameter sweeps
    const char * sweepTypeNames[5] = {
        "none", "linear sweep", "logarithmic sweep", "linear search", "logarithmic search"};
    for (int iSweep = maxSweeps-1; iSweep >= 0; iSweep--) {
        if (d->sweepsUsed & (1 << iSweep)) {
            snprintf(text, sizeof(text), "\n# %s of %s:"
                "\n# sweepStartValue = %.3G, sweepEndValue = %.3G",
                sweepTypeNames[d->sweepType[iSweep]],
                sweepParameterList[d->sweepParameter[iSweep]].name,
                d->sweepStartValue[iSweep], d->sweepEndValue[iSweep]);
            if (d->sweepType[iSweep] == 1) {
                snprintf(text + strlen(text), sizeof(text), ", step %.3G", d->sweepStep[iSweep]);
            }
            if (d->sweepType[iSweep] == 2) {
                snprintf(text + strlen(text), sizeof(text), ", %.1f steps per decade", d->sweepStep[iSweep]);
            }
            outputFile.write(text);
        }
    }
    // make column headings
    outputFile.write("\n\n");
    // names of sweep parameters
    for (int iSweep = maxSweeps-1; iSweep >= 0; iSweep--) {
        if (d->sweepsUsed & (1 << iSweep)) {
            snprintf(text, sizeof(text), "%s", sweepParameterList[d->sweepParameter[iSweep]].name);
            // replace spaces
            char * s = 0;
            do {
                s = strchr(text, ' ');
                if (s) *s = '_';            
            }
            while (s);
            outputFile.write(text);
            outputFile.write(", ");
        }
    }
    if (d->bOutOptions & 2) outputFile.write("generations, ");
    if (d->bOutOptions & 4) outputFile.write("inhabited_islands, ");
    if (d->bOutOptions & 8) { // gene fractions for all loci
        text[0] = 0;
        for (int i = 0; i < maxLoci; i++) {
            if (d->locusUsed[i]) {
                snprintf(text + strlen(text), sizeof(text), "%s_gene_fraction, ", d->sLocusName[i]);
            }
        }
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x10) outputFile.write("altrusts, ");
    if (d->bOutOptions & 0x20) outputFile.write("mutations, ");
    if (d->bOutOptions & 0x40) outputFile.write("migrants, ");
    if (d->bOutOptions & 0x100) outputFile.write("altruism_demes, ");
    if (d->bOutOptions & 0x200) outputFile.write("extinctions, ");
    if (d->bOutOptions & 0x400) outputFile.write("result, ");
    if (d->bOutOptions & 0x800) {  // model specific variable names
        for (int i = 0; i < d->nUser; i++) {
            const char * uname = "?";
            if (d->userDataNames) uname = d->userDataNames[i];
            snprintf(text, sizeof(text), "%s, ", uname);
            outputFile.write(text);
        }
    }
    if (d->bOutOptions & 0x10000) outputFile.write("limit1, limit2, ");
    outputFile.flush();
}

void Worker::fileOutGeneration() {
    // write single generation to data output file
    if (!d->makeOutputFile) return;
    // interval for output
    if (d->fileOutInterval <= 0) return;
    if (d->generations % d->fileOutInterval != 0) return;
    // make sure file is open
    if (!outputFile.isOpen()) fileOutStart();
    if (!outputFile.isOpen()) return;
    char text[256];
    outputFile.write("\n");

    // write any sweep parameters
    for (int iSweep = maxSweeps-1; iSweep >= 0; iSweep--) {
        if (d->sweepsUsed & (1 << iSweep)) {
            snprintf(text, sizeof(text), "%.5G, ", lastParameter[iSweep]);
            outputFile.write(text);
        }
    }
    if (d->bOutOptions & 2) {          // generations
        snprintf(text, sizeof(text), "%u, ", (uint32_t)d->generations);
        outputFile.write(text);
    }
    if (d->bOutOptions & 4) {          // inhabited islands
        snprintf(text, sizeof(text), "%i, ", d->inhabitedDemes);
        outputFile.write(text);
    }
    if (d->bOutOptions & 8) {          // gene fractions for all loci
        text[0] = 0;
        for (int i = 0; i < maxLoci; i++) {
            if (d->locusUsed[i]) {
                snprintf(text + strlen(text), sizeof(text), "%.3f, ", d->geneFraction[i]);
            }
        }
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x10) {       // fraction of phenotypic altruists
        snprintf(text, sizeof(text), "%.3f, ", d->fAltru);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x20) {       // mutations
        snprintf(text, sizeof(text), "%i, ", d->sumMutations);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x40) {       // migrants
        snprintf(text, sizeof(text), "%u, ", d->migrantsTot);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x100) {      // altruism demes
        snprintf(text, sizeof(text), "%u, ", d->altruistDemes);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x200) {      // extinctions
        snprintf(text, sizeof(text), "%u, ", d->demesDied);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x400) {      // simulation result not available yet
        outputFile.write("0, ");
    }
    if (d->bOutOptions & 0x800) {      // model specific parameter
        fileOutModelSpecificResults();
    }
    if (d->bOutOptions & 0x10000) {    // search results (limit1, limit2)
        outputFile.write(", , ");
    }
    outputFile.flush();
}

void Worker::fileOutModelSpecificResults() {
    // write all model specific result variables
    char text[32];
    for (int i = 0; i < d->nUser; i++) {
        if (d->bUserDataType & 1 << i) {
            snprintf(text, sizeof(text), "%.5G, ", d->userData[i].f);
        }
        else {
            snprintf(text, sizeof(text), "%u, ", d->userData[i].i);
        }
        outputFile.write(text);
    }
}

void Worker::fileOutSimulationFinished() {
    // write result of a single simulation run (generation loop)
    if (!d->makeOutputFile) return;
    // make sure file is open
    if (!outputFile.isOpen()) fileOutStart();
    if (!outputFile.isOpen()) return;
    char text[256];
    outputFile.write("\n");

    // write any sweep parameters
    for (int iSweep = maxSweeps-1; iSweep >= 0; iSweep--) {
        if (d->sweepsUsed & (1 << iSweep)) {
            snprintf(text, sizeof(text), "%.5G, ", lastParameter[iSweep]);
            outputFile.write(text);
        }
    }
    if (d->bOutOptions & 2) {          // generations
        snprintf(text, sizeof(text), "%u, ", (uint32_t)d->generations);
        outputFile.write(text);
    }
    if (d->bOutOptions & 4) {          // inhabited islands
        snprintf(text, sizeof(text), "%i, ", d->inhabitedDemes);
        outputFile.write(text);
    }
    if (d->bOutOptions & 8) {          // gene fractions for all loci
        text[0] = 0;
        for (int i = 0; i < maxLoci; i++) {
            if (d->locusUsed[i]) {
                snprintf(text + strlen(text), sizeof(text), "%.3f, ", d->geneFraction[i]);
            }
        }
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x10) {       // fraction of phenotypic altruists
        snprintf(text, sizeof(text), "%.3f, ", d->fAltru);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x20) {       // mutations
        snprintf(text, sizeof(text), "%i, ", d->sumMutations);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x40) {       // migrants
        snprintf(text, sizeof(text), "%G, ", double(d->sumMigrants));
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x100) {      // altruism demes
        snprintf(text, sizeof(text), "%u, ", d->altruistDemes);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x200) {      // extinctions
        snprintf(text, sizeof(text), "%u, ", d->sumExtinctions);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x400) {      // simulation result not available yet
        snprintf(text, sizeof(text), "%u, ", d->simulationResult);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x800) {      // model specific statistics
        fileOutModelSpecificResults();
    }
    if (d->bOutOptions & 0x10000) {    // search results (limit1, limit2)
        outputFile.write(", , ");
    }
    // write time used
    if (d->bOutOptions & 0x2000) {    
        snprintf(text, sizeof(text), "\n# Time used: %.3G seconds\n", double(d->timeUsed) * 0.001);
        outputFile.write(text);
    }
    if (d->sweepsUsed) outputFile.flush();
    else outputFile.close();
}


void Worker::fileOutSweepNext() {
    // parameter sweep or search
    if (!d->makeOutputFile) return;
    // make sure file is open
    if (!outputFile.isOpen()) fileOutStart();
    if (!outputFile.isOpen()) return;
    char text[256];                    // text buffer
    outputFile.write("\n");
    // write sweep parameters
    for (int iSweep = maxSweeps-1; iSweep >= 0; iSweep--) {
        if (d->sweepsUsed & (1 << iSweep)) {
            snprintf(text, sizeof(text), "%.5G, ", lastParameter[iSweep]);
            outputFile.write(text);
        }
    }
    if (d->bOutOptions & 2) {          // generations
        snprintf(text, sizeof(text), "%u, ", (uint32_t)d->generations);
        outputFile.write(text);
    }
    if (d->bOutOptions & 4) {          // inhabited islands
        snprintf(text, sizeof(text), "%i, ", d->inhabitedDemes);
        outputFile.write(text);
    }
    if (d->bOutOptions & 8) {          // gene fractions for all loci
        text[0] = 0;
        for (int i = 0; i < maxLoci; i++) {
            if (d->locusUsed[i]) {
                snprintf(text + strlen(text), sizeof(text), "%.3f, ", d->geneFraction[i]);
            }
        }
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x10) {       // fraction of phenotypic altruists
        snprintf(text, sizeof(text), "%.3f, ", d->fAltru);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x20) {       // mutations
        snprintf(text, sizeof(text), "%i, ", d->sumMutations);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x40) {       // migrants
        snprintf(text, sizeof(text), "%G, ", double(d->sumMigrants));
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x100) {      // altruism demes
        snprintf(text, sizeof(text), "%u, ", d->altruistDemes);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x200) {      // extinctions
        snprintf(text, sizeof(text), "%u, ", d->demesDied);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x400) {      // simulation result
        snprintf(text, sizeof(text), "%u, ", lastResult);
        outputFile.write(text);
    }
    if (d->bOutOptions & 0x800) {      // model specific statistics
        fileOutModelSpecificResults();
    }
    if (d->bOutOptions & 0x10000) {    // search results (limit1, limit2)
        if (loops[0].type == 2 && loopsFinished && limitListNum > 0) {
            snprintf(text, sizeof(text), "%.4G, %.4G, ", limitList[limitListNum-1].limit1, limitList[limitListNum-1].limit2);
            outputFile.write(text);
        }
        else {                         // no search results       
            outputFile.write(", , ");
        }
    }
    outputFile.flush();
}

void Worker::fileOutSweepFinished() {
    // parameter sweep finished. write statistics
    fileOutSweepNext();
    // write time used
    char text[256];                    // text buffer
    if (d->bOutOptions & 0x2000) {     // write time used
        if (loopsFinished == (1 << maxSweeps) - 1) { // all loops finished
            snprintf(text, sizeof(text), "\n# Time used total: %.3fs", double(d->timeUsedLoops) * 0.001);
            outputFile.write(text);
            QDateTime date = QDateTime::currentDateTime();
            QString datestring = date.toString("yyyy-MM-dd hh:mm:ss");
            snprintf(text, sizeof(text), "\n# Time finished: %s\n",
                datestring.toUtf8().data());
            outputFile.write(text);
        }
    }
    if (loopsFinished == (1 << maxSweeps) - 1) outputFile.close(); // all loops finished
    else outputFile.flush();
}
