/****************************  altruist.cpp   *********************************
* Author:        Agner Fog
* Date created:  2023-03-30
* Last modified: 2023-12-31
* Version:       3.00.00
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the general user interface with menus and dialogs.
* Each simulation model is defined in a separate C++ file.
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

// main application entry
int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    Altruist win;
    win.show();
    return app.exec();
}

// global list of implemented models
ModelDescriptorList models;

// global error reporter
ErrorReporter errors;

// To trace an error, set a breakpoint here and trace the call stack:
void ErrorReporter::reportError(const char * text) {       // report an error
    if (num++ == 0) strncpy(errorText, text, 1024);        // preserve only first error text
}

void ErrorReporter::clear() {                              // clear error message
    num = 0;
    errorText[0] = '?'; errorText[1] = 0;
}

// constructor for main Altruist class
Altruist::Altruist(QWidget *parent) : 
    QMainWindow(parent),
        menuFile("&File"),
    menuParameters("&Parameters"),
    menuRun("&Run")
{
    ui.setupUi(this);                                      // set up user interface
    memset(&d, 0, sizeof(d));                              // clear all data
    initializeData();                                      // initialize data
    worker = new Worker(&d);                               // create worker object
    worker->moveToThread(&workerThread);                   // put worker in its own thread
    connect(&workerThread, &QThread::finished, worker, &QObject::deleteLater);
    connect(this, &Altruist::doWorkSignal, worker, &Worker::doWorkSlot); // signal from main to worker
    connect(worker, &Worker::resultReadySignal, this, &Altruist::resultReadySlot); // signal from worker to main
    workerThread.start();                                  // start message loop in worker thread

    readParameterFile(lastParameterFile);                  // read last parameters from last.altru
    setupMenus();                                          // make menus
    graphicsView = new AltruistView(this);

    // call update in message loop
    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &Altruist::update);
    timer->start();
}

Altruist::~Altruist() {
    writeParameterFile(lastParameterFile);                 // write parameters to file last.altru
    d.runState = state_idle;  d.sweepState = state_idle;
    workerThread.quit();                                   // stop worker thread
    thread()->msleep(1);
    if (workerThread.isRunning()) {                        // make sure it stops even if it is too busy for the message loop
        thread()->msleep(100);
        if (workerThread.isRunning()) workerThread.terminate();
        if (workerThread.isRunning()) workerThread.wait();
    }
    // free allocated memory
    if (d.demeData) delete[] d.demeData;
    for (int i = 0; i < numExtraBuffers; i++) {
        if (d.extraBuffer[i]) delete[] d.extraBuffer[i];
    }

    //delete worker; // deleted automatically by parent workerThread
    // delete graphicsView;
}

void Altruist::initializeData() {
    models.sort();                     // sort installed models by name

    d.bGeographyParametersUsed = -1;
    d.totArea = 1000;
    d.bImmigrationPattern = 0b111111;
    d.emiPatName[0] = "constant";
    d.emiPatName[1] = "proportional w. population";
    d.emiPatName[2] = "prop. w. excess population";
    d.emiPatName[3] = "prop. w. group fitness * population";
    d.immiPatName[0] = "from common gene pool";
    d.immiPatName[1] = "from random neighbor";
    d.immiPatName[2] = "prop. w. population, from neighbor";
    d.immiPatName[3] = "prop. w. vacant capacity, from neighbor";
    d.immiPatName[4] = "from all neighbor groups";
    d.immiPatName[5] = "from random inhabited island";
    d.immiPatName[6] = "from strong neighbor group";
    d.bColonPat = 0b1110011;
    d.colonPatName[0] = "from a common gene pool";        
    d.colonPatName[1] = "from random neighbor group";
    d.colonPatName[2] = "?";
    d.colonPatName[3] = "?";
    d.colonPatName[4] = "from all neighbor groups";
    d.colonPatName[5] = "from a random inhabited island";
    d.colonPatName[6] = "from strong neighbor group";
    d.bTopology = 0b01111111;
    d.migTopName[0] = "linear. 2 neighbors";
    d.migTopName[1] = "two rows. 3 neighbors";
    d.migTopName[2] = "quadratic. 4 neighbors";
    d.migTopName[3] = "honeycomb. 6 neighbors";
    d.migTopName[4] = "octal. 8 neighbors";
    d.migTopName[5] = "random. from any group";
    d.migTopName[6] = "common gene pool";
    d.migTopName[7] = "dynamic territories";
    d.selModelName[0] = "fecundity selection";
    d.selModelName[1] = "positive viability sel.";
    d.selModelName[2] = "negative viability sel.";
    d.selModelName[3] = "independent viability sel.";
    d.selModelName[4] = "minimum viability sel.";
    d.extinctionPatternName[0] = "depends on own group";
    d.extinctionPatternName[1] = "depends on random neighbor";
    d.extinctionPatternName[2] = "depends on strongest neighbor";
    d.warPatternName[0] = "war against all";
    d.warPatternName[1] = "depends on shared border length";
    d.warPatternName[2] = "depends on difference in strength";

    d.nLoci = 3;     // number of loci
    d.sLocusName[locusAltruism] = "altruism";
    d.sLocusName[locusEndogamy] = "endogamy";
    d.sLocusName[locusConformity] = "conformity";

    d.sGeneName[locusAltruism][0] = "E";
    d.sGeneName[locusAltruism][1] = "A";
    d.sGeneName[locusEndogamy][0] = "Exo";
    d.sGeneName[locusEndogamy][1] = "Endo";
    d.sGeneName[locusConformity][0] = "NonConf";
    d.sGeneName[locusConformity][1] = "Conform";

    d.sPhenotypeName[locusAltruism][0] = "egoist";
    d.sPhenotypeName[locusAltruism][1] = "altruist";
    d.sPhenotypeName[locusEndogamy][0] = "exogamist";
    d.sPhenotypeName[locusEndogamy][1] = "endogamist";
    d.sPhenotypeName[locusConformity][0] = "nonconformist";
    d.sPhenotypeName[locusConformity][1] = "conformist";

    d.fitRowLocus = locusAltruism;
    d.fitColLocus = locusAltruism;
    d.fitSupColLocus = locusConformity;
    d.fitConditionName = 0;

    d.bSelModels = -1;
    d.bIndividualPropertiesUsed = 0xF;
    d.bExtinctionPatterns = -1;
    d.bWarPatterns = 0;
    d.bGroupPropertiesUsed = -1;
    d.bFitFunc = -1;

    if (parameterFilePath.isEmpty()) parameterFilePath = QDir::currentPath();
    if (lastParameterFile.isEmpty()) lastParameterFile = QDir::currentPath() + "/last.altru";    
}

void Altruist::update() {
    // update user interface. This must be done inside GUI thread
    // check global error reporter
    static int numErrors = 0;
    if (errors.getNum()) {
        pause();                                 // pause current run
        errorMessage(errors.getMessage());       // make message box with error messag
        errors.clear();
        if (++numErrors > 4) close();            // stop infinitely repeating error popup messages
    }
    else {
        // draw graphics and show status bar
        char line[200];
        if (d.requestUpdate) {
            if (d.sweepsUsed == 0) {
                setStateTitle();
                sprintf_s(line, "generation %i", (int)d.generations);
                statusBar()->showMessage(line);
                graphicsView->draw();
            }
            else {
                // parameter sweeps are running. show state for all loops
                setStateTitle();

                const char * format = "%X";
                bool written = false;
                line[0] = 0;
                // show outermost loop first
                for (int iloop = maxSweeps - 1; iloop >= 0; iloop--) {
                    if (worker->loops[iloop].state != loopStateUnused) { // loop is used
                        if (written) {
                            snprintf(line + strlen(line), 100, ",  "); // comma after previous parameter
                        }
                        if (worker->loops[iloop].type == loopSearch && worker->loops[iloop].state > loopStateRunning) {
                            // parameter search finished. show results
                            snprintf(line + strlen(line), 100, "limit1 = %.4f, limit2 = %.4f",
                                worker->loops[iloop].limit1, worker->loops[iloop].limit2);                        
                        }
                        else {
                            int index = worker->loops[iloop].parIndex;      // get name of loop parameter
                            snprintf(line + strlen(line), 100, "%s = ", sweepParameterList[index].name);
                            // find appropriate format for value
                            if (sweepParameterList[index].varType <= varInt32) format = "%.0f";
                            else if (worker->loops[iloop].type == 2) format = "%.3f";  // search
                            else if (worker->loops[iloop].logarithmic) format = "%.3f";  // logarithmic
                            else if (fabs(worker->loops[iloop].increment) <= 0.0051f) format = "%.3f";
                            else if (fabs(worker->loops[iloop].increment) <= 0.051f) format = "%.2f";
                            else if (fabs(worker->loops[iloop].increment) <= 0.51f) format = "%.1f";
                            else format = "%.2G";
                            float value;
                            if (worker->lastResult != resultUnknown) value = worker->loops[iloop].lastValue;
                            else value = worker->loops[iloop].parameterValue; // no result yet                               
                            snprintf(line + strlen(line), 100, format, value);
                        }
                        written = true;
                    }
                }
                const char * result;
                switch (worker->lastResult) {    // write result
                case resultAborted:              // max generations reached or aborted
                    result = " Aborted";  break;
                case resultEgoism:               // wild type fixation
                    result = " E";  break;
                case resultPolymorphism:         // polymorphism
                    result = " P";  break;
                case resultAltruism:             // mutant fixation
                    result = " A";  break;
                case resultDied:                 // all died
                    result = " Died";  break;
                default:
                    result = " _";
                }
                snprintf(line + strlen(line), 100, result);

                statusBar()->showMessage(line);
                graphicsView->draw();
            }
        }
    }
    d.requestUpdate = false;
}


// List of supported models
// Each model source file must call models.addModel to register itself

void ModelDescriptorList::addModel(ModelDescriptor const & m) {
    // we cannot make sure that the constructor for the model list is called before
    // the constructors of each model. To solve this problem, we make nModels initialized
    // on first use
    static int n = 0;
    if (n >= maxModels) { // model list is full
        messageBox("Error: Number of models exceeds maxModels");
    }
    else {
        models[n++] = m;
        nModels = n;
        sorted = false;
    }
}

uint32_t ModelDescriptorList::getNum() {
    return nModels;
}

ModelDescriptor * ModelDescriptorList::getModel(uint32_t i) {
    if (i >= nModels) return nullptr;
    else return models+i;
}

// constructor. Data members of this class are deliberately not 
// initialized in the constructor, but in the first call to addModel
ModelDescriptorList::ModelDescriptorList() {
}

void ModelDescriptorList::sort() {
    // Sort models alphabetically. Simple Bubble sort
    int32_t i, j;
    ModelDescriptor temp, * p1, * p2;
    if (!sorted) {
        for (i = 0; i < nModels; i++) {
            for (j = 0; j < nModels - i - 1; j++) {
                p1 = models + j;
                p2 = models + j + 1;
                if (strcasecmp(p1->name, p2->name) > 0) {
                    temp = *p1;  *p1 = *p2;  *p2 = temp; // swap records
                }
            }
        }
        sorted = true;
    }
}

// constructor for each model source file to add the model to the list
Construct::Construct (ModelDescriptor const & m) {
    models.addModel(m);
}
