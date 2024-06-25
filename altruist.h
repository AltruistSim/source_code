/****************************  altruist.h   ***********************************
* Author:        Agner Fog
* Date created:  1997-10-09
* Last modified: 2024-06-25
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines common data structures, classes, and functions
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#pragma once
#include "random.h"

#ifdef _MSC_VER 
#define strcasecmp _stricmp
#endif

#define ARR_LEN(x) (sizeof(x)/sizeof(x[0]))      // length of an array

// Altruist version number
const int altruistMajorVersion = 3;
const int altruistMinorVersion = 1;

// maximum number of gene loci
const int maxLoci = 4;

// maximum number of model-specific parameters in dialog box
const int maxModelSpecificParameters = 10;

// maximum number of parameter sweep loops
const int maxSweeps = 3;


/******************************************************************************
*                        named values
******************************************************************************/

// varType values
const int varInt8            = 1;      // 8 bit signed integer
const int varInt16           = 2;      // 16 bit signed integer
const int varInt32           = 3;      // 32 bit signed integer
const int varInt64           = 4;      // 64 bit signed integer
const int varFloat           = 9;      // 32 bit float
const int varDouble          = 10;     // 64 bit double precision float

// runState values
const int state_idle         = 0;      // idle
const int state_start        = 1;      // start simulation
const int state_run          = 2;      // simulation running
const int state_single_step  = 3;      // single stepping one generation at a time
const int state_pause_next   = 4;      // pause when current generation finished
const int state_pause        = 5;      // paused
const int state_stop         = 6;      // stop and abort
const int state_error        = 7;      // error occurred during simulation

// state for model initialization function in each model
const int model_initialize   = 1;      // initialize model
const int model_parameters_changed = 2;// parameters changed

// values for Mendelian dominance of a gene
const int recessive          = 0;      // mutant gene is recessive
const int dominant           = 1;      // mutant gene is dominant
const int incompleteDominant = 2;      // mutant gene is partially dominant

// values for simulationResult
const int resultUnknown      = 0;      // result not available yet
const int resultEgoism       = 1;      // egoism or wild type fixation
const int resultPolymorphism = 2;      // stable polymorphism
const int resultAltruism     = 3;      // altruism or mutant gene fixation
const int resultDied         = 4;      // all died
const int resultAborted      = 5;      // simulation was not finished

// values for stop criterion
const int stopCritNone       = 0;      // none
const int stopIfGenesUniform = 1;      // stop if genes uniform at all loci
const int stopIfPhenoUniform = 2;      // stop if phenotype uniform at all loci
const int stopIfGene0Uniform = 3;      // stop if genes uniform at locus 0
const int stopIfGene0Mutant  = 4;      // stop if mutant fixated at locus 0
const int stopIfGene1Uniform = 5;      // stop if genes uniform at locus 1
const int stopIfGene1Mutant  = 6;      // stop if mutant fixated at locus 1
const int stopIfGene2Uniform = 7;      // stop if genes uniform at locus 2
const int stopIfGene2Mutant  = 8;      // stop if mutant fixated at locus 2

// values for stop cause
const int stopCauseNone      = 0;      // not stoppen
const int stopCauseAbort     = 1;      // aborted
const int stopCauseDied      = 2;      // all died
const int stopCauseMaxGen    = 3;      // maximum number of generations reached
const int stopCauseCriterion = 4;      // stop criterion reached

// Don't change these values unless you want to change all parameter files:
// values for immigrationPattern
const int immigrationCommonPool = 0;   // immigration from a common gene pool with the same composition as the entire metapopulation
const int immigrationNeighbor = 1;     // immigration from one neighbor, randomly selectiond
const int immigrationProportional = 2; // immigration from one neighbor, proportional with own population
const int immigrationVacant = 3;       // immigration from one neighbor, proportional with own vacant capacity
const int immigrationAllNeighbors = 4; // immigration from all neighbors
const int immigrationRandomGroup = 5;  // immigration from a random group
const int immigrationStrongNeighbor = 6; // immigration from neighbor group selected with probability proportional w strength, as decided by emigration pattern

// values for colonizationPattern
const int colonizationCommonPool = immigrationCommonPool; // colonists come from a common gene pool with the same composition as the entire metapopulation
const int colonizationNeighbor = immigrationNeighbor; // colonists all come from the same random neighbor group
// values 2 and 3 unused
const int colonizationAllNeighbors = immigrationAllNeighbors;// colonists come from all neighbor groups
const int colonizationRandomGroup = immigrationRandomGroup;// colonists all come from the same random island, not necessarily a neighbor island
const int colonizationStrongNeighbor = 6; // colonists from neighbor group selected with probability proportional w strength, as decided by emigration pattern

// values for emigrationPattern
const int emigrationConstant = 0;      // emigration independent of population
const int emigrationProportional = 1;  // emigration proportional to population
const int emigrationExcess = 2;        // emigration proportional to excess population
const int emigrationGroupfit = 3;      // emigration proportional to group fitness * population

// values for migrationTopology:
const int topologyLinear = 0;          // linear organization, 2 neighbors
const int topology3neighb = 1;         // two rows, 3 neighbors
const int topologyQuadratic = 2;       // quadratic organization, 4 neighbors
const int topologyHoneycomb = 3;       // honeycomb organization, 6 neighbors
const int topologyOctal = 4;           // octal organization, 8 neighbors
const int topologyRandom = 5;          // random migration, no organization
const int topologyCommonPool = 6;      // migration from a common gene pool with the same composition as the entire metapopulation
const int topologyFloating = 7;        // territories with floating boundaries formed by conflict

// values for extinctionPattern:
const int extinctionOwn = 0;           // extinction probability depend on own group fitness
const int extinctionNeighbor = 1;      // extinction probability depend on fitness relative to random neighbor
const int extinctionStrongest = 2;     // extinction probability depend on fitness relative to strongest neighbor

// values for warPattern:
const int warAgainstAll = 0;           // attack all neighbor groups every generation
const int warAgainstBorder = 1;        // attack one neighbor group. probability by length of shared border
const int warAgainstWeak = 2;          // attack one neighbor group. probability by difference in strenght

// values for selectionModel
const int selectionFecundity = 0;
const int selectionPositiveViability = 1;
const int selectionNegativeViability = 2;
const int selectionIndependentViability = 3;
const int selectionMinimumViability = 4;

// locus types (index into locusUsed[] etc.)
const int locusAltruism = 0;           // altruism locus
const int locusEndogamy = 1;           // endogamy locus
const int locusConformity = 2;         // conformity locus


// Preliminary declarations

struct ModelDescriptor;                // defined below
class  Worker;                         // defined below
class  AltruistView;                   // defined in graphics.h

union UserData {                       // model-specific variables stored in common AltruData
    int32_t i;                         // integer
    float f;                           // floating point
};
    
const int maxUserData = 32;            // maximum number of model-specific variables

const int numExtraBuffers = 5;         // number of buffers for extra model-specific data


/******************************************************************************
*                               model parameter data
******************************************************************************/

struct AltruData {
    int8_t * demeData;                 // buffer for storing demes
    int64_t bufferSize;                // size of allocated data
    // Extra buffers can be used by model modules for allocating data structures for arbitrary use.
    // The pointers must be null when not allocated and non-zero when allocated.
    // The buffers will be deallocated if nonzero by operator delete[] when the program exits.
    // If the model module deallocates a buffer, it must set the pointer to null.
    int8_t * extraBuffer[numExtraBuffers];     // extra buffers for arbitrary use by model modules
    int32_t  extraBufferSize[numExtraBuffers]; // size of extra allocated buffers (in bytes)

    RandomVariates * ran;              // random number generator

    //**************** Parameters: Model *********************
    int32_t iModel;                    // model index
    int modelVersionMajor;             // model version
    int modelVersionMinor;             // model sub-version
    int modelDemeStructureSize;        // size of model-specific deme structure
    int modelPopOffset;                // offset to population count in model-specific deme structure

    //**************** Parameters: Geography and migration *********************
    uint32_t bGeographyParametersUsed; // define which geography parameters are used:
    // 1: totArea
    // 2: maxIslands
    // 4: nMaxPerDeme
    // 0x10: carryingCapacity under egoism
    // 0x20: carryingCapacity under altruism (if different)
    // 0x40: territorySizeMax
    // 0x80: territorySizeMin
    // 0x100: carryingCapacityStandardDeviation
    // 0x200: colonySize
    // 0x400:
    // 0x800:
    // 0x1000: migrationRate[0]
    // 0x2000: migrationRate[1]
    // 0x10000: emigrationPattern
    // 0x20000: (emigrationPattern 2)
    // 0x40000: immigrationPattern
    // 0x80000: (immigrationPattern 2)
    // 0x100000: colonizationPattern
    // 0x200000: migrationTopology

    int32_t maxIslands;                // maximum number of islands or demes
    int32_t nMaxPerDeme;               // max number of individuals per deme
    int32_t totArea;                   // total area of habitat (floating territories)
    int32_t lastArea;                  // value of totArea when simulation started
    int32_t rowLength;                 // length of rows in organization of islands
    int32_t rowLengthTerri;            // length of rows in territorial map
    int32_t numRows;                   // number of rows in territorial map or island map
    float   carryingCapacity[2];       // carrying capacity per island or area unit under egoism/altruism
    float   carryingCapacityStandardDeviation; // standard deviation if carryingCapacity is random
    int32_t territorySizeMax;          // max area of deme territory (floating territories)
    int32_t territorySizeMin;          // min area of deme territory (floating territories)
    int32_t minGroupSize;              // minimum group size
    int32_t colonySize;                // recolonization group size
    float   migrationRate[2];          // interdemic migration coefficient, normal and under endogamy

    int32_t emigrationPattern;         // emigration pattern:
                                       // 0 = constant, 1 = proportional with population, excess production
    int32_t bEmigrationPattern;        // define which emigration patterns supported:
                                       // 1: constant, 2: proportional, 4: excess,

    int32_t immigrationPattern;        // immigration pattern:
                                       // immigrationCommonPool, immigrationNeighbor, immigrationProportional, 
                                       // immigrationVacant, immigrationAllNeighbors, immigrationRandomGroup
    int32_t bImmigrationPattern;       // define which immigration patterns supported: one bit for each pattern

    int32_t colonizationPattern;       // recolonization pattern: 
                                       // colonizationNeighbor, colonizationStrongNeighbor, colonizationAllNeighbors,
                                       // colonizationRandomGroup, colonizationCommonPool
    int32_t bColonPat;                 // define which recolonization patterns supported: one bit for each pattern

    int32_t migrationTopology;         // migration topology:
                                       // topologyLinear, topology3neighb, topologyQuadratic, topologyHoneycomb,
                                       // topologyOctal, topologyRandom, topologyCommonPool, topologyFloating
    uint32_t bTopology;                // define which migration topologies supported: one bit for each pattern

    const char * emiPatName[4];        // names of emigration patterns
    const char * immiPatName[7];       // names of immigration patterns
    const char * colonPatName[7];      // names of colonization patterns
    const char * migTopName[8];        // names of migration topologies  


    //**************** Parameters: Loci and mutation *********************
    int32_t nLoci;                     // maximum number of loci for the chosen model
    const char * sLocusName[maxLoci];  // names of loci
    const char * sGeneName[maxLoci][2];// names of genes
    const char * sPhenotypeName[maxLoci][2]; // names of phenotypes
    bool    locusUsed[maxLoci];        // true if locus is used
    int8_t  dominance[maxLoci];        // mutant gene dominance:
                                       // recessive, dominant, half dominant
    float fg0[maxLoci];                // initial gene fraction for mutant
    float murate[maxLoci][2];          // murate[i][0] = forward mutation rate
                                       // murate[i][1] = back mutation rate


//**************** Parameters: Individual fitness *********************
    int32_t bIndividualPropertiesUsed; // set bits for each parameter used:
                                       // bit 0-7: fit[0-7]
                                       // bit 8-11: fit2[0-3]
                                       // bit 16: growthRate
    float  fit[8];                     // intrademic reproductive fitness of:
                                       // fit[0]: egoists   among egoists
                                       // fit[1]: egoists   among altruists
                                       // fit[2]: altruists among egoists
                                       // fit[3]: altruists among altruists
                                       // fit[4] - fit[7]: same, under conformity
    float  fit2[4];                    // relative fitness of mutants at the second and subsequent loci
    int8_t fitRowLocus;                // locus represented in fit[] rows in fitness dialog box
    int8_t fitColLocus;                // locus represented in fit[] columns in fitness dialog box
    int8_t fitSupColLocus;             // locus represented in fit[] super-columns in fitness dialog box
    int8_t fitConditionName;           // add word to column names in fit[] dialog box: 0: none, 1: "among", 2: "under"

    float growthRate;                  // average fertility
    int32_t selectionModel;            // define which selection model is used:
                                       // selectionFecundity, selectionPositiveViability
                                       // selectionNegativeViability, selectionIndependentViability
                                       // selectionMinimumViability
    const char * selModelName[5];      // names of selection models
                                       // "fecundity selection", "positive viability sel.", "negative viability sel.";
                                       // "independent viability sel.", "minimum viability sel."
    int32_t bSelModels;                // define which selection models supported, one bit for each

    //**************** Parameters: group fitness *********************
    int32_t extinctionPattern;         // extinction pattern:
                                       // extinctionOwn, extinctionNeighbor, extinctionStrongest
    const char * extinctionPatternName[3]; // names of extinction patterns
    int32_t bExtinctionPatterns;       // define which extinction patterns supported, one bit for each

    int32_t warPattern;                // war pattern:
    const char * warPatternName[3];    // names of war patterns:
                                       // warAgainstAll, warAgainstBorder, warAgainstWeak
    int32_t bWarPatterns;              // define which war patterns supported, one bit for each

    int32_t bGroupPropertiesUsed;      // set bits for each parameter used:
                                       // 1: exr
                                       // 2: exr[0-1]
                                       // 4: exr[0+2]
                                       // 6: exr[0-3]
                                       // 0x10: surviv
                                       // 0x20: extinctionPattern
                                       // 0x40: growthRate
                                       // 0x80: warIntensity
                                       // 0x100: haystackPeriod
                                       // 0x200: mixingPeriod
                                       // 0x1000: fitExpo
                                       // 0x10000: leaderAdvantage
                                       // 0x20000: leaderSelection

    float extinctionRate[4];           // extinction rates for demes (island model)
                                       // extinctionRate[0] = small deme with egoists
                                       // extinctionRate[1] = big   deme with egoists
                                       // extinctionRate[2] = small deme with altruists
                                       // extinctionRate[3] = big   deme with altruists

    float surviv;                      // survival degree for extinct islands or for defeated group in war
    float warIntensity;                // intensity of territorial war
    float leaderAdvantage;             // relative advantage of being a leader in regality model
    float leaderSelection;             // preference for leader being a regalists in regality model

    int32_t bFitFunc;                  // define which group fitness functions supported:
                                       // 1: one is enough, 
                                       // 2: convex(0.5), 4: linear, 8: concave(2.0), 
                                       // 0x10: all or nothing,
                                       // 0x100: fitExpo can have any value

    float fitExpo;                     // group fitness exponent
    int32_t fitfunc;                   // group fitness as function of altruists:
                                       // fitfunc:      fitExpo:        group fitness function:
                                       //    0              0           one altruist is enough
                                       //    1             0-1 (0.5)    convex
                                       //    2              1           linear
                                       //    3             > 1 (2)      concave
                                       //    4             inf          all or nothing

    uint32_t haystackPeriod;           // period in haystacks in haystack model
    uint32_t mixingPeriod;             // period outside haystacks in haystack model


    //**************** Model-specific parameters *********************
    union {
        int32_t modelspec_i[maxModelSpecificParameters];
        float modelspec_f[maxModelSpecificParameters];
    };

    //**************** Parameters for Run control *********************
    int32_t demeStructureSize;         // memory size per deme
    int32_t seed;                      // random number seed
    int32_t minimumGenerations;        // minimum number of generations
    int32_t maximumGenerations;        // maximum number of generations
    int32_t stopCriterion;             // stop criterion:
                                       // stopCritNone, stopIfGenesUniform, stopIfPhenoUniform, stopIfGene0Uniform, ...
    int32_t bStopCriterionUsed;        // set bits for each stop criterion supported
    float   stopCriterionDegree;       // degree for stop criterion in btype (for example 0.99)
    int32_t delayms;                   // delay for screen update, milliseconds
    int32_t sweepsUsed;                // set one bit for each parameter loop used
    int32_t sweepType[maxSweeps];      // 0: none, 1: linear, 2: logarithmic, 3: linear search, 4: log search
    int32_t sweepParameter[maxSweeps]; // parameter to sweep, index into sweepParameterList
    float   sweepStartValue[maxSweeps];// start value for parameter sweep
    float   sweepEndValue[maxSweeps];  // end value
    float   sweepStep[maxSweeps];      // step

    //******   statistics output from model generation function      ******
    int32_t nIslands;                  // number of islands or territories
    int32_t mutations[maxLoci][2];     // number of mutations in this generation
    int32_t demesDied;                 // number of demes extinguished in this generation
    int32_t newColonies;               // number of new colonies formed by splitting or by colonizing vacant islands
    int32_t migrantsTot;               // total number of migrants in this generation

    //*****    statistics that can be calculated by run.cpp:    ****
    int64_t genePool[maxLoci][2];      // global gene pool
    int64_t totalPopulation;           // total number of individuals
    int64_t totalPhenotypes[maxLoci*2];// total number of individuals of various fenotypes. Use of this array is model-specific
    int32_t inhabitedDemes;            // number of inhabited islands
    int32_t altruistDemes;             // number of demes with altruism
    //int32_t egoistDemes;               // number of demes with egoism
    float   geneFraction[maxLoci];     // resulting gene fraction for mutant in each locus
    float   fAltru;                    // overall fraction of phenotypic altruists

    //******          statistics sum over all generations           ******
    double   fgsum[maxLoci];           // sum of mutant gene fractions over generations in steady state
    uint32_t nSum;                     // number of items in sums in fgsum
    uint32_t nSteady;                  // number of generations before steady state statistics
    uint32_t sumMutations;             // sum number of mutations
    uint32_t sumExtinctions;           // sum number of extinctions
    int64_t  sumMigrants;              // sum number of migrants

    //******          run model                                        ******
    int32_t volatile runState;         // state of a single simulation: state_idle, state_start, etc.
    int32_t volatile sweepState;       // state of a parameter sweep: state_idle, state_start, etc.
    int32_t stopCause;                 // stop cause:stopCauseNone, stopCauseAbort, stopCauseDied,
                                       // stopCauseMaxGen, stopCauseCriterion
    // statistics
    int64_t generations;               // count number of generations
    int64_t timeUsed;                  // time used for last simulation, in milliseconds
    int64_t timeUsedLoops;             // time used in a series of simulations, in milliseconds
    ModelDescriptor * currentModel;    // descriptor of current model
    int32_t simulationResult;          // resultAborted, resultEgoism, resultPolymorphism, resultAltruism, resultDied
    bool    parametersChanged;         // parameters have been changed
    bool volatile requestUpdate;       // update graphics display or file output

    //******         output control                   ******
    bool    makeOutputFile;            // make data output file
    int32_t graphicsTypeForModel;      // type of graphic display for current model
    int32_t graphicsType;              // type of graphic display, changed if parameter map

    int32_t fileOutInterval;           // file output interval
    int32_t fileOutSteadyState;        // file output steady state after
    int32_t bOutOptions;               // data output options:
                                       // 1: initial parameters
                                       // 2: generations
                                       // 4: inhabited islands
                                       // 8: gene fractions
                                       // 0x10: phenotypic altruists
                                       // 0x20: mutations
                                       // 0x40: migrants
                                       // 0x80:
                                       // 0x100: altruism demes
                                       // 0x200: extinctions
                                       // 0x400: simulation result
                                       // 0x800: model-specific statistics
                                       // 0x1000: steady state
                                       // 0x2000: time consumption
                                       // 0x10000: search results

    QString outputTitle;               // title on data output
    QString outputFilePath;            // output file drive and path
    QString outputFileName;            // output file name

    //******          space for model-specific variables                 ******
    int32_t nUser;                     // number of model-specific variables
    UserData userData[maxUserData];    // model-specific output variables
    int32_t bUserDataType;             // userData type. one bit for each element. 0 means integer, 1 means float
    const char * * userDataNames;      // list of names of userData
};


/******************************************************************************
*                               main class: Altruist
******************************************************************************/
class Altruist : public QMainWindow {
    Q_OBJECT

public:
    Altruist(QWidget *parent = nullptr);
    ~Altruist();

    AltruData d;                       // all simulation parameters and data
    void initializeData();             // initialize parameters to useful values

    // functions for reading and writing edit fields in dialog boxes
    void writeIntField(QLineEdit & field, int value);
    void writeFloatField(QLineEdit & field, float value);
    static int readIntField(QLineEdit & field);
    static double readFloatField(QLineEdit & field);

protected:
    // functions for reading and writing files
    void readParameterFile(QString filename);
    void writeParameterFile(QString filename);

    // functions for user interface
    void setupMenus();                 // set up menus and dialog boxes
    void update();                     // update graphics output

    // run control
    void run();                        // start running model
    void pause();                      // pause run
    void resume();                     // resume after pause
    void singleStep();                 // run single generation
    void stop();                       // stop running
    void checkIfParametersChanged();   // check if parameters changed while running
    void setStateTitle();              // update the main window title to indicate run state

private:
    Ui::AltruistClass ui;
    char textBuffer[64];               // temporary text buffer

    // menus
    QMenu menuFile;
    QAction * load_action;
    QAction * save_action;
    QAction * exit_action;

    QAction * menuModel;
    QMenu     menuParameters;
    QAction * geography_action;
    QAction * loci_action;
    QAction * fitness_action;
    QAction * outputfile_action1;
    QAction * outputfile_action2;
    QAction * run_control_action1;

    QMenu     menuRun;
    QAction * start_action;
    QAction * pause_action;
    QAction * resume_action;
    QAction * singlestep_action;
    QAction * run_control_action2;
    QAction * abort_action;
    QAction * zoom_in_action;
    QAction * zoom_out_action;

    QAction * menuHelp;

public:
    // file names
    QString parameterFilePath;
    QString parameterFileName;
    QString lastParameterFile;         // last.altru: path and filename for saving last parameters when closing
    Worker * worker;                   // object for workerThread

private:
    QThread workerThread;              // thread for doing heavy calculations
    AltruistView * graphicsView;       // contains functions and data for drawing graphics

public slots:
    void resultReadySlot(int state);   // signal received from worker

signals:
    void doWorkSignal(int state);      // activate worker
};


/******************************************************************************
*                   List of parameters in parameter file
******************************************************************************/
struct ParameterDef {
    int32_t type;                      // 0: header, 1: int, 2: float, 
                                       // 3: int used as bool, 4: enumerated model names
    int32_t num;                       // array size if parameter is an array
    int32_t offset;                    // offset into AltruData (this works as a member pointer)
    const char * name;                 // name of parameter
};

extern ParameterDef const parameterDefinitions[];// list of parameter definitions
extern uint32_t const parameterDefinitionsSize;  // size of parameterDefinitions

// Get offset of a data member of AltruData. 
// This works like a member pointer for a data member of arbitrary type
#define altruDataOffset(x) (int32_t((int8_t*)&((AltruData*)0)->x - (int8_t*)0))

// Get offset of a data member of Deme structure or arbitrary structure. 
#define demeFieldOffset(D,x) (int32_t((int8_t*)&((D*)0)->x - (int8_t*)0))



/******************************************************************************
*                        List of implemented models
******************************************************************************/

const int maxModels = 32; // maximum number of defined models

// Descriptor of Deme structure.
// A deme is a geographically limited territory containing a group of organisms.
// An array of Deme data structures contains gene counts and group proporties for 
// each Deme. The Deme data structure may be different for each simulation model.
// The Deme field descriptor gives information about each field in this structure 
// to the main program:
struct DemeFieldDescriptor {
    int32_t type;            // 0: end of list, 
                             // 1: size of Deme structure, 
                             // 2: carrying capacity (max individuals)
                             // 3: population (number of individuals)
                             // 4: gene count, 
                             // 5: genotype or phenotype count
                             // 8: group property
    int32_t varType;         // 1: 8-bit integer, 2: 16-bit integer, 3: 32-bit integer, 4: 64-bit integer
                             // 8: float, 9: double
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


// Model descriptor in list of models:
// Each simulation model is defined in its own .cpp file.
// Each simulation model has a model descriptor containing information needed 
// by the main program
struct ModelDescriptor {
    const char * name;                                     // name of model
    const char * description;                              // describing text
    const ParameterDef * extraParameters;                  // extra fields in parameters dialog box
    const DemeFieldDescriptor * demeFields;                // list of fields in Deme structure
    void (*initFunction)(AltruData * a, int state);        // initialization function pointer
    void (*generationFunction)(AltruData * a, int mode);   // generation function pointer
};

// Make list of supported models
class ModelDescriptorList {
public:
    ModelDescriptorList();
    void addModel(ModelDescriptor const & m);
    uint32_t getNum();
    ModelDescriptor * getModel(uint32_t i);
    void sort();
protected:
    // These variables are deliberately left uninitialized because the initialization order for
    // global constructors is uncertain. They are initialized on first call to addModel instead:
    uint32_t nModels;
    bool sorted;
    ModelDescriptor models[maxModels];
};

extern ModelDescriptorList models;

// An object of class Construct in each model source file makes sure the model is added to the
// models list
class Construct {
public:
    Construct(ModelDescriptor const & m);
};


// Global error message handler.
// Message boxes can only be made in the GUI thread.
// This class stores error messages to be reported from the GUI thread
class ErrorReporter {
public:
    ErrorReporter() {  // constructor
        clear();
    }
    void reportError(const char * text);         // report an error
    int getNum() {                               // get number of errors
        return num;
    }
    const char * getMessage() {                  // get first error message
        return errorText;
    }
    void clear();                                // clear error message
protected:
    int num;
    char errorText[1024];
};

extern ErrorReporter errors;                     // global ErrorReporter object



/***********************************************************************
  Functions applied to simulation of evolution in file run.cpp
***********************************************************************/

// initialize statics variables before simulation start
void statisticsInit0(AltruData * d);

// initialize statics variables before each generation
void statisticsInit1(AltruData * d);

// update statistics after each generation
void statisticsUpdate(AltruData * d);

// check if stop criterion reached
void checkStopCriterion(AltruData * d);

// stochastic combination of a gene pool into genotypes for a biallelic locus
void combineGenes(int32_t genePool[2], int32_t genotypes[3], RandomVariates * ran);

// differential growth under fecundity selection
void differentialGrowth(int32_t genePool[2], double growthRate[3], int32_t genotypes[3], RandomVariates * ran);

// growth and selection when number of offspring is known
void growthAndSelection(int32_t genes[2], int32_t numChildren, int dominance, float relativeFitness, RandomVariates * ran);

// find all neighbors to a deme, according to migration topology
int findNeighbors(AltruData * d, int32_t deme, int32_t neighbors[]);



/******************************************************************************
*                        Threads and synchronization
******************************************************************************/

// The worker class is doing the heavy calculations in a separate thread
// in order to keep the user interface responsive
class Worker : public QObject {
    Q_OBJECT
public: 
    // constructor
    Worker(AltruData * d) : QObject() , outputFile() {
        this->d = d;
        d->ran = &randomGenerator;
        loopsPaused = false;
        lastResult = 0;
        limitListNum = 0;
        lastLimitListNum = 0;
        resultSets = 0;
    }
    AltruData * d;                               // pointer to data structure with all parameters etc.
    void doWholeSimulation();                    // do a simulation run as part of a parameter sweep or search
    void doMultipleSimulations();                // a parameter sweep or search involvinb multiple simulation runs
    void waitForScreen();                        // wait for curent graphic screen to be updated before changing islands etc.

public slots:
    void doWorkSlot(int state);                  // time-consuming calculations initiated here

signals:
    void resultReadySignal(int state);           // simulation results are ready

protected:
    RandomVariates randomGenerator;              // random number generator
    QElapsedTimer timer;                         // measure time
    // Functions for data file output. This is done in the worker thread to make sure data do not change before they have been written
    void fileOutStart();                         // start writing data output file
    void fileOutGeneration();                    // write single generation to data output file
    void fileOutSimulationFinished();            // write result of a single simulation run (generation loop)
    void fileOutSweepNext();                     // parameter sweep or search
    void fileOutSweepFinished();                 // parameter sweep or search finished. write statistics
    void fileOutModelSpecificResults();          // write model specific results
    QFile outputFile;                            // data output file
public:
    int  resultSets;                             // number of result sets in output file
    // parameter sweeps
    bool loopsPaused;                            // loops are interrupted, may resume
    int loopsFinished;                           // one bit for each loop in state loopStateFinished
    int lastResult;                              // result of last finished simulation
    int limitListNum;                            // number of records in limitList
    int lastLimitListNum;                        // last value of limitListNum drawn in graphics output. chack if map neds to be updated
    float lastParameter[maxSweeps];              // parameters is last finished parameter loop
    ParameterLoop loops[maxSweeps];              // parameter loops defined in parameterloop.cpp
    ResultSet resultList[resultListLength];      // list of simulation results used in parameter search
    ParameterLimits limitList[limitListLength];  // list of search results. contains x limits for each y value
};
