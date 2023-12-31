/*******************************  menus.h   ***********************************
* Author:        Agner Fog
* Date created:  2023-08-08
* Last modified: 2023-12-31
* Version:       3.00.00
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines menus and dialog boxes in the user interface
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#pragma once

/******************************************************************************
*                              Dialog boxes
******************************************************************************/

// make message box. Can only be called from GUI thread
void messageBox(const char * text, const char * title = " ");

// make message box. Can only be called from GUI thread
void errorMessage(QString text);

// dialog box "select model"
class ModelDialogBox : public QDialog {
public:
    ModelDialogBox(Altruist *parent);
private:
    QGridLayout dialogLayout;
    QComboBox modelComboBox;
    QLabel modelLabel;
    QPushButton okButton;
    QPushButton cancelButton;
    QPushButton descriptionButton;
};

// dialog box "geography and migration"
class GeographyDialogBox : public QDialog {
public:
    GeographyDialogBox(Altruist *parent);
private:
    QGridLayout dialogLayout;
    QLabel labelTotArea;               // total area of habitat (floating territories)
    QLineEdit editTotArea;
    QLabel labelMaxIslands;            // maximum number of islands or demes
    QLineEdit editMaxIslands;
    QLabel labelMaxPerDeme;            // max number of individuals per deme
    QLineEdit editMaxPerDeme;

    QLabel labelPopulationDensity0;    // carrying capacity per area unit under egoism
    QLineEdit editPopulationDensity0;
    QLabel labelPopulationDensity1;    // carrying capacity per area unit under altruism
    QLineEdit editPopulationDensity1;

    QLabel labelTerritorySizeMax;      // max area of deme territory (floating territories)
    QLineEdit editTerritorySizeMax;
    QLabel labelTerritorySizeMin;      // min area of deme territory (floating territories)
    QLineEdit editTerritorySizeMin;
    QLabel labelCapacityStdDev;        //  minimum group size
    QLineEdit editCapacityStdDev;
    QLabel labelRecolGroupSize;        //  recolonization group size
    QLineEdit editRecolGroupSize;

    QLabel labelMix;                   // interdemic migration coefficient (float)
    QLineEdit editMix;
    QLabel labelMix2;                  // interdemic migration coefficient under endogamy (float)
    QLineEdit editMix2;

    QLabel labelEmi_p;                 // emigration pattern 
    QComboBox comboEmi_p;
    QLabel labelImmi_p;                // immigration pattern 
    QComboBox comboImmi_p;
    QLabel labelRecol_p;               // recolonization pattern 
    QComboBox comboRecol_p;
    QLabel labelMigTop;                // migration topology
    QComboBox comboMigTop;

    QPushButton okButton;              // OK button
    QPushButton cancelButton;          // cancel button

protected:
    uint32_t fieldsUsed1;              // which fields are used
};


// dialog box "Loci and mutation"
class LociDialogBox : public QDialog {
public:
    LociDialogBox(Altruist *parent);
private:
    QGridLayout dialogLayout;
    QLabel labelLocusName[maxLoci];              // name of locus
    QLabel labelLocusUsed[maxLoci];              // label for checkbox locus used
    QLabel labelGeneNames[maxLoci];              // label for gene names
    QLabel labelDominance[maxLoci];              // label for mutant gene dominance
    QLabel labelInitialFraction[maxLoci];        // label for initial fraction of mutant gene
    QLabel labelForwardMutation[maxLoci];        // label for forward mutation rate
    QLabel labelBackMutation[maxLoci];           // label for back mutation rate
    QCheckBox checkboxLocusUsed[maxLoci];        // checkbox for locus used
    QComboBox comboDominance[maxLoci];           // combo box for mutant gene dominance
    QLineEdit editInitFraction[maxLoci];         // edit field for initial fraction of mutant gene
    QLineEdit editForwardMutation[maxLoci];      // edit field for forward mutation rate
    QLineEdit editBackMutation[maxLoci];         // edit field for forward mutation rate
    QPushButton okButton;
    QPushButton cancelButton;
    void checkGrey();                            // set disabled fields grey
};

// dialog box "Individual and group fitness"
class FitnessDialogBox : public QDialog {
public:
    FitnessDialogBox(Altruist *parent);
private:
    QGridLayout dialogLayout;
    QLabel labelFitnessFor;
    QLabel labelFitRow[2];
    QLabel labelFitColumn[4];
    QLabel labelFit2[2];
    QLabel labelSelModel;
    QLabel labelGrowthRate;
    QLabel labelLeaderAdvantage;
    QLabel labelLeaderSelection;
    QLabel labelGroupProperties;
    QLabel labelExtinctionPattern;
    QLabel labelWarPattern;
    QLabel labelSurviv;
    QLabel labelwarIntensity;    
    QLabel labelExtinctionRateFor;
    QLabel labelExtinctionRate0;
    QLabel labelExtinctionRate1;
    QLabel labelExtinctionRate2;
    QLabel labelExtinctionRate3;
    QLabel labelHaystackPeriod;
    QLabel labelMixingPeriod;
    QLabel labelDependGroupSuccessOn;
    QLabel labelGroupFitnessFactor[6];
    QLabel labelGroupFitnessFunction;
    QLabel labelGroupFitnessExponent;

    QLineEdit editIndividualFitness[8];     // edit field for individual fitness
    QLineEdit editIndividualFitness2[4];    // edit field for individual fitness at extra loci
    QComboBox comboSelModel;
    QLineEdit editGrowthRate;
    QLineEdit editLeaderAdvantage;
    QLineEdit editLeaderSelection;

    QComboBox comboExtinctionPattern;
    QComboBox comboWarPattern;
    QLineEdit editSurviv;
    QLineEdit editwarIntensity;    
    QLineEdit editHaystackPeriod;
    QLineEdit editMixingPeriod;

    QLineEdit editExtinctionRates[4];
    QLineEdit editGroupFitnessFactors[6];
    QComboBox comboGroupFitnessFunction;
    QLineEdit editGroupFitnessExponent;

    QPushButton okButton;
    QPushButton cancelButton;

    QLabel labelModelSpecificParameters;
    QLabel labelModelSpecific[maxModelSpecificParameters];
    QLineEdit editModelSpecific[maxModelSpecificParameters];
    QCheckBox checkboxModelSpecific[maxModelSpecificParameters];

protected:
    // lists of cathegorical variables in combo box
    int  groupfitFunctions[5];     // index translation for group fitness function
    int  groupfitFunctionsRev[5];  // reverse index translation for group fitness function
    int  groupfitIndex;            // current selection in combo box
    bool groupfitModified;         // group fitness index or exponent modified
};


// dialog box "Run control"
class RunControlDialogBox : public QDialog {
public:
    RunControlDialogBox(Altruist *parent);
private:
    AltruData * d;
    int32_t enableSweepParameter;
    QGridLayout dialogLayout;
    QLabel labelRandomSeed;
    QLabel labelMinGenerations;
    QLabel labelMaxGenerations;
    QLabel labelStopCriterion;
    QLabel labelCriterionDegree;
    QLabel labelDelay;
    QLabel labelParameterSweeps;
    QLabel labelSweepType[maxSweeps];
    QLabel labelParameterToSweep[maxSweeps];
    QLabel labelStartValue[maxSweeps];
    QLabel labelEndValue[maxSweeps];
    QLabel labelStep[maxSweeps];
    QLineEdit editRandomSeed;
    QLineEdit editMinGenerations;
    QLineEdit editMaxGenerations;
    QComboBox comboStopCriterion;
    QLineEdit editDelay;
    QLineEdit editCriterionDegree;
    QComboBox comboSweepType[maxSweeps];
    QComboBox comboParameterToSweep[maxSweeps];
    QLineEdit editStartValue[maxSweeps];
    QLineEdit editEndValue[maxSweeps];
    QLineEdit editStep[maxSweeps];
    QPushButton okButton;
    QPushButton cancelButton;
    void checkGrey();  // set disabled fields grey
};


// dialog box "Data file output"
class DataOutputDialogBox : public QDialog {
public:
    DataOutputDialogBox(Altruist *parent);
private:
    QGridLayout dialogLayout;
    QLabel labelMakeOutputFile;
    QLabel labelTitle;
    QLabel labelOutputInterval;
    QLabel labelSteadyStateAfter;
    QLabel labelInitialParameters;
    QLabel labelGenerations;
    QLabel labelInhabitedIslands;
    QLabel labelGeneFractions;
    QLabel labelAltruists;
    QLabel labelMutations;
    QLabel labelMigrants;
    QLabel labelAltruismDemes;
    QLabel labelExtinctions;
    QLabel labelResult;
    QLabel labelModelSpecific;
    QLabel labelSteadyState;
    QLabel labelTimeConsumption;
    QLabel labelSearchResult;

    QLineEdit editTitle;
    QLineEdit editOutputInterval;
    QLineEdit editSteadyStateAfter;

    QCheckBox checkboxMakeOutputFile;
    QCheckBox checkboxInitialParameters;
    QCheckBox checkboxGenerations;
    QCheckBox checkboxInhabitedIslands;
    QCheckBox checkboxGeneFractions;
    QCheckBox checkboxAltruists;
    QCheckBox checkboxMutations;
    QCheckBox checkboxMigrants;
    QCheckBox checkboxAltruismDemes;
    QCheckBox checkboxExtinctions;
    QCheckBox checkboxResult;
    QCheckBox checkboxModelSpecific;
    QCheckBox checkboxSteadyState;
    QCheckBox checkboxTimeConsumption;
    QCheckBox checkboxSearchResult;
    QPushButton changeButton;
    QPushButton okButton;
    QPushButton cancelButton;
};
