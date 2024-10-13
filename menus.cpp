/****************************  menus.cpp   ************************************
* Author:        Agner Fog
* Date created:  2023-08-07
* Last modified: 2024-10-13
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines the general user interface with menus and dialogs.
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

// general function to show a message box, mainly for debugging purposes
// Can only be called from GUI thread
void messageBox(const char * text, const char * title) {
    QMessageBox mb(title, text, QMessageBox::NoIcon,
        QMessageBox::Ok | QMessageBox::Default,
        QMessageBox::NoButton, QMessageBox::NoButton);
    mb.exec();
}

// general function to show a message box with an error message
// Can only be called from GUI thread
void errorMessage(QString text) {
    QMessageBox mb("Error", text, QMessageBox::NoIcon,
        QMessageBox::Ok | QMessageBox::Default,
        QMessageBox::NoButton, QMessageBox::NoButton);
    mb.exec();
}


void Altruist::setupMenus() {
    // set up menus and dialog boxes
    // menuFile:
    load_action = menuFile.addAction("&Load parameter file");
    save_action = menuFile.addAction("&Save parameter file");
    menuFile.addSeparator() ;
    outputfile_action1 = menuFile.addAction("Set &output file");
    menuFile.addSeparator() ;
    exit_action = menuFile.addAction("E&xit");

    // menuParameters:
    geography_action = menuParameters.addAction("&Geography and migration");
    loci_action = menuParameters.addAction("&Loci and mutation");
    fitness_action = menuParameters.addAction("&Individual and group fitness");
    outputfile_action2 = menuParameters.addAction("&File output");
    run_control_action1 = menuParameters.addAction("Run &control");

    // menuRun:
    run_control_action2 = menuRun.addAction("Run &control");
    menuRun.addSeparator();
    start_action = menuRun.addAction("&Start\tF5");
    pause_action = menuRun.addAction("&Pause\tF6");
    resume_action = menuRun.addAction("&Resume\tF7");
    singlestep_action = menuRun.addAction("Single step\tF8");
    abort_action = menuRun.addAction("&Stop\tF10");
    menuRun.addSeparator();
    zoom_in_action = menuRun.addAction("&Zoom in\tCtrl +");
    zoom_out_action = menuRun.addAction("Zoom &out\tCtrl -");
    
    start_action->setShortcut(Qt::Key_F5);
    pause_action->setShortcut(Qt::Key_F6);
    resume_action->setShortcut(Qt::Key_F7);
    singlestep_action->setShortcut(Qt::Key_F8);
    abort_action->setShortcut(Qt::Key_F10);

    zoom_in_action->setShortcut(Qt::ControlModifier | Qt::Key_Plus);
    zoom_out_action->setShortcut(Qt::ControlModifier | Qt::Key_Minus);

    // menuBar() is defined in parent QMainWindow
    menuBar()->addMenu(&menuFile);
    menuModel = menuBar()->addAction("&Model");
    menuBar()->addMenu(&menuParameters);
    menuBar()->addMenu(&menuRun);
    menuHelp = menuBar()->addAction("&Help");

    // slot connections
    connect(load_action, &QAction::triggered,
        [this](bool checked) {
            QString filename = QFileDialog::getOpenFileName(this, 
                "Read parameter file", parameterFilePath, "*.altru");
            if (filename != "") {
                QFileInfo fi(filename);
                parameterFileName = fi.fileName();
                parameterFilePath = fi.absolutePath();
                readParameterFile(filename);
                d.parametersChanged = true;
                (*(models.getModel(d.iModel)->initFunction))(&d, model_initialize); // initialize new model
            }
        });

    connect(save_action, &QAction::triggered,
        [this](bool checked) {
            QString title = "Save parameter file";
            if (parameterFileName != "") {
                title = title + " (" + parameterFileName + ")"; // show last filename in title
            }
            QString filename = QFileDialog::getSaveFileName(this, 
                title, parameterFilePath, "*.altru");
            if (filename != "") {
                QFileInfo fi(filename);
                parameterFileName = fi.fileName();
                parameterFilePath = fi.absolutePath();
                writeParameterFile(filename);
            }
        });

    connect(outputfile_action1, &QAction::triggered,
        [this](bool checked) {
            auto m = DataOutputDialogBox(this);
            m.show();  m.exec();
        });

    connect(outputfile_action2, &QAction::triggered,
        [this](bool checked) {
            auto m = DataOutputDialogBox(this);
            m.show();  m.exec();
        });

    connect(exit_action, &QAction::triggered,
        [this](bool checked) {
          this->close();
        });

    connect(menuModel, &QAction::triggered,
        [this](bool checked) {
            auto m = ModelDialogBox(this);
            m.show();  m.exec();
        });

    connect(geography_action, &QAction::triggered,
        [this](bool checked) {
            auto m = GeographyDialogBox(this);
            m.show();  m.exec();
            //d.parametersChanged = true;
        });

    connect(loci_action, &QAction::triggered,
        [this](bool checked) {
            auto m = LociDialogBox(this);
            m.show();  m.exec();
            //d.parametersChanged = true;
            (*(models.getModel(d.iModel)->initFunction))(&d, model_initialize); // initialize modified model
        });

    connect(fitness_action, &QAction::triggered,
        [this](bool checked) {
            auto m = FitnessDialogBox(this);
            m.show();  m.exec();
            //d.parametersChanged = true;
        });

    connect(run_control_action1, &QAction::triggered,
        [this](bool checked) {
            auto m = RunControlDialogBox(this);
            m.show();  m.exec();
        });


    connect(start_action, &QAction::triggered,
        [this](bool checked) {
            run();
        });

    connect(pause_action, &QAction::triggered,
        [this](bool checked) {
            pause();
        });
    connect(resume_action, &QAction::triggered,
        [this](bool checked) {
            resume();
        });
    connect(singlestep_action, &QAction::triggered,
        [this](bool checked) {
            singleStep();
        });
    connect(abort_action, &QAction::triggered,
        [this](bool checked) {
            stop();
        });

    connect(run_control_action2, &QAction::triggered,
        [this](bool checked) {
            auto m = RunControlDialogBox(this);
            m.show();  m.exec();
        });

    connect(zoom_in_action, &QAction::triggered,
        [this](bool checked) {
            graphicsView->zoomIn();
        });

    connect(zoom_out_action, &QAction::triggered,
        [this](bool checked) {
            graphicsView->zoomOut();
        });

    connect(menuHelp, &QAction::triggered,
        [this](bool checked) {
            QString helpText("Please see the <a href=\"https://github.com/AltruistSim/documentation/blob/main/altruist_manual.pdf\">Altruist manual</a> for help.");
            QString versionText("<br /><br />Altruist version %1.%2.");
            versionText = versionText.arg(altruistMajorVersion).arg(altruistMinorVersion, 3, 10, QChar('0'));
            //versionText = versionText.arg(altruistMajorVersion).arg(altruistMinorVersion);

            QMessageBox mb("Help", helpText + versionText, QMessageBox::NoIcon,
            QMessageBox::Ok | QMessageBox::Default,
            QMessageBox::NoButton, QMessageBox::NoButton);
            mb.exec();



//            messageBox("Please see altruist_manual.pdf for instructions", "Help");
//            QMessageBox hh(0, QString("Help"), text);
        });
}


// functions for reading and writing edit fields in dialog boxes
void  Altruist::writeIntField(QLineEdit & field, int value) {
    snprintf(textBuffer, sizeof(textBuffer), "%i", value);
    field.setText(textBuffer);    
}

void Altruist::writeFloatField(QLineEdit & field, float value) {
    snprintf(textBuffer, sizeof(textBuffer), "%.6G", value);
    field.setText(textBuffer);    
}

int Altruist::readIntField(QLineEdit & field) {
    return atoi(field.text().toUtf8());
}

double Altruist::readFloatField(QLineEdit & field) {
    return atof(field.text().toUtf8());
}

int expandIndex(int reducedIndex, uint32_t bEnabled) {
    // translate reduced index to full index when some options disabled
    int skip = 0;
    for (int i=0; i<32; i++) {
        if ((bEnabled & 1 << i) == 0) skip++;
        if (i == reducedIndex + skip) return i;
    }
    return 0; // error, not found
}

int compressIndex(int fullIndex, uint32_t bEnabled) {
    // translate full index to reduced index when some options disabled
    int skip = 0;
    for (int i=0; i<32; i++) {
        if ((bEnabled & 1 << i) == 0) skip++;
        if (i == fullIndex) return i - skip;
    }
    return 0; // error, not found
}


/*************************************************************************
*                          make dialog boxes
*************************************************************************/

ModelDialogBox::ModelDialogBox(Altruist *parent) :
    modelLabel(tr("&Model:")),
    okButton(tr("OK"), this),
    cancelButton(tr("Cancel"), this),
    descriptionButton(tr("See description"), this)
{
    models.sort();

    // combo box with list of model names
    for (int i = 0; i < models.getNum(); i++) {
        auto m = models.getModel(i);
        modelComboBox.addItem(tr(m->name), i);
    }
    modelComboBox.setCurrentIndex(parent->d.iModel);
    modelLabel.setBuddy(&modelComboBox);

    dialogLayout.addWidget(&modelLabel, 0, 0);
    dialogLayout.addWidget(&modelComboBox, 0, 1, Qt::AlignRight);
    dialogLayout.addWidget(&descriptionButton, 1, 1);
    modelFileName.setText(parent->parameterFileName); // name of last model file
    dialogLayout.addWidget(&modelFileName, 2, 0);    

    dialogLayout.addWidget(&okButton, 3, 0);
    dialogLayout.addWidget(&cancelButton, 3, 1);
    setLayout(&dialogLayout);

    setWindowTitle(tr("Select model"));

    connect(&descriptionButton, &QPushButton::clicked, [this]() {
        // description button. show description of chosen model
        int i = modelComboBox.currentIndex();
        const char * text = models.getModel(i)->description;
        const char * title = models.getModel(i)->name;
        messageBox(text, title);
        });

    connect(&okButton, &QPushButton::clicked, [this, parent]() {
        // OK button. save model
        parent->d.iModel = modelComboBox.currentIndex();
        close();
        parent->d.parametersChanged = true;
        //parent->initializeData();            // this will reset all parameters when changing model
        (*(models.getModel(parent->d.iModel)->initFunction))(&parent->d, model_initialize); // initialize new model
        });

    connect(&cancelButton, &QPushButton::clicked, [this]() {
        // Cancel button. close without saving
        close();
        });
}

GeographyDialogBox::GeographyDialogBox(Altruist *parent) :
    labelTotArea(tr("&Total area:")),
    labelMaxIslands("Maximum number of groups"),

    labelMaxPerGroup("Maximum individuals per group"),
    labelPopulationDensity0((parent->d.bGeographyParametersUsed & 0x20) ?        
        "Carrying capacity under egoism" : "Carrying capacity"), // carrying capacity per area unit (floating territories)
    labelPopulationDensity1("Carrying capacity under altruism"), // carrying capacity per area unit (floating territories)
    labelTerritorySizeMax("Max. territory area"),                // max area of group territory (floating territories)
    labelTerritorySizeMin("Min. territory area"),                // min area of group territory (floating territories)
    labelCapacityStdDev("Carrying capacity std. deviation"),     // Carrying capacity standard deviation
    labelRecolGroupSize("Colonization group size"),              //  recolonization group size
    labelMix("Migration rate"),
    labelMix2("Migration rate under endogamy"),

    labelEmi_p("Emigration pattern"),
    labelImmi_p("Immigration pattern"),
    labelRecol_p("Colonization pattern"),
    labelMigTop("Migration topology"),

    editTotArea(this),
    editMaxIslands(this),
    editMaxPerGroup(this),
    editPopulationDensity0(this),
    editPopulationDensity1(this),
    editTerritorySizeMax(this),
    editTerritorySizeMin(this),
    editCapacityStdDev(this),
    editRecolGroupSize(this),
    comboEmi_p(this),
    comboImmi_p(this),
    comboRecol_p(this),
    comboMigTop(this),

    okButton("OK", this),
    cancelButton("Cancel", this)
{
    fieldsUsed1 = parent->d.bGeographyParametersUsed;
    //fieldsUsed2 = parent->d.bMigrationParametersUsed;
    int row = 1;      // row in dialogLayout grid
    int i, j; // temporary indexes

    if (fieldsUsed1 & 1) {
        labelTotArea.setBuddy(&editTotArea);
        //editTotArea.setInputMask("9000000000");  // allow 1-10 digits 0-9
        dialogLayout.addWidget(&labelTotArea, row, 0);
        dialogLayout.addWidget(&editTotArea, row++, 1, Qt::AlignRight);
        parent->writeIntField(editTotArea, parent->d.totArea);
    }
    else {
        labelTotArea.setVisible(false);
        editTotArea.setVisible(false);
    }
    if (fieldsUsed1 & 2) {
        labelMaxIslands.setBuddy(&editMaxIslands);
        dialogLayout.addWidget(&labelMaxIslands, row, 0);
        dialogLayout.addWidget(&editMaxIslands, row++, 1, Qt::AlignRight);
        parent->writeIntField(editMaxIslands, parent->d.maxIslands);
    }
    else {
        labelMaxIslands.setVisible(false);
        editMaxIslands.setVisible(false);
    }
    if (fieldsUsed1 & 4) {
        labelMaxPerGroup.setBuddy(&editMaxPerGroup);
        dialogLayout.addWidget(&labelMaxPerGroup, row, 0);
        dialogLayout.addWidget(&editMaxPerGroup, row++, 1, Qt::AlignRight);
        parent->writeIntField(editMaxPerGroup, parent->d.nMaxPerGroup);
    }
    else {
        labelMaxPerGroup.setVisible(false);
        editMaxPerGroup.setVisible(false);
    }
    if (fieldsUsed1 & 0x10) {
        labelPopulationDensity0.setBuddy(&editPopulationDensity0);  // float: allow 0-9+-.E
        dialogLayout.addWidget(&labelPopulationDensity0, row, 0);
        dialogLayout.addWidget(&editPopulationDensity0, row, 1, Qt::AlignRight);
        parent->writeFloatField(editPopulationDensity0, parent->d.carryingCapacity[0]);
    }
    else {
        labelPopulationDensity0.setVisible(false);
        editPopulationDensity0.setVisible(false);
    }
    if (fieldsUsed1 & 0x20) {
        labelPopulationDensity1.setBuddy(&editPopulationDensity1);  // float: allow 0-9+-.E
        dialogLayout.addWidget(&labelPopulationDensity1, row, 2);
        dialogLayout.addWidget(&editPopulationDensity1, row, 3/*, Qt::AlignRight*/);
        parent->writeFloatField(editPopulationDensity1, parent->d.carryingCapacity[1]);
    }
    else {
        labelPopulationDensity1.setVisible(false);
        editPopulationDensity1.setVisible(false);
    }
    if (fieldsUsed1 & 0x30) row++;
    if (fieldsUsed1 & 0x40) {
        labelTerritorySizeMax.setBuddy(&editTerritorySizeMax);
        dialogLayout.addWidget(&labelTerritorySizeMax, row, 0);
        dialogLayout.addWidget(&editTerritorySizeMax, row, 1, Qt::AlignRight);
        parent->writeIntField(editTerritorySizeMax, parent->d.territorySizeMax);
    }
    else {
        labelTerritorySizeMax.setVisible(false);
        editTerritorySizeMax.setVisible(false);
    }
    if (fieldsUsed1 & 0x80) {
        labelTerritorySizeMin.setBuddy(&editTerritorySizeMin);
        dialogLayout.addWidget(&labelTerritorySizeMin, row, 2);
        dialogLayout.addWidget(&editTerritorySizeMin, row++, 3, Qt::AlignRight);
        parent->writeIntField(editTerritorySizeMin, parent->d.territorySizeMin);
    }
    else {
        labelTerritorySizeMin.setVisible(false);
        editTerritorySizeMin.setVisible(false);
    }
    if (fieldsUsed1 & 0x100) {
        labelCapacityStdDev.setBuddy(&editCapacityStdDev);
        dialogLayout.addWidget(&labelCapacityStdDev, row, 0);
        dialogLayout.addWidget(&editCapacityStdDev, row++, 1, Qt::AlignRight);
        parent->writeFloatField(editCapacityStdDev, parent->d.carryingCapacityStandardDeviation);
    }
    else {
        labelCapacityStdDev.setVisible(false);
        editCapacityStdDev.setVisible(false);
    }
    if (fieldsUsed1 & 0x200) {
        labelRecolGroupSize.setBuddy(&editRecolGroupSize);
        dialogLayout.addWidget(&labelRecolGroupSize, row, 0);
        dialogLayout.addWidget(&editRecolGroupSize, row++, 1, Qt::AlignRight);
        parent->writeIntField(editRecolGroupSize, parent->d.colonySize);
    }
    else {
        labelRecolGroupSize.setVisible(false);
        editRecolGroupSize.setVisible(false);
    }

    if (fieldsUsed1 & 0x1000) {
        labelMix.setBuddy(&editMix);
        dialogLayout.addWidget(&labelMix, row, 0);
        dialogLayout.addWidget(&editMix, row++, 1, Qt::AlignRight);
        parent->writeFloatField(editMix, parent->d.migrationRate[0]);
    }
    else {
        labelMix.setVisible(false);
        editMix.setVisible(false);
    }
    if (fieldsUsed1 & 0x2000) {
        labelMix2.setBuddy(&editMix2);
        dialogLayout.addWidget(&labelMix2, row-1, 2);
        dialogLayout.addWidget(&editMix2, row-1, 3, Qt::AlignRight);
        parent->writeFloatField(editMix2, parent->d.migrationRate[1]);
    }
    else {
        labelMix2.setVisible(false);
        editMix2.setVisible(false);
    }

    if (fieldsUsed1 & 0x10000) {
        for (i = 0, j = 0; i < ARR_LEN(parent->d.emiPatName); i++) {
            if (parent->d.bEmigrationPattern & 1 << i) {  // pattern (1 << i) supported
                comboEmi_p.addItem(parent->d.emiPatName[i], j);
            }
        }
        comboEmi_p.setCurrentIndex(compressIndex(parent->d.emigrationPattern, parent->d.bEmigrationPattern));
        labelEmi_p.setBuddy(&comboEmi_p);
        dialogLayout.addWidget(&labelEmi_p, row, 0);
        dialogLayout.addWidget(&comboEmi_p, row++, 1);
    }
    else {
        labelEmi_p.setVisible(false);
        comboEmi_p.setVisible(false);
    }

    if (fieldsUsed1 & 0x40000) {
        for (i = 0, j = 0; i < ARR_LEN(parent->d.immiPatName); i++) {
            if (parent->d.bImmigrationPattern & 1 << i) {  // pattern (1 << i) supported
                comboImmi_p.addItem(parent->d.immiPatName[i], j);
            }
        }
        comboImmi_p.setCurrentIndex(compressIndex(parent->d.immigrationPattern, parent->d.bImmigrationPattern));
        labelImmi_p.setBuddy(&comboImmi_p);
        dialogLayout.addWidget(&labelImmi_p, row, 0);
        dialogLayout.addWidget(&comboImmi_p, row++, 1/*, Qt::AlignRight*/);
    }
    else {
        labelImmi_p.setVisible(false);
        comboImmi_p.setVisible(false);
    }

    if (fieldsUsed1 & 0x100000) {
        for (i = 0, j = 0; i < ARR_LEN(parent->d.colonPatName); i++) {
            if (parent->d.bColonPat & 1 << i) {  // pattern (1 << i) supported
                comboRecol_p.addItem(parent->d.colonPatName[i], j);
            }
        }
        comboRecol_p.setCurrentIndex(compressIndex(parent->d.colonizationPattern, parent->d.bColonPat));
        labelRecol_p.setBuddy(&comboRecol_p);
        dialogLayout.addWidget(&labelRecol_p, row, 0);
        dialogLayout.addWidget(&comboRecol_p, row++, 1/*, Qt::AlignRight*/);
    }
    else {
        labelRecol_p.setVisible(false);
        comboRecol_p.setVisible(false);
    }

    if (fieldsUsed1 & 0x200000) {
        for (i = 0, j = 0; i < ARR_LEN(parent->d.migTopName); i++) {
            if (parent->d.bTopology & 1 << i) {  // pattern (1 << i) supported
                comboMigTop.addItem(parent->d.migTopName[i], j++);                
            }
        }
        comboMigTop.setCurrentIndex(compressIndex(parent->d.migrationTopology, parent->d.bTopology));
        labelMigTop.setBuddy(&comboMigTop);
        dialogLayout.addWidget(&labelMigTop, row, 0);
        dialogLayout.addWidget(&comboMigTop, row++, 1/*, Qt::AlignRight*/);
    }
    else {
        labelMigTop.setVisible(false);
        comboMigTop.setVisible(false);
    }

    dialogLayout.addWidget(&okButton, row, 0);
    dialogLayout.addWidget(&cancelButton, row, 1);
    setLayout(&dialogLayout);

    setWindowTitle(tr("Geography and migration"));

    connect(&cancelButton, &QPushButton::clicked, [this]() {
        // Cancel button. close without saving
        close();
        });

    connect(&okButton, &QPushButton::clicked, [this, parent]() {
        // OK button. save parameters

        if ((fieldsUsed1 & 1) && editTotArea.isModified()) {
            parent->d.totArea = parent->readIntField(editTotArea);
        }
        if ((fieldsUsed1 & 2) && editMaxIslands.isModified()) {
            parent->d.maxIslands = parent->readIntField(editMaxIslands);
        }
        if ((fieldsUsed1 & 4) && editMaxPerGroup.isModified()) {
            parent->d.nMaxPerGroup = parent->readIntField(editMaxPerGroup);
        }
        if ((fieldsUsed1 & 0x10) && editPopulationDensity0.isModified()) {
            parent->d.carryingCapacity[0] = parent->readFloatField(editPopulationDensity0);
        }
        if ((fieldsUsed1 & 0x20) && editPopulationDensity1.isModified()) {
            parent->d.carryingCapacity[1] = parent->readFloatField(editPopulationDensity1);
        }
        if ((fieldsUsed1 & 0x40) && editTerritorySizeMax.isModified()) {
            parent->d.territorySizeMax = parent->readIntField(editTerritorySizeMax);
        }
        if ((fieldsUsed1 & 0x80) && editTerritorySizeMin.isModified()) {
            parent->d.territorySizeMin = parent->readIntField(editTerritorySizeMin);
        }
        if ((fieldsUsed1 & 0x100) && editCapacityStdDev.isModified()) {
            parent->d.carryingCapacityStandardDeviation = parent->readFloatField(editCapacityStdDev);
        }
        if ((fieldsUsed1 & 0x200) && editRecolGroupSize.isModified()) {
            parent->d.colonySize = parent->readIntField(editRecolGroupSize);
        }
        if ((fieldsUsed1 & 0x1000) && editMix.isModified()) {
            parent->d.migrationRate[0] = parent->readFloatField(editMix);
        }
        if ((fieldsUsed1 & 0x2000) && editMix2.isModified()) {
            parent->d.migrationRate[1] = parent->readFloatField(editMix2);
        }
        if (fieldsUsed1 & 0x10000) {parent->d.emigrationPattern = expandIndex(comboEmi_p.currentIndex(), parent->d.bEmigrationPattern);}

        if (fieldsUsed1 & 0x40000) {
            parent->d.immigrationPattern = expandIndex(comboImmi_p.currentIndex(), parent->d.bImmigrationPattern);
        }
        if (fieldsUsed1 & 0x100000) {
            parent->d.colonizationPattern = expandIndex(comboRecol_p.currentIndex(), parent->d.bColonPat);
        }
        if (fieldsUsed1 & 0x200000) {
            parent->d.migrationTopology = expandIndex(comboMigTop.currentIndex(), parent->d.bTopology);
        }

        close();
    });
}

LociDialogBox::LociDialogBox(Altruist *parent) :
    okButton("OK", this),
    cancelButton("Cancel", this)
{
    int i;         // locus index
    int row = 0;   // layout row

    // initialize locus labels and fields
    for (i = 0; i < maxLoci; i++) {
        if (i < parent->d.nLoci) {
            QString locusName = QString("Locus ") + parent->d.sLocusName[i];
            labelLocusName[i].setText(locusName);
            dialogLayout.addWidget(&labelLocusName[i], row, 0);
            labelLocusUsed[i].setText("Locus used:");
            dialogLayout.addWidget(&labelLocusUsed[i], row, 2);
            QString geneNames = QString("Wild gene: ") + parent->d.sGeneName[i][0] + ", mutant: " + parent->d.sGeneName[i][1];
            labelGeneNames[i].setText(geneNames);
            dialogLayout.addWidget(&labelGeneNames[i], row, 1);            
            QString geneDominant = QString("Mutant gene ") + parent->d.sGeneName[i][1] + " dominance:";
            labelDominance[i].setText(geneDominant);
            dialogLayout.addWidget(&labelDominance[i], row+1, 0);
            labelInitialFraction[i].setText("Initial fraction: ");
            dialogLayout.addWidget(&labelInitialFraction[i], row+1, 2);
            labelForwardMutation[i].setText("Forward mutation rate: ");
            dialogLayout.addWidget(&labelForwardMutation[i], row+2, 0);
            labelBackMutation[i].setText("Back mutation rate: ");
            dialogLayout.addWidget(&labelBackMutation[i], row+2, 2);

            dialogLayout.addWidget(&checkboxLocusUsed[i], row, 3);
            labelLocusUsed[i].setBuddy(&checkboxLocusUsed[i]);
            checkboxLocusUsed[i].setChecked(parent->d.locusUsed[i]);
            comboDominance[i].addItem("recessive");
            comboDominance[i].addItem("dominant");
            comboDominance[i].addItem("incomplete");
            comboDominance[i].setCurrentIndex(parent->d.dominance[i]);
            dialogLayout.addWidget(&comboDominance[i], row+1, 1);
            labelDominance[i].setBuddy(&comboDominance[i]);
            dialogLayout.addWidget(&editInitFraction[i], row+1, 3);
            labelInitialFraction[i].setBuddy(&editInitFraction[i]);
            parent->writeFloatField(editInitFraction[i], parent->d.fg0[i]);
            dialogLayout.addWidget(&editForwardMutation[i], row+2, 1);
            labelForwardMutation[i].setBuddy(&editForwardMutation[i]);
            parent->writeFloatField(editForwardMutation[i], parent->d.murate[i][0]);
            dialogLayout.addWidget(&editBackMutation[i], row+2, 3);
            labelBackMutation[i].setBuddy(&editBackMutation[i]);
            parent->writeFloatField(editBackMutation[i], parent->d.murate[i][1]);
            dialogLayout.setRowMinimumHeight(row+3, 20);
            row += 4;
        }
        else {
            labelLocusName[i].setVisible(false);
            labelLocusUsed[i].setVisible(false);
            labelGeneNames[i].setVisible(false);
            labelDominance[i].setVisible(false);
            labelInitialFraction[i].setVisible(false);
            labelForwardMutation[i].setVisible(false);
            labelBackMutation[i].setVisible(false);
            checkboxLocusUsed[i].setVisible(false);
            comboDominance[i].setVisible(false);
            editInitFraction[i].setVisible(false);
            editForwardMutation[i].setVisible(false);
            editBackMutation[i].setVisible(false);
        }
    }

    dialogLayout.addWidget(&okButton, row, 1);
    dialogLayout.addWidget(&cancelButton, row, 2);
    
    setLayout(&dialogLayout);

    setWindowTitle(tr("Loci and mutation"));

    for (int locus = 0; locus < maxLoci; locus++) {
        connect(&checkboxLocusUsed[locus], &QCheckBox::stateChanged, [this]() {
            checkGrey();  // grey out disabled fields
            });
    }

    connect(&cancelButton, &QPushButton::clicked, [this]() {
        // Cancel button. close without saving
        close();
        });

    connect(&okButton, &QPushButton::clicked, [this, parent]() {
        // OK button. save parameters
        int i;
        for (i = 0; i < maxLoci; i++) {
            if (i < parent->d.nLoci) {
                parent->d.locusUsed[i] = checkboxLocusUsed[i].isChecked();
                parent->d.dominance[i] = comboDominance[i].currentIndex();
                if (editInitFraction[i].isModified()) {
                    parent->d.fg0[i] = parent->readFloatField(editInitFraction[i]);
                }
                if (editForwardMutation[i].isModified()) {
                    parent->d.murate[i][0] = parent->readFloatField(editForwardMutation[i]);
                }
                if (editBackMutation[i].isModified()) {
                    parent->d.murate[i][1] = parent->readFloatField(editBackMutation[i]);
                }
            }
            else {
                parent->d.locusUsed[i] = false;
            }
        }
        close();
        parent->d.parametersChanged = true;
    });

    checkGrey();
}

void LociDialogBox::checkGrey() {  // set disabled fields grey
    for (int locus = 0; locus < maxLoci; locus++) {
        bool active = checkboxLocusUsed[locus].isChecked();
        comboDominance[locus].setEnabled(active);
        editInitFraction[locus].setEnabled(active);
        editForwardMutation[locus].setEnabled(active);
        editBackMutation[locus].setEnabled(active);        
    }
}


FitnessDialogBox::FitnessDialogBox(Altruist *parent) :
    labelFitnessFor("Fitness for:"),
    /* These label texts depend on model. They are set below:
    labelEgoNonConform(parent->d.locusUsed[locusConformity] ? "among egoists, nonconform" : "among egoists"),
    labelAltruNonConform(parent->d.locusUsed[locusConformity] ? "among altruists, nonconform" : "among altruists"),
    labelEgoConform(parent->d.locusUsed[locusConformity] ? "among egoists under conformity" : "among egoists"),
    labelAltruConform(parent->d.locusUsed[locusConformity] ? "among altruists under conformity" : "among altruists"),
    labelForEgo("Egoists: "),
    labelForAltru("Altruists: "),*/
    labelGrowthRate("Growth rate"),
    labelLeaderAdvantage("Leader advantage"),
    labelLeaderSelection("Leader selection"),
    labelSelModel("Selection model"),
    labelSurviv("Survival rate"),
    labelwarIntensity("War intensity"),
    labelGroupProperties("Group properties:"),
    labelExtinctionPattern("Extinction pattern"),
    labelWarPattern("War pattern"),
    labelExtinctionRateFor("Extinction rate for:"),
    labelExtinctionRate0("small group with egoists:"),
    labelExtinctionRate1("big group with egoists:"),
    labelExtinctionRate2("small group with altruists:"),
    labelExtinctionRate3("big group with altruists:"),
    labelHaystackPeriod("haystack period"),
    labelMixingPeriod("mixing period"),
    labelDependGroupSuccessOn("Dependence of group success on:"),
    labelGroupFitnessFunction("Group fitness function:"),
    labelGroupFitnessCurvature("curvature exponent:"),

    okButton("OK", this),
    cancelButton("Cancel", this)
{
    int row = 0;         // row in dialogLayout grid
    int column = 0;      // column in dialogLayout grid
    int i, j;            // temporary indexes
    AltruData * d = &parent->d;
    bool conformityUsed = (d->bIndividualPropertiesUsed & 0x10) && d->locusUsed[locusConformity]; // conformity used

    setWindowTitle("Individual and group fitness");

    // set label texts that depend on model, fit[]
    for (i = 0; i < 2; i++) {
        labelFitRow[i].setText(QString(d->sPhenotypeName[d->fitRowLocus][i]) + "s:");
    }
    QString coltext1 = "";
    switch(d->fitConditionName) {
    case 1: 
        coltext1 = "among ";  break;
    case 2: 
        coltext1 = "under ";  break;
    }

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            QString coltext2 = d->sPhenotypeName[d->fitColLocus][i];
            if (d->fitConditionName == 1) coltext2 = coltext2 + "s";
            QString coltext3 = "";
            if (conformityUsed) coltext3 = QString(", ") + d->sPhenotypeName[d->fitSupColLocus][j];
            labelFitColumn[i+2*j].setText(coltext1 + coltext2 + coltext3);
        }
    }
    dialogLayout.addWidget(&labelFitnessFor, row, 0);
    dialogLayout.addWidget(&labelFitColumn[0], row, 1);
    dialogLayout.addWidget(&labelFitColumn[1], row, 2);
    if (conformityUsed) {
        dialogLayout.addWidget(&labelFitColumn[2], row, 3);
        dialogLayout.addWidget(&labelFitColumn[3], row, 4);
    }
    else {
        labelFitColumn[2].setVisible(false);
        labelFitColumn[3].setVisible(false);
    }
    dialogLayout.addWidget(&labelFitRow[0], ++row, 0);
    dialogLayout.addWidget(&labelFitRow[1], row+1, 0);

    for (i = 0; i < 4; i++) {
        if (d->bIndividualPropertiesUsed & 1 << i) {
            parent->writeFloatField(editIndividualFitness[i], d->fit[i]);
            dialogLayout.addWidget(&editIndividualFitness[i], row + (i >> 1), (i & 1) + 1);
        }
        else {
            editIndividualFitness[i].setVisible(false);
        }
    }
    for (; i < 8; i++) {
        if (conformityUsed) {
            parent->writeFloatField(editIndividualFitness[i], d->fit[i]);
            dialogLayout.addWidget(&editIndividualFitness[i], row + ((i & 2) >> 1), (i & 1) + 3);
        }
        else {
            editIndividualFitness[i].setVisible(false);
        }
    }
    row += 2;

    // combo box for selection model
    if (d->bSelModels) {
        int reducedIndex = compressIndex(d->selectionModel, d->bSelModels);
        for (i = 0, j = 0; i < 5; i++) {
            if (d->bSelModels & 1 << i) {  // pattern (1 << i) supported
                comboSelModel.addItem(d->selModelName[i], j++);
            }
        }
        comboSelModel.setCurrentIndex(reducedIndex);
        labelSelModel.setBuddy(&comboSelModel);
        dialogLayout.addWidget(&labelSelModel, row, 0);
        dialogLayout.addWidget(&comboSelModel, row, 1);
        if (d->bGroupPropertiesUsed & 0x40) {
            dialogLayout.addWidget(&labelGrowthRate, row, 2);
            dialogLayout.addWidget(&editGrowthRate, row, 3);
            labelGrowthRate.setBuddy(&editGrowthRate);
            parent->writeFloatField(editGrowthRate, d->growthRate);
            editGrowthRate.setEnabled(d->selectionModel != selectionFecundity); // no growth rate for fecundity selection
        }
    }
    else {
        labelSelModel.setVisible(false);
        comboSelModel.setVisible(false);
        if (d->bGroupPropertiesUsed & 0x40) {
            dialogLayout.addWidget(&labelGrowthRate, row, 0);
            dialogLayout.addWidget(&editGrowthRate, row, 1);
            labelGrowthRate.setBuddy(&editGrowthRate);
            editGrowthRate.setEnabled(true);
        }
    }
    row++; column = 0;

    // set label texts that depend on model, fit2[]
    for (int i = 0; i < 2; i++) {
        if (d->bIndividualPropertiesUsed & 0x300 << (i*2)) {
            labelFit2[i].setText(QString("Relative fitness for ") + d->sPhenotypeName[i+1][1] + "s:");
            dialogLayout.addWidget(&labelFit2[i], row, 0);
            for (int j = 0; j < 2; j++) {
                if (d->bIndividualPropertiesUsed & 0x100 << (i*2+j)) {
                    parent->writeFloatField(editIndividualFitness2[i*2+j], d->fit2[i*2+j]);
                    dialogLayout.addWidget(&editIndividualFitness2[i*2+j], row, j+1);
                }
                else {
                    editIndividualFitness2[i*2+j].setVisible(false);
                }
            }
            row++;
        }
        else {
            labelFit2[i].setVisible(false);
        }
    }

    if (d->bGroupPropertiesUsed & 0x10000) { // leader advantage
        dialogLayout.addWidget(&labelLeaderAdvantage, row, column);
        dialogLayout.addWidget(&editLeaderAdvantage, row, column+1);
        labelLeaderAdvantage.setBuddy(&editLeaderAdvantage);
        parent->writeFloatField(editLeaderAdvantage, d->leaderAdvantage);
        column += 2;
    }
    else {
        labelLeaderAdvantage.setVisible(false);
        editLeaderAdvantage.setVisible(false);
    }
    if (d->bGroupPropertiesUsed & 0x20000) { // leader leaderSelection
        dialogLayout.addWidget(&labelLeaderSelection, row, column);
        dialogLayout.addWidget(&editLeaderSelection, row, column+1);
        labelLeaderSelection.setBuddy(&editLeaderSelection);
        parent->writeFloatField(editLeaderSelection, d->leaderSelection);
        column += 2;
    }
    else {
        labelLeaderSelection.setVisible(false);
        editLeaderSelection.setVisible(false);
    }
    if (column > 0) {
        row++; column = 0;
    }

    // Group properties:
    dialogLayout.setRowMinimumHeight(row++, 16);
    dialogLayout.addWidget(&labelGroupProperties, row++, 0);

    if (d->bFitFunc) {
        dialogLayout.addWidget(&labelGroupFitnessFunction, row, 0);
        dialogLayout.addWidget(&comboGroupFitnessFunction, row, 1);
        labelGroupFitnessFunction.setBuddy(&comboGroupFitnessFunction);
        dialogLayout.addWidget(&labelGroupFitnessCurvature, row, 2);
        dialogLayout.addWidget(&editGroupFitnessCurvature, row++, 3);
        labelGroupFitnessCurvature.setBuddy(&editGroupFitnessCurvature);
        static const char * groupFitNames[5] = {
            "one is enough", "convex", "linear", "concave", "all or nothing"};
        for (i = 0, j = 0; i < 5; i++) {
            groupfitFunctions[i] = 0;
            if (d->bFitFunc & 1 << i) {  // pattern (1 << i) supported
                comboGroupFitnessFunction.addItem(groupFitNames[i], j);
                groupfitFunctionsRev[i] = j;
                groupfitFunctions[j++] = i;
            }
            else {
                groupfitFunctionsRev[i] = 1;
            }
        }
        if (d->groupFitCurvature <= 0.0f) groupfitIndex = 0;
        else if (d->groupFitCurvature < 1.0f) groupfitIndex = 1;
        else if (d->groupFitCurvature == 1.0f) groupfitIndex = 2;
        else if (d->groupFitCurvature < 1.0E30) groupfitIndex = 3;
        else groupfitIndex = 4;
        comboGroupFitnessFunction.setCurrentIndex(groupfitFunctionsRev[groupfitIndex]);
        parent->writeFloatField(editGroupFitnessCurvature, d->groupFitCurvature);
    }
    else {
        labelGroupFitnessFunction.setVisible(false);
        comboGroupFitnessFunction.setVisible(false);
        labelGroupFitnessCurvature.setVisible(false);
        editGroupFitnessCurvature.setVisible(false);
    }
    groupfitModified = false;

    // combo box for extinction pattern
    if (d->bExtinctionPatterns) {
        int reducedIndex = compressIndex(d->extinctionPattern, d->bExtinctionPatterns);
        for (i = 0, j = 0; i < ARR_LEN(d->extinctionPatternName); i++) {
            if (d->bExtinctionPatterns & 1 << i) {  // pattern (1 << i) supported
                comboExtinctionPattern.addItem(d->extinctionPatternName[i], j++);
            }
        }
        comboExtinctionPattern.setCurrentIndex(reducedIndex);
        labelExtinctionPattern.setBuddy(&comboExtinctionPattern);
        dialogLayout.addWidget(&labelExtinctionPattern, row, 0);
        dialogLayout.addWidget(&comboExtinctionPattern, row, 1);
        column = 2;
    }
    else {
        labelExtinctionPattern.setVisible(false);
        comboExtinctionPattern.setVisible(false);
        column = 0;
    }

    // combo box for war pattern
    if (d->bWarPatterns) {
        int reducedIndex = compressIndex(d->warPattern, d->bWarPatterns);
        for (i = 0, j = 0; i < ARR_LEN(d->warPatternName); i++) {
            if (d->bWarPatterns & 1 << i) {  // pattern (1 << i) supported
                comboWarPattern.addItem(d->warPatternName[i], j++);
            }
        }
        comboWarPattern.setCurrentIndex(reducedIndex);
        labelWarPattern.setBuddy(&comboWarPattern);
        dialogLayout.addWidget(&labelWarPattern, row, column);
        dialogLayout.addWidget(&comboWarPattern, row, column+1);
        column += 2;
    }
    else {
        labelWarPattern.setVisible(false);
        comboWarPattern.setVisible(false);
    }
    if (column > 2) {
        row++; column = 0;
    }
    
    if (d->bGroupPropertiesUsed & 0x10) { // survival degree used
        dialogLayout.addWidget(&labelSurviv, row, column);
        dialogLayout.addWidget(&editSurviv, row, column+1);
        labelSurviv.setBuddy(&editSurviv);
        parent->writeFloatField(editSurviv, d->surviv);
        column += 2;
    }
    else {
        labelSurviv.setVisible(false);
        editSurviv.setVisible(false);
    }
    if (column > 2) {
        row++; column = 0;
    }

    if (d->bGroupPropertiesUsed & 0x80) { // war intensity used
        dialogLayout.addWidget(&labelwarIntensity, row, column);
        dialogLayout.addWidget(&editwarIntensity, row, column+1);
        labelwarIntensity.setBuddy(&editwarIntensity);
        parent->writeFloatField(editwarIntensity, d->warIntensity);
        column += 2;
    }
    else {
        labelwarIntensity.setVisible(false);
        editwarIntensity.setVisible(false);
    }
    if (column) {
        row++; column = 0;
    }

    if (d->bGroupPropertiesUsed & 7) { // group extinction rates
        dialogLayout.addWidget(&labelExtinctionRateFor, row++, 0);
        dialogLayout.addWidget(&labelExtinctionRate0, row, 0);
        dialogLayout.addWidget(&editExtinctionRates[0], row, 1);
        dialogLayout.addWidget(&labelExtinctionRate1, row, 2);
        dialogLayout.addWidget(&editExtinctionRates[1], row++, 3);
        dialogLayout.addWidget(&labelExtinctionRate2, row, 0);
        dialogLayout.addWidget(&editExtinctionRates[2], row, 1);
        dialogLayout.addWidget(&labelExtinctionRate3, row, 2);
        dialogLayout.addWidget(&editExtinctionRates[3], row++, 3);
        for (i = 0; i < 4; i++) {
            parent->writeFloatField(editExtinctionRates[i], d->extinctionRate[i]);
        }
    }
    else {
        labelExtinctionRateFor.setVisible(false);
        labelExtinctionRate0.setVisible(false);
        labelExtinctionRate1.setVisible(false);
        labelExtinctionRate2.setVisible(false);
        labelExtinctionRate3.setVisible(false);
        for (i = 0; i < 4; i++) {
            editExtinctionRates[i].setVisible(false);
        }
    }
    if (column) {
        row++; column = 0;
    }

    if (d->bGroupPropertiesUsed & 0x100) { // haystackPeriod used
        dialogLayout.addWidget(&labelHaystackPeriod, row, column++);
        dialogLayout.addWidget(&editHaystackPeriod, row, column++);
        labelHaystackPeriod.setBuddy(&editHaystackPeriod);
        parent->writeIntField(editHaystackPeriod, d->haystackPeriod);
    }
    else {
        labelHaystackPeriod.setVisible(false);
        editHaystackPeriod.setVisible(false);
    }
    if (d->bGroupPropertiesUsed & 0x200) { // mixingPeriod used
        dialogLayout.addWidget(&labelMixingPeriod, row, column++);
        dialogLayout.addWidget(&editMixingPeriod, row, column++);
        labelMixingPeriod.setBuddy(&editMixingPeriod);
        parent->writeIntField(editMixingPeriod, d->mixingPeriod);
    }
    else {
        labelMixingPeriod.setVisible(false);
        editMixingPeriod.setVisible(false);
    }
    if (column > 2) {
        row++; column = 0;
    }

    int m;     // count model-specific parameter fields
    // clear unused model-specific parameter widgets
    labelModelSpecificParameters.setVisible(false);
    for (m = 0; m < maxModelSpecificParameters; m++) {
        labelModelSpecific[m].setVisible(false);
        editModelSpecific[m].setVisible(false);
        checkboxModelSpecific[m].setVisible(false);
    }
    // add fields for model-specific parameters
    int imodel = d->iModel; 
    if (imodel >= models.getNum()) imodel = 0;
    ModelDescriptor * md = models.getModel(imodel);
    ParameterDef const * pardef = md->extraParameters; // parameter definition record
    m = 0;
    // loop through model-specific parameter list
    for (; pardef->num != 0 || pardef->name != nullptr; pardef++) {
        if (pardef->name == nullptr || pardef->type == 0) continue;
        if (m == 0) {  // first model-specific parameter
            labelModelSpecificParameters.setText("Model-specific parameters:");
            labelModelSpecificParameters.setVisible(true);
            dialogLayout.addWidget(&labelModelSpecificParameters, row++, 0);
        }
        // copy pardef->name and replace underscores with spaces
        char parname[40];
        for (int i = 0; i < sizeof(parname); i++) {
            char c = pardef->name[i];
            if (c == '_') c = ' ';
            parname[i] = c;
            if (c == 0) break;        
        }
        parname[sizeof(parname)-1] = 0;
        labelModelSpecific[m].setText(parname);
        labelModelSpecific[m].setVisible(true);
        dialogLayout.addWidget(&labelModelSpecific[m], row, column++);
        if (pardef->type == 3) {
            checkboxModelSpecific[m].setVisible(true);
            dialogLayout.addWidget(&checkboxModelSpecific[m], row, column++);
            labelModelSpecific[m].setBuddy(&checkboxModelSpecific[m]);
            checkboxModelSpecific[m].setChecked(d->modelspec_i[m]);
        }
        else {
            editModelSpecific[m].setVisible(true);
            dialogLayout.addWidget(&editModelSpecific[m], row, column++);
            labelModelSpecific[m].setBuddy(&editModelSpecific[m]);
            if (pardef->type == 2) {
                parent->writeFloatField(editModelSpecific[m], d->modelspec_f[m]);
            }
            else {
                parent->writeIntField(editModelSpecific[m], d->modelspec_i[m]);
            }
        }
        if (column > 3) {
            row++; column = 0;
        }
        if (++m > maxModelSpecificParameters) break;
    }
    if (column > 0) {
        column = 0; row++;
    }
    dialogLayout.addWidget(&okButton, row, 1);
    dialogLayout.addWidget(&cancelButton, row, 2);

    setLayout(&dialogLayout);

    connect(&comboSelModel, &QComboBox::activated, [this, parent]() {
        // selection model combo box index changed
        int selectionModel = expandIndex(comboSelModel.currentIndex(), parent->d.bSelModels);
        editGrowthRate.setEnabled(selectionModel != 0);
        });

    connect(&comboGroupFitnessFunction, &QComboBox::activated, [this, parent]() {
        // group fitness combo box index changed. change editGroupFitnessCurvature accordingly
        groupfitModified = true;
        int j = groupfitFunctions[comboGroupFitnessFunction.currentIndex()];
        float gf;
        union {
            uint32_t i = 0x7F800000;
            float f;
        } infinity;
        if (j != groupfitIndex) {
            switch (j) {
            case 0: gf = 0.0f;  break;
            case 1: gf = 0.5f;  break;
            case 2: gf = 1.0f;  break;
            case 3: gf = 2.0f;  break;
            default: gf = infinity.f;
            }
            parent->writeFloatField(editGroupFitnessCurvature, gf);
            groupfitIndex = j;
        }
        });

    connect(&editGroupFitnessCurvature, &QLineEdit::editingFinished, [this, parent]() {
        groupfitModified = true;
        float gf = parent->readFloatField(editGroupFitnessCurvature);

        if (gf <= 0.0f) groupfitIndex = 0;
        else if (gf < 1.0f) groupfitIndex = 1;
        else if (gf == 1.0f) groupfitIndex = 2;
        else if (gf < 1.0E30) groupfitIndex = 3;
        else groupfitIndex = 4;
        comboGroupFitnessFunction.setCurrentIndex(groupfitFunctionsRev[groupfitIndex]);
        });

    connect(&cancelButton, &QPushButton::clicked, [this]() {
        // Cancel button. close without saving
        close();
        });

    connect(&okButton, &QPushButton::clicked, [this, parent, d]() {
        // OK button. save parameters
        int i;
        for (i = 0; i < 8; i++) {
            if (editIndividualFitness[i].isModified()) {
                d->fit[i] = parent->readFloatField(editIndividualFitness[i]);
            }
        }
        for (i = 0; i < 4; i++) {
            if (editIndividualFitness2[i].isModified()) {
                d->fit2[i] = parent->readFloatField(editIndividualFitness2[i]);
            }
        }
        if (d->bSelModels) {
            d->selectionModel = expandIndex(comboSelModel.currentIndex(), d->bSelModels);
        }
        if (editGrowthRate.isModified()) {
            d->growthRate = parent->readFloatField(editGrowthRate);
        }
        if (editLeaderAdvantage.isModified()) {
            d->leaderAdvantage = parent->readFloatField(editLeaderAdvantage);
        }
        if (editLeaderSelection.isModified()) {
            d->leaderSelection = parent->readFloatField(editLeaderSelection);
        }

        if (d->bExtinctionPatterns) {        
            d->extinctionPattern = expandIndex(comboExtinctionPattern.currentIndex(), d->bExtinctionPatterns);
        }

        if (d->bWarPatterns) {        
            d->warPattern = expandIndex(comboWarPattern.currentIndex(), d->bWarPatterns);
        }

        if (editSurviv.isModified()) {
            d->surviv = parent->readFloatField(editSurviv);
        }
        if (editwarIntensity.isModified()) {
            d->warIntensity= parent->readFloatField(editwarIntensity);
        }

        if (editHaystackPeriod.isModified()) {
            d->haystackPeriod = parent->readIntField(editHaystackPeriod);
        }
        if (editMixingPeriod.isModified()) {
            d->mixingPeriod = parent->readIntField(editMixingPeriod);
        }

        if (d->bGroupPropertiesUsed & 7) { // group extinction rates
            for (i = 0; i < 4; i++) {
                d->extinctionRate[i] = parent->readFloatField(editExtinctionRates[i]);
            }
        }
        if (groupfitModified) {
            d->groupFitCurvature = parent->readFloatField(editGroupFitnessCurvature);
            d->fitfunc = groupfitIndex;
        }

        // save model-specific parameters
        int imodel = d->iModel;
        if (imodel >= models.getNum()) imodel = 0;
        ModelDescriptor * md = models.getModel(imodel);
        ParameterDef const * pardef = md->extraParameters; // parameter definition record
        int m = 0;
        // loop through model-specific parameter list
        for (; pardef->num != 0 || pardef->name != nullptr; pardef++) {
            if (pardef->name == nullptr || pardef->type == 0) continue;
            switch (pardef->type) {
            case 1:  // int
                if (editModelSpecific[m].isModified()) {
                    d->modelspec_i[m] = parent->readIntField(editModelSpecific[m]);
                }
                break;
            case 2:  // float
                if (editModelSpecific[m].isModified()) {
                    d->modelspec_f[m] = parent->readFloatField(editModelSpecific[m]);
                }
                break;
            case 3:  // bool
                d->modelspec_i[m] = checkboxModelSpecific[m].isChecked();
                break;
            }
            m++;
        }

        close();
    });
}


RunControlDialogBox::RunControlDialogBox(Altruist *parent) :
    labelRandomSeed("Random generator seed:"),
    labelMinGenerations("Minimum generations:"),
    labelMaxGenerations("Maximum generations:"),
    labelStopCriterion("Stop criterion:"),
    labelDelay("Delay (ms):"),
    labelCriterionDegree("Criterion degree:"),
    labelParameterSweeps("Parameter sweeps:"),
    okButton("OK", this),
    cancelButton("Cancel", this)
{
    d = &parent->d;
    enableSweepParameter = -1;
    int row = 0;
    setWindowTitle("Run control");

    dialogLayout.addWidget(&labelRandomSeed, row, 0);
    dialogLayout.addWidget(&editRandomSeed, row++, 1);
    labelRandomSeed.setBuddy(&editRandomSeed);
    parent->writeIntField(editRandomSeed, parent->d.seed);
    dialogLayout.addWidget(&labelMinGenerations, row, 0);
    dialogLayout.addWidget(&editMinGenerations, row, 1);
    labelMinGenerations.setBuddy(&editMinGenerations);
    parent->writeIntField(editMinGenerations, parent->d.minimumGenerations);
    dialogLayout.addWidget(&labelMaxGenerations, row, 2);
    dialogLayout.addWidget(&editMaxGenerations, row++, 3);
    labelMaxGenerations.setBuddy(&editMaxGenerations);
    parent->writeIntField(editMaxGenerations, parent->d.maximumGenerations);
    dialogLayout.addWidget(&labelStopCriterion, row, 0);
    dialogLayout.addWidget(&comboStopCriterion, row, 1);
    labelStopCriterion.setBuddy(&comboStopCriterion);
    int i = 0;
    if (d->bStopCriterionUsed & 1) comboStopCriterion.addItem("none", i++);
    if (d->bStopCriterionUsed & 2) comboStopCriterion.addItem("genotypes uniform", i++);
    if (d->bStopCriterionUsed & 4) comboStopCriterion.addItem("phenotypes uniform", i++);

    for (int locus = 0; locus < d->nLoci; locus++) {
        if (d->bStopCriterionUsed & 1 << (locus*2 + 3)) {
            comboStopCriterion.addItem(QString(d->sLocusName[locus]) + " uniform", i++);
        }
        if (d->bStopCriterionUsed & 1 << (locus*2 + 4)) {
            comboStopCriterion.addItem(QString(d->sGeneName[locus][1]) + " fixated", i++);
        }
    }
    comboStopCriterion.setCurrentIndex(compressIndex(d->stopCriterion, d->bStopCriterionUsed));
    dialogLayout.addWidget(&labelCriterionDegree, row, 2);
    dialogLayout.addWidget(&editCriterionDegree, row++, 3);
    labelCriterionDegree.setBuddy(&editCriterionDegree);
    parent->writeFloatField(editCriterionDegree, d->stopCriterionDegree);
    dialogLayout.addWidget(&labelDelay, row, 0);
    dialogLayout.addWidget(&editDelay, row++, 1);
    labelDelay.setBuddy(&editDelay);
    parent->writeIntField(editDelay, d->delayms);

    dialogLayout.setRowMinimumHeight(row++, 20);
    dialogLayout.addWidget(&labelParameterSweeps, row++, 0);

    // check which sweep parameters to enable
    SweepParameter * sw = sweepParameterList;
    i = 0;
    enableSweepParameter = 0;
    do {
        if (sw->bbGeographyParametersUsed & d->bGeographyParametersUsed || sw->bbGroupPropertiesUsed & d->bGroupPropertiesUsed) {
            enableSweepParameter |= 1 << i;            // enable sweep for this parameter
        }
        sw++;  i++;
    } while (sw->offset != 0);


    for (int isweep = 0; isweep < maxSweeps; isweep++) {
        labelSweepType[isweep].setText("Sweep type:");
        labelParameterToSweep[isweep].setText("Parameter to sweep:");
        labelStartValue[isweep].setText("Start value:");
        labelEndValue[isweep].setText("End value:");
        labelStep[isweep].setText("Step:");

        dialogLayout.addWidget(&labelSweepType[isweep], row, 0);
        dialogLayout.addWidget(&comboSweepType[isweep], row, 1);
        comboSweepType[isweep].addItem("none", 0);
        comboSweepType[isweep].addItem("linear", 1);
        comboSweepType[isweep].addItem("logarithmic", 2);
        if (isweep == 0) {
            comboSweepType[isweep].addItem("linear search", 3);
            comboSweepType[isweep].addItem("logarithmic search", 4);
        }
        comboSweepType[isweep].setCurrentIndex(d->sweepType[isweep]);
        dialogLayout.addWidget(&labelParameterToSweep[isweep], row, 2);
        dialogLayout.addWidget(&comboParameterToSweep[isweep], row++, 3);
        // loop through sweepParameterList to list possible sweep parameters
        sw = sweepParameterList;
        i = 0; int j = 0;
        do {
            if (enableSweepParameter & 1 << i) {
                comboParameterToSweep[isweep].addItem(sw->name, j++);
            }
            sw++;  i++;
        } while (sw->offset != 0);
        comboParameterToSweep[isweep].setCurrentIndex(compressIndex(d->sweepParameter[isweep], enableSweepParameter));
        dialogLayout.addWidget(&labelStartValue[isweep], row, 0);
        dialogLayout.addWidget(&editStartValue[isweep], row, 1);
        parent->writeFloatField(editStartValue[isweep], d->sweepStartValue[isweep]);
        dialogLayout.addWidget(&labelEndValue[isweep], row, 2);
        dialogLayout.addWidget(&editEndValue[isweep], row, 3);
        parent->writeFloatField(editEndValue[isweep], d->sweepEndValue[isweep]);
        dialogLayout.addWidget(&labelStep[isweep], row, 4);
        dialogLayout.addWidget(&editStep[isweep], row++, 5);
        parent->writeFloatField(editStep[isweep], d->sweepStep[isweep]);
        dialogLayout.setRowMinimumHeight(row++, 20);
    }

    dialogLayout.addWidget(&okButton, row, 1);
    dialogLayout.addWidget(&cancelButton, row, 2);

    setLayout(&dialogLayout);

    for (int isweep = 0; isweep < maxSweeps; isweep++) {
        connect(&comboSweepType[isweep], &QComboBox::currentIndexChanged, [this, isweep]() {
            checkGrey();  // grey out disabled fields
            });
    }

    connect(&comboStopCriterion, &QComboBox::currentIndexChanged, [this]() {
        checkGrey();  // grey out disabled fields
        });

    connect(&cancelButton, &QPushButton::clicked, [this]() {
        // Cancel button. close without saving
        close();
        });

    connect(&okButton, &QPushButton::clicked, [this, parent]() {
        // OK button. save parameters
        if (editRandomSeed.isModified()) {
            parent->d.seed = parent->readIntField(editRandomSeed);
        }
        if (editMinGenerations.isModified()) {
            parent->d.minimumGenerations = parent->readIntField(editMinGenerations);
        }
        if (editMaxGenerations.isModified()) {
            parent->d.maximumGenerations = parent->readIntField(editMaxGenerations);
        }
        parent->d.stopCriterion = expandIndex(comboStopCriterion.currentIndex(), parent->d.bStopCriterionUsed);
        if (editCriterionDegree.isModified()) {
            parent->d.stopCriterionDegree = parent->readFloatField(editCriterionDegree);
        }
        if (editDelay.isModified()) {
            parent->d.delayms = parent->readIntField(editDelay);
        }

        parent->d.sweepsUsed = 0;
        for (int isweep = 0; isweep < maxSweeps; isweep++) {
            int s = comboSweepType[isweep].currentIndex();
            parent->d.sweepType[isweep] = s;
            parent->d.sweepsUsed |= int(s != loopUnused) << isweep;
            parent->d.sweepParameter[isweep] = expandIndex(comboParameterToSweep[isweep].currentIndex(), enableSweepParameter);

            if (editStartValue[isweep].isModified()) {
                parent->d.sweepStartValue[isweep] = parent->readFloatField(editStartValue[isweep]);
            }
            if (editEndValue[isweep].isModified()) {
                parent->d.sweepEndValue[isweep] = parent->readFloatField(editEndValue[isweep]);
            }
            if (editStep[isweep].isModified()) {
                parent->d.sweepStep[isweep] = parent->readFloatField(editStep[isweep]);
            }
            // check for duplicate parameter
            for (int j = 0; j < isweep; j++) {
                if ((parent->d.sweepsUsed & 1 << isweep) && parent->d.sweepParameter[isweep] == parent->d.sweepParameter[j]) {
                    errorMessage("Multiple sweeps with same parameter");
                    return;
                }
            }
            // check for logarithm of non-positive
            if (s == 2 || s == 4) { // logarithmic, must be positive
                if (parent->d.sweepStartValue[isweep] <= 0.f || parent->d.sweepEndValue[isweep] <= 0.f) {
                    errorMessage("parameter values must be positive in logarithmic sweep");
                    return;
                }
            }
        }

        close();
        d->parametersChanged = true;
        });
    checkGrey();  // grey out disabled fields
}

void RunControlDialogBox::checkGrey() {  // set disabled fields grey
    for (int isweep = 0; isweep < maxSweeps; isweep++) {
        bool active = comboSweepType[isweep].currentIndex() != 0;
        editStartValue[isweep].setEnabled(active);
        editEndValue[isweep].setEnabled(active);
        editStep[isweep].setEnabled(active);
        comboParameterToSweep[isweep].setEnabled(active);
    }

    // integer variables can be used for sweep but not for search.
    // grey out integer variables for search
    bool search0 = comboSweepType[0].currentIndex() >= 3;            // sweep 0 is a search
    for (int i = 0; i < comboParameterToSweep[0].count(); i++) {     // loop through items in combo box
        int sweepParameter = expandIndex(i, enableSweepParameter);   // translate reduced index to full index
        bool isFloat = sweepParameterList[sweepParameter].varType >= varFloat; // this variable is float
        bool enableItem = isFloat || !search0;
        if (comboParameterToSweep[0].currentIndex() == i && !enableItem) {
            // a disabled item is selected. change index
            comboParameterToSweep[0].setCurrentIndex(1);
        }
        auto * model = qobject_cast<QStandardItemModel*>(comboParameterToSweep[0].model());
        if (model == 0) break;
        auto * item = model->item(i);
        if (item == 0) break;
        item->setEnabled(enableItem);            // grey out or enable item in combo box
    }

    bool criterionActive = expandIndex(comboStopCriterion.currentIndex(), d->bStopCriterionUsed) != stopCritNone;
    editCriterionDegree.setEnabled(criterionActive);
}



DataOutputDialogBox::DataOutputDialogBox(Altruist *parent) :
    labelMakeOutputFile("Make data output file:"),
    labelTitle("Title:"),
    labelOutputInterval("Output interval:"),
    labelSteadyStateAfter("Steady state after:"),
    labelInitialParameters("Initial parameters:"),
    labelGenerations("Generations:"),
    labelInhabitedIslands("Inhabited islands"),
    labelGeneFractions("Gene fractions:"),
    labelAltruists("Phenotypic altruists:"),
    labelModelSpecific("Model-specific data:"),
    labelMutations("Mutations:"),
    labelMigrants("Migrants:"),
    labelAltruismGroups("Altruism groups:"),
    labelExtinctions("Extinctions:"),
    labelResult("Simulation result:"),
    labelSteadyState("Steady state:"),
    labelTimeConsumption("Time consumption:"),
    labelSearchResult("Search results:"),
    changeButton("Change file name", this),
    okButton("OK", this),
    cancelButton("Cancel", this)
{
    int row = 0;
    setWindowTitle("Data file output");

    dialogLayout.addWidget(&labelMakeOutputFile, row, 0);
    dialogLayout.addWidget(&checkboxMakeOutputFile, row, 1);
    labelMakeOutputFile.setBuddy(&checkboxMakeOutputFile);
    checkboxMakeOutputFile.setChecked(parent->d.makeOutputFile);

    dialogLayout.addWidget(&changeButton, row++, 2);
    dialogLayout.addWidget(&labelTitle, row, 0);
    dialogLayout.addWidget(&editTitle, row++, 1);
    labelTitle.setBuddy(&editTitle);
    if (!parent->d.outputTitle.isEmpty()) editTitle.setText(parent->d.outputTitle);
    dialogLayout.addWidget(&labelOutputInterval, row, 0);
    dialogLayout.addWidget(&editOutputInterval, row, 1);
    labelOutputInterval.setBuddy(&editOutputInterval);
    parent->writeIntField(editOutputInterval, parent->d.fileOutInterval);
    dialogLayout.addWidget(&labelSteadyStateAfter, row, 2);
    dialogLayout.addWidget(&editSteadyStateAfter, row++, 3);
    labelSteadyStateAfter.setBuddy(&editSteadyStateAfter);
    parent->writeIntField(editSteadyStateAfter, parent->d.fileOutSteadyState);

    dialogLayout.addWidget(&labelInitialParameters, row, 0);
    dialogLayout.addWidget(&checkboxInitialParameters, row, 1);
    labelInitialParameters.setBuddy(&checkboxInitialParameters);
    checkboxInitialParameters.setChecked(parent->d.bOutOptions & 1);

    dialogLayout.addWidget(&labelGenerations, row, 2);
    dialogLayout.addWidget(&checkboxGenerations, row++, 3);
    labelGenerations.setBuddy(&checkboxGenerations);
    checkboxGenerations.setChecked(parent->d.bOutOptions & 2);

    dialogLayout.addWidget(&labelInhabitedIslands, row, 0);
    dialogLayout.addWidget(&checkboxInhabitedIslands, row, 1);
    labelInhabitedIslands.setBuddy(&checkboxInhabitedIslands);
    checkboxInhabitedIslands.setChecked(parent->d.bOutOptions & 4);

    dialogLayout.addWidget(&labelGeneFractions, row, 2);
    dialogLayout.addWidget(&checkboxGeneFractions, row++, 3);
    labelGeneFractions.setBuddy(&checkboxGeneFractions);
    checkboxGeneFractions.setChecked(parent->d.bOutOptions & 8);

    dialogLayout.addWidget(&labelAltruists, row, 0);
    dialogLayout.addWidget(&checkboxAltruists, row, 1);
    labelAltruists.setBuddy(&checkboxAltruists);
    checkboxAltruists.setChecked(parent->d.bOutOptions & 0x10);

    dialogLayout.addWidget(&labelMutations, row, 2);
    dialogLayout.addWidget(&checkboxMutations, row++, 3);
    labelMutations.setBuddy(&checkboxMutations);
    checkboxMutations.setChecked(parent->d.bOutOptions & 0x20);

    dialogLayout.addWidget(&labelMigrants, row, 0);
    dialogLayout.addWidget(&checkboxMigrants, row, 1);
    labelMigrants.setBuddy(&checkboxMigrants);
    checkboxMigrants.setChecked(parent->d.bOutOptions & 0x40);

    dialogLayout.addWidget(&labelAltruismGroups, row, 2);
    dialogLayout.addWidget(&checkboxAltruismGroups, row++, 3);
    labelAltruismGroups.setBuddy(&checkboxAltruismGroups);
    checkboxAltruismGroups.setChecked(parent->d.bOutOptions & 0x100);

    dialogLayout.addWidget(&labelExtinctions, row, 0);
    dialogLayout.addWidget(&checkboxExtinctions, row, 1);
    labelExtinctions.setBuddy(&checkboxExtinctions);
    checkboxExtinctions.setChecked(parent->d.bOutOptions & 0x200);

    dialogLayout.addWidget(&labelResult, row, 2);
    dialogLayout.addWidget(&checkboxResult, row++, 3);
    labelResult.setBuddy(&checkboxResult);
    checkboxResult.setChecked(parent->d.bOutOptions & 0x400);

    dialogLayout.addWidget(&labelModelSpecific, row, 0);
    dialogLayout.addWidget(&checkboxModelSpecific, row, 1);
    labelModelSpecific.setBuddy(&checkboxModelSpecific);
    checkboxModelSpecific.setChecked(parent->d.bOutOptions & 0x800);

    dialogLayout.addWidget(&labelSteadyState, row, 2);
    dialogLayout.addWidget(&checkboxSteadyState, row++, 3);
    labelSteadyState.setBuddy(&checkboxSteadyState);
    checkboxSteadyState.setChecked(parent->d.bOutOptions & 0x1000);

    dialogLayout.addWidget(&labelTimeConsumption, row, 0);
    dialogLayout.addWidget(&checkboxTimeConsumption, row, 1);
    labelTimeConsumption.setBuddy(&checkboxTimeConsumption);
    checkboxTimeConsumption.setChecked(parent->d.bOutOptions & 0x2000);

    dialogLayout.addWidget(&labelSearchResult, row, 2);
    dialogLayout.addWidget(&checkboxSearchResult, row++, 3);
    labelSearchResult.setBuddy(&checkboxSearchResult);
    checkboxSearchResult.setChecked(parent->d.bOutOptions & 0x10000);

    dialogLayout.setRowMinimumHeight(row++, 20);

    dialogLayout.addWidget(&okButton, row, 1);
    dialogLayout.addWidget(&cancelButton, row, 2);
    setLayout(&dialogLayout);

    connect(&changeButton, &QPushButton::clicked, [this, parent]() {
        // change file name
        if (parent->d.outputFilePath.isEmpty()) {
            parent->d.outputFilePath = parent->parameterFilePath;
        }            
        QString filename = QFileDialog::getSaveFileName(this,
            "Data output file", parent->d.outputFilePath, "*.csv;*.txt");
        QFileInfo fi(filename);
        parent->d.outputFileName = fi.fileName();
        parent->d.outputFilePath = fi.absolutePath();
        parent->worker->resultSets = 0;       // count result sets in output file
        });

    connect(&cancelButton, &QPushButton::clicked, [this]() {
        // Cancel button. close without saving
        close();
        });

    connect(&okButton, &QPushButton::clicked, [this, parent]() {
        // OK button. save parameters
        Worker * worker = parent->worker;
        if (editTitle.isModified()) {
            parent->d.outputTitle = editTitle.text();
        }
        if (editOutputInterval.isModified()) {
            parent->d.fileOutInterval = parent->readIntField(editOutputInterval);
        }
        if (editSteadyStateAfter.isModified()) {
            parent->d.fileOutSteadyState = parent->readIntField(editSteadyStateAfter);
        }
        parent->d.bOutOptions &= ~0xFFFFF;
        parent->d.bOutOptions |= (int)checkboxInitialParameters.isChecked();
        parent->d.bOutOptions |= (int)checkboxGenerations.isChecked() << 1;
        parent->d.bOutOptions |= (int)checkboxInhabitedIslands.isChecked() << 2;
        parent->d.bOutOptions |= (int)checkboxGeneFractions.isChecked() << 3;
        parent->d.bOutOptions |= (int)checkboxAltruists.isChecked() << 4;
        parent->d.bOutOptions |= (int)checkboxMutations.isChecked() << 5;
        parent->d.bOutOptions |= (int)checkboxMigrants.isChecked() << 6;
        parent->d.bOutOptions |= (int)checkboxAltruismGroups.isChecked() << 8;
        parent->d.bOutOptions |= (int)checkboxExtinctions.isChecked() << 9;
        parent->d.bOutOptions |= (int)checkboxResult.isChecked() << 10;
        parent->d.bOutOptions |= (int)checkboxModelSpecific.isChecked() << 11;
        parent->d.bOutOptions |= (int)checkboxSteadyState.isChecked() << 12;
        parent->d.bOutOptions |= (int)checkboxTimeConsumption.isChecked() << 13;
        parent->d.bOutOptions |= (int)checkboxSearchResult.isChecked() << 16;

        if (parent->d.outputFileName.isEmpty() && checkboxMakeOutputFile.isChecked()) {
            errors.reportError("Filename missing");
            parent->d.makeOutputFile = false;
        }
        else {
            parent->d.makeOutputFile = checkboxMakeOutputFile.isChecked();
        }
        close();
        });
}
