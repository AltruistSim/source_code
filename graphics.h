/*******************************  graphics.h   ********************************
* Author:        Agner Fog
* Date created:  2023-08-07
* Last modified: 2023-12-31
* Version:       3.001
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines graphic display of islands, territories, etc.
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#pragma once

// values for graphicsType
const int graphicsNone       = 0;                // no graphics display
const int graphicsIslands    = 1;                // show 2-dimensionsl map of islands
const int graphicsTerritories= 2;                // show 2-dimensionsl map of group territories with floating boundaries
const int graphicsLimits     = 3;                // show x-y map of limiting parameter ranges for altruism/egoism/polymorphism

// other constants
const int maxGeographicRows = 32;                // maximum number of rows to show in geographic map
const int maxGeographicColumns = 32;             // maximum number of columns to show in geographic map
const int maxIslandsShown = maxGeographicRows * maxGeographicColumns; // maximum number of islands to show in geographic map
const int gridUnit = 32;                         // distance between islands on graphic display
const int maxPolygons = 1024;                    // maximum number of polygons to draw
const int maxTextItems = 64;                     // maximum number of texts to draq
const int xAxisLength = 256;                     // length of x axis in parameter map
const int yAxisLength = 256;                     // length of y axis in parameter map

class AltruistView : public QGraphicsView {
Q_OBJECT

public:
    AltruistView(Altruist * a);                  // constructor
    void draw();                                 // draw graphics
    void drawIslands();                          // draw map of islands
    void drawTerritories();                      // draw map of territories
    void drawLimitMap();                         // draw x-y map of limiting parameter values
    void clear();                                // clear drawing

public slots:
    void zoomIn();                               // zoom in graphics view
    void zoomOut();                              // zoom out graphics view
    void mousePressEvent(QMouseEvent * event);   // handle mouse click

protected:
    QGraphicsScene scene;                        // container for graphics items
    Altruist * main;                             // pointer to main Altruist object
    int squaresUsed;                             // number of squared in use
    int polygonsUsed;                            // number of polygons in use
    int textItemsUsed;                           // number of text items in use

    QGraphicsRectItem squares[maxIslandsShown];  // islands to draw
    QGraphicsPolygonItem polygons[maxPolygons];  // polygons to draw
    QGraphicsTextItem textItems[maxTextItems];   // text to draw
    void mouseHandlerIslands(int x, int y);      // handle mouse click on islands graphic
    void mouseHandlerTerritories(int x, int y);  // handle mouse click on territories graphic
    void mouseHandlerLimitMap(int x, int y);     // handle mouse click on limits map graphic
};

void demeDescriptionMessageBox(int id, AltruData * d); // make a message box showing the contents of a deme record
