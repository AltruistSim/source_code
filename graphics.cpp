/****************************  graphics.cpp   *********************************
* Author:        Agner Fog
* Date created:  2023-08-07
* Last modified: 2024-10-13
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines functions for drawing a geograpic map with islands
* or territories and graphic presentations of parameter ranges.
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

AltruistView::AltruistView(Altruist * a) : QGraphicsView(), scene(this) {
    // constructor
    main = a;
    main->setCentralWidget(this);
    setScene(&scene);
    for (int i = 0; i < maxIslandsShown; i++) {
        scene.addItem(&squares[i]);
    }
    for (int i = 0; i < maxPolygons; i++) {
        scene.addItem(&polygons[i]);
    }
    for (int i = 0; i < maxTextItems; i++) {
        scene.addItem(&textItems[i]);
    }

    squaresUsed = 0;  polygonsUsed = 0;  textItemsUsed = 0 ;
    // Add dummy invisible rectangle. 
    // The program hangs if there is nothing on the scene, I don't know why?
    QGraphicsRectItem * rect = new QGraphicsRectItem(0,0,10,10);
    rect->setVisible(false);
    scene.addItem(rect);
    //clear();
}

void AltruistView::zoomIn() {
    // zoom in graphics view
    scale(1.25, 1.25); 
}

void AltruistView::zoomOut() {
    // zoom out graphics view
    scale(1. / 1.25, 1. / 1.25); 
}

void AltruistView::draw() {
    if (main->worker->lastResult == resultUnknown && main->d.runState == state_start) {
        clear();
    }
    switch (main->d.graphicsType) {
    case graphicsNone:
        clear();
        break;
    case graphicsIslands:
        drawIslands();
        break;
    case graphicsTerritories:
        drawTerritories();
        break;
    case graphicsLimits:
        drawLimitMap();
        break;
    }
}

void AltruistView::clear() {
    // clear drawing
    int i;
    // hide unused squares
    for (i = 0; i < squaresUsed && i < maxIslandsShown; i++) {
        squares[i].setRect(0, 0, 0, 0);
        squares[i].setVisible(false);
    }
    // hide unused polygons
    for (i = 0; i < polygonsUsed && i < maxPolygons; i++) {
        polygons[i].setPolygon(QPolygonF());
        polygons[i].setVisible(false);
    }
    // hide unused text
    for (i = 0; i < maxTextItems; i++) {
        textItems[i].setPos(0, 0);
        textItems[i].setRotation(0);
        textItems[i].setVisible(false);
    }
    squaresUsed = 0;  polygonsUsed = 0;
}

void AltruistView::drawIslands() {
    if (polygonsUsed) clear();                   // hide unused polygons

    // draw islands or territories
    AltruData * d = &main->d;                    // point to main data
    if (d == 0 || d->currentModel == 0) return;  // not initialized yet

    // number of rows and columns to show
    int32_t numColumns = d->rowLength;
    int32_t numRows = (d->nIslands + numColumns - 1) / numColumns;
    d->numRows = numRows;
    int32_t row, column;
    int32_t index = 0;
    int8_t * group;                              // point to current group
    int maxSquare = gridUnit - 4;                // max size of squares in graphics
    if (maxSquare < 2) maxSquare = 2;

    // inspect group descriptor to find offset for size and group count
    int32_t nn = 0;                              // group population
    int32_t nref = d->nMaxPerGroup;               // max group population
    int32_t ga0 = 0;                             // gene count for mutant at primary locus. determines hue or first color
    int32_t ga1 = 0;                             // gene count for mutant at secondary locus. determines second color
    int32_t ga2 = 0;                             // gene count for mutant at third locus. determines third color
    float groupProperty = 0;                     // group fitness or other group property
    int onn = -1;                                // offset of population count in group data structure
    int oga0 = 0;                                // offset of locus 0 mutant gene in group data structure
    int oga1 = 0;                                // offset of locus 1 mutant gene in group data structure
    int oga2 = 0;                                // offset of locus 2 mutant gene in group data structure
    int oGroupProperty = 0;                      // offset of groupProperty in group data structure
    int varType = 0;
    int oddMove = 0;                             // shift of odd rows in honeycomb organization
    int g;      // graphics field type
    int i = 0;  // index into groupFields
    do {    // search through group descriptor list
        g = d->currentModel->groupFields[i].graphics;
        switch (g) {
        case 3:   // field contains population size
            onn = d->currentModel->groupFields[i].offset;
            varType = d->currentModel->groupFields[i].varType; 
            break;
        case 10:  // field contains first gene count used for selecting color
            oga0 = d->currentModel->groupFields[i].offset;
            break;
        case 11:  // field contains second gene count used for selecting color
            if (d->locusUsed[1]) oga1 = d->currentModel->groupFields[i].offset;
            break;
        case 12:  // field contains third gene count used for selecting color
            if (d->locusUsed[2]) oga2 = d->currentModel->groupFields[i].offset;
            break;
        case 20:  // field contains group property
            if (d->currentModel->groupFields[i].varType == varFloat) { // float            
                oGroupProperty = d->currentModel->groupFields[i].offset;
            }
            break;
        }
    } while (d->currentModel->groupFields[i++].type != 0);  // group descriptor list ends with 0 

    // loop through rows and columns
    for (row = 0; row < numRows; row++) {
        if (d->migrationTopology == topologyHoneycomb && (row & 1)) {
            oddMove = gridUnit >> 1;
        }
        else {
            oddMove = 0;
        }
        for (column = 0; column < numColumns; column++) {
            index = row * numColumns + column;
            if (index >= d->nIslands || index >= maxIslandsShown) break;  // stop if last row is incomplete or too many islands
            group = d->groupData + index * d->groupStructureSize; // point to current group
            // extract data for this group according to structure defined in group descriptors:
            if (varType == varInt16) {
                if (onn >= 0) nn = *(int16_t*)(group + onn);
                if (oga0) ga0 = *(int16_t*)(group + oga0);
                if (oga1) ga1 = *(int16_t*)(group + oga1);
                if (oga2) ga2 = *(int16_t*)(group + oga2);
            }
            else if (varType == varInt32) {
                if (onn >= 0) nn = *(int32_t*)(group + onn);
                if (oga0) ga0 = *(int32_t*)(group + oga0);
                if (oga1) ga1 = *(int32_t*)(group + oga1);
                if (oga2) ga2 = *(int32_t*)(group + oga2);
            }

            // use sqrt to make area proportional with population
            int32_t size = lround(sqrtf(nn / float(nref)) * float(maxSquare)); // size of square
            QColor color(QColor(0xFF, 0xFF, 0xFF));        // white if nn = 0
            if (nn > 0) {                                  // avoid division by 0
                if (oga1 || oga2) {                        // multiple loci. red, green, blue for each locus
                    int red   = (ga0 * 0xFF + nn) / (nn*2);// red   represents fraction of mutant at first locus
                    int blue  = (ga1 * 0xFF + nn) / (nn*2);// green represents fraction of mutant at second locus
                    int green = (ga2 * 0xFF + nn) / (nn*2);// blue  represents fraction of mutant at third locus
                    color = QColor(red, green, blue);
                }
                else {                                     // single locus. express gene fraction with hue
                    int hue = (nn * 2 - ga0) * 240 / (2 * nn);// 240 = blue, 0 = red
                    color = QColor::fromHsv(hue, 0xFF, 0xFF); // get color with hue from blue to red
                }
            }
            if (d->nIslands == 1) squares[index].setRect(0, 0, size, size); // one big square
            else squares[index].setRect(column * gridUnit - (size >> 1) + oddMove, row * gridUnit - (size >> 1), size, size);
            squares[index].setBrush(QBrush(color));
        }
        if (index >= maxIslandsShown) break;
    }
    int textIndex = 0;
    int islandsShown = d->nIslands;
    if (islandsShown > maxIslandsShown) {
        islandsShown = maxIslandsShown;
        // tell that only an excerpt is shown
        textItems[0].setPlainText("... only excerpt shown ...");
        textItems[0].setPos(QPointF(0, (row+1) * gridUnit));
        textItems[0].setRotation(0);
        textItems[0].setVisible(true);
        textIndex++;
    }

    if (squaresUsed != islandsShown) {
        squaresUsed = islandsShown;
        // show used squares, hide unused squares
        for (index = 0; index < islandsShown; index++) {
            squares[index].setVisible(true);
        }
        for (; index < maxIslandsShown; index++) {
            squares[index].setRect(0, 0, 0, 0);
            squares[index].setVisible(false);
        }
    }    
    if (textIndex != textItemsUsed) { // hide unused text items
        for (i = textIndex; i < maxTextItems; i++) {
            textItems[i].setPos(0, 0);
            textItems[i].setRotation(0);
            textItems[i].setVisible(false);
        }
        textItemsUsed = textIndex;
    }
}

void formatAxis(float start, float end, int maxLabels, float * unit, int * firstTick, int * lastTick, 
    int * ticsPerLabel, char formatString[]) {
    // find appropriate format for a linear axes in a coordinate system
    // inputs: 
    // start, end: first and last data point
    // maxLabels: maximum number of numeric labels that the axis has space for
    // outputs:
    // unit: value difference between tics on axis
    // firstTick: first axis tick position is unit * firstTick
    // lastTick:  last axis tick position is unit * lastTick
    // ticsPerLabel: number of ticks between numeric labels
    // formatString: format string for numeric labels
    if (end < start) { // make sure start < end
        float temp = end;  end = start;  start = temp;    
    }
    float length = end - start;                            // length of axis in value units
    int maxTicks = 25;                                     // maximum number of ticks, you may change this
    float unit1 = length / maxTicks;                       // unit before rounding
    if (unit1 <= 0.f) unit1 = 1.f;

    // round unit1 to 1, 2, or 5 * 10^exponent
    float exponent = floorf(log10f(unit1) + 0.01f);        // decimal exponent of unit
    float fraction1 = unit1 * powf(10.f, -exponent);       // unit1 = decimals * 10^exponent
    float fraction;                                        // set fraction to 1, 2, or 5
    if (fraction1 < 1.41f) fraction = 1.f;
    else if (fraction1 < 3.16f) fraction = 2.f;
    else fraction = 5.f;
    *unit = fraction * powf(10.f, exponent);               // value increment per tick, rounded
    *firstTick = (int)lround(start / *unit);               // position of first tick, in 'unit' units
    if (*firstTick * *unit < start - 0.01f * *unit) (*firstTick)++;
    *lastTick = (int)lround(end / *unit);                  // position of last tick, in 'unit' units
    if (*lastTick * *unit > end + 0.01f * *unit) (*lastTick)--;
    int numTicks = *lastTick - *firstTick;                 // number of ticsk - 1
    int extraDecimal = 0;

    // place numeric labels at round intervals
    if (fraction == 1.f) {
        if (numTicks <= maxLabels * 10) {
            *ticsPerLabel = 10;  //extraDecimal = 1;
        }
        else if (numTicks <= maxLabels * 20) {
            *ticsPerLabel = 20;  //extraDecimal = 1;
        }
        else {
            *ticsPerLabel = 50;  //extraDecimal = 1;
        }
    }
    else if (fraction == 2.f) {
        int etick = abs(*firstTick);
        if (abs(*lastTick) > etick) etick = abs(*lastTick);
        if (etick % 5 == 0 && etick % 10 != 0) {
            // avoid fraction = 2 if start or end value is a multiple of 5, not divisible by 2
            fraction = 1.f;
            *ticsPerLabel = 5;
            exponent += 1.f;
        }
        else if (numTicks <= maxLabels * 2) {
            *ticsPerLabel = 2;
        }
        else if (numTicks <= maxLabels * 4) {
            *ticsPerLabel = 4;  //extraDecimal = 1;
        }
        else {
            *ticsPerLabel = 5;
        }
    }
    else {  // decimals = 5
        if (numTicks <= maxLabels * 5) {
            *ticsPerLabel = 5;
        }
        else if (numTicks <= maxLabels * 10) {
            *ticsPerLabel = 10;
            exponent += 1.f;
        }
        else {
            *ticsPerLabel = 50;  //extraDecimal = 1;  
            exponent += 1.f;
        }
    }
    float startEnd = fabs(start);                          // numerically largest point
    if (fabs(end) > startEnd) startEnd = fabs(end);
    int exponent1 = (int)exponent;                         // 10's exponent of unit
    int exponent2 = (int)floorf(log10f(startEnd) + 0.01f); // 10's exponent of largest value
    // find number of significant digits in numeric labels
    int formatSignificantDigits = exponent2 - exponent1 - (*ticsPerLabel >= 10) + extraDecimal;
    int digits = 0;                                        // number of digits required after decimal point
    // make appropriate format string
    if (exponent2 >= 4) {
        snprintf(formatString, 16, "%%.%iE", formatSignificantDigits);        
    }
    else if (exponent2 >= 0) {
        if (exponent1 < 0) {
            int dig = -exponent1 + extraDecimal;
            snprintf(formatString, 16, "%%.%if", dig);
        }
        else {
            snprintf(formatString, 16, "%%.0f");
        }
    }
    else if (exponent2 >= -3) {
        int dig = - exponent1 + extraDecimal;
        if (dig < 0) dig = 0;
        snprintf(formatString, 16, "%%.%if", dig);
    }
    else {
        snprintf(formatString, 16, "%%.%iE", formatSignificantDigits);
    }
}

void printLogNumber (char * text, int n) {
    // print pow(10,n) to string
    // used for logarithmic coordinate system ticks in drawLimitMap()
    switch (n) {
    case -3:  strncpy(text, "0.001", 10);  break;
    case -2:  strncpy(text, "0.01", 10);   break;
    case -1:  strncpy(text, "0.1", 10);    break;
    case  0:  strncpy(text, "1", 10);      break;
    case  1:  strncpy(text, "10", 10);     break;
    case  2:  strncpy(text, "100", 10);    break;
    case  3:  strncpy(text, "1000", 10);   break;
    default:
        snprintf(text, 20, "1.E%i", n);
    }
}

void AltruistView::drawLimitMap() {
    // draw x-y map of limiting parameter values from 2-dimensional search for parameter limits

    ParameterLoop *loops = main->worker->loops;            // parameter loops
    ParameterLimits *limitList = main->worker->limitList;  // list of search results
    
    int textIndex = 0;                                     // index to textItem
    int squaresIndex = 0;                                  // index to squares
    int i;                                                 // loop counter

    // number of points in limitList
    int limitListNum = main->worker->limitListNum;
    if (limitListNum == main->worker->lastLimitListNum) return;  // no change  
    main->worker->lastLimitListNum = limitListNum;

    if (loops[0].type != 2 || loops[1].type != 1) {        // inner loop must be search, next loop must be sweep
        errorMessage("Wrong graphics view type");  return;
    }

    // parameter ranges for x and y
    float xStart = loops[1].startValue;
    float xEnd = loops[1].endValue;
    float yStart = loops[0].startValue;
    float yEnd = loops[0].endValue;
    bool xDescending = xStart > xEnd;             // x axis direction
    bool yDescending = yStart > yEnd;             // y axis direction
    if (xDescending) {
        float temp = xStart;                      // make sure start < end
        xStart = xEnd;
        xEnd = temp;    
    }
    if (yDescending) {
        float temp = yStart;
        yStart = yEnd;
        yEnd = temp;    
    }
    // linearized parameter ranges for x and y
    float xStartLin = loops[1].p2x(xStart);
    float xEndLin = loops[1].p2x(xEnd);
    float yStartLin = loops[0].p2x(yStart);
    float yEndLin = loops[0].p2x(yEnd);

    // lambda function to return linearized x value to graphics coordinate
    auto xCoord = [xStartLin, xEndLin](float x) {
        return (int)lround((x-xStartLin) * xAxisLength / (xEndLin - xStartLin));
    };

    // lambda function to return linearized y value to graphics coordinate
    auto yCoord = [yStartLin, yEndLin](float y) {
        return yAxisLength - (int)lround((y-yStartLin) * yAxisLength / (yEndLin - yStartLin));
    };

    // draw rectangle for coordinate system
    squares[squaresIndex].setRect(0, 0, xAxisLength, yAxisLength);
    squares[squaresIndex].setBrush(QBrush(QColor(0xFF, 0xFF, 0xFF))); // white
    squares[squaresIndex].setVisible(true);
    squaresIndex++;

    // draw x axis
    int x, y;                                              // coordinates
    int tickSize = 5;                                      // size of tick marks on axes
    float xx, yy;                                          // x and y value
    char text[32];                                         // text string
    char format[16];                                       // format string
    int width, height;                                     // size of text string on screen

    if (loops[1].logarithmic) {
        for (xx = roundf(xStartLin); xx <= xEndLin + 0.01f; xx++) {
            if (loops[1].x2p(xx + 1.E-5) < xStart) continue;
            x = xCoord(xx);
            squares[squaresIndex].setRect(x, yAxisLength, 0, tickSize);
            squares[squaresIndex++].setVisible(true);
            // add label
            int xi = lround(xx);
            if (xi == lround(xStartLin) || xi == lround(xEndLin) || xi == 0 || xStartLin - xEndLin < 6) {
                printLogNumber(text, xi);
                textItems[textIndex].setPlainText(text);
                width = lround(textItems[textIndex].boundingRect().width());
                height = lround(textItems[textIndex].boundingRect().height());
                textItems[textIndex].setPos(QPointF(x-width/2, yAxisLength+height/2));
                textItems[textIndex].setRotation(0);
                textItems[textIndex++].setVisible(true);
            }
        }
    }
    else {  // loops[1] linear
        float unit;
        int firstTick, lastTick, ticsPerLabel;
        formatAxis(loops[1].startValue, loops[1].endValue, 5, &unit, &firstTick, &lastTick, &ticsPerLabel, format);
        for (int ix = firstTick; ix <= lastTick; ix++) {
            bool hasLabel = ix % ticsPerLabel == 0;
            xx = ix * unit;
            x = xCoord(xx); 
            squares[squaresIndex].setRect(x, yAxisLength, 0, tickSize + hasLabel * tickSize);
            squares[squaresIndex++].setVisible(true);
            if (hasLabel) {   // add label
                snprintf(text, sizeof(text), format, xx);
                textItems[textIndex].setPlainText(text);
                width  = lround(textItems[textIndex].boundingRect().width());
                height = lround(textItems[textIndex].boundingRect().height());
                textItems[textIndex].setPos(QPointF(x-width/2, yAxisLength+height/2));
                textItems[textIndex].setRotation(0);
                textItems[textIndex++].setVisible(true);
            }
        }
    }
    const char * xName = sweepParameterList[loops[1].parIndex].name; // name of x axis parameter
    textItems[textIndex].setPlainText(xName);
    width = lround(textItems[textIndex].boundingRect().width());
    height = lround(textItems[textIndex].boundingRect().height());
    textItems[textIndex].setPos(QPointF(xAxisLength-width, yAxisLength+height+height/2));
    textItems[textIndex].setRotation(0);
    textItems[textIndex++].setVisible(true);

    // draw y axis
    if (loops[0].logarithmic) {
        for (yy = roundf(yStartLin); yy <= yEndLin + 0.01f; yy++) {
            if (loops[0].x2p(yy + 1.E-5) < yStart) continue;
            y = yCoord(yy);
            squares[squaresIndex].setRect(-tickSize, y, tickSize, 0);
            squares[squaresIndex++].setVisible(true);
            // add label
            int yi = lround(yy);
            if (yi == lround(yStartLin) || yi == lround(yEndLin) || yi == 0 || yStartLin - yEndLin < 8) {
                printLogNumber(text, yi);
                textItems[textIndex].setPlainText(text);
                width = lround(textItems[textIndex].boundingRect().width());
                height = lround(textItems[textIndex].boundingRect().height());
                textItems[textIndex].setPos(-tickSize - 2 - width, y-height/2);
                textItems[textIndex].setRotation(0);
                textItems[textIndex++].setVisible(true);
            }
        }
    }
    else { // loops[0] linear
        float unit;
        int firstTick, lastTick, ticsPerLabel;
        formatAxis(loops[0].startValue, loops[0].endValue, 10, &unit, &firstTick, &lastTick, &ticsPerLabel, format);
        for (int iy = firstTick; iy <= lastTick; iy++) {
            bool hasLabel = iy % ticsPerLabel == 0;
            yy = iy * unit;
            y = yCoord(yy);
            int tickSize1 = tickSize + hasLabel * tickSize;
            squares[squaresIndex].setRect(-tickSize1, y, tickSize1, 0);
            squares[squaresIndex++].setVisible(true);
            if (hasLabel) {   // add label
                snprintf(text, sizeof(text), format, yy);
                textItems[textIndex].setPlainText(text);
                width = lround(textItems[textIndex].boundingRect().width());
                height = lround(textItems[textIndex].boundingRect().height());
                textItems[textIndex].setPos(-tickSize - 2 - width, y-height/2);
                textItems[textIndex].setRotation(0);
                textItems[textIndex++].setVisible(true);
            }
        }
    }

    const char * yName = sweepParameterList[loops[0].parIndex].name; // name of y axis parameter
    textItems[textIndex].setPlainText(yName);
    int width1 = width;
    width = lround(textItems[textIndex].boundingRect().width());
    height = lround(textItems[textIndex].boundingRect().height());
    textItems[textIndex].setPos(QPointF(-height-width1-tickSize, width));
    textItems[textIndex].setRotation(-90);
    textItems[textIndex++].setVisible(true);

    // polygons to enclose parameter areas
    QPolygonF altruRange, polymorphismRange, egoismRange;
    bool  direction = false;
    int   allResults = 0;
    float limit1;                                          // limit between altruim and polymorphism
    float limit2;                                          // limit between polymorphism and egoism
    float limit1Lin, limit2Lin;                            // limit1 and limit2 linearized if logarithmic scale
    int   y0, y1, y2, y3;                                  // coordinates for y limits
    float sValue;                                          // value of swept parameter
    float xLin;                                            // xValue linearized if logarithmic scale
    int   x0;                                              // x coordinate
    int   xIndex;                                          // index into list of points

    // loop through result points to get lower side points of polygons
    for (xIndex = 0; xIndex < limitListNum; xIndex++) {
        limit1 = limitList[xIndex].limit1;
        limit2 = limitList[xIndex].limit2;
        direction = limitList[xIndex].result & 1;
        allResults |= limitList[xIndex].result;
        limit1Lin = loops[0].p2x(limit1);
        limit2Lin = loops[0].p2x(limit2);
        y0 = yCoord(yStartLin);
        y1 = yCoord(limit1Lin);
        y2 = yCoord(limit2Lin);
        y3 = yCoord(yEndLin);
        sValue = limitList[xIndex].sValue;
        xLin = loops[1].p2x(sValue);
        x0 = xCoord(xLin);

        // add points to polygons
        if (xIndex == 0) { // add left edge
            if (direction) {
                egoismRange << QPointF(x0, y1);
                polymorphismRange << QPointF(x0, y2);
                altruRange << QPointF(x0, y3);
            }
            else {
                altruRange << QPointF(x0, y1);
                polymorphismRange << QPointF(x0, y2);
                egoismRange << QPointF(x0, y3);
            }
        }
        // add lower side
        if (direction) {
            egoismRange << QPointF(x0, y0);
            polymorphismRange << QPointF(x0, y1);
            altruRange << QPointF(x0, y2);
        }
        else {
            altruRange << QPointF(x0, y0);
            polymorphismRange << QPointF(x0, y1);
            egoismRange << QPointF(x0, y2);
        }
    }
    // loop back to get upper side points of polygons
    for (xIndex = limitListNum-1; xIndex >= 0; xIndex--) {
        limit1 = limitList[xIndex].limit1;
        limit2 = limitList[xIndex].limit2;
        direction = limitList[xIndex].result & 1;
        limit1Lin = loops[0].p2x(limit1);
        limit2Lin = loops[0].p2x(limit2);
        y1 = yCoord(limit1Lin);
        y2 = yCoord(limit2Lin);
        y3 = yCoord(yEndLin);
        sValue = limitList[xIndex].sValue;
        xLin = loops[1].p2x(sValue);
        x0 = xCoord(xLin);
        // add right edge and upper side
        if (direction) {
            egoismRange << QPointF(x0, y1);
            polymorphismRange << QPointF(x0, y2);
            altruRange << QPointF(x0, y3);
        }
        else {
            altruRange << QPointF(x0, y1);
            polymorphismRange << QPointF(x0, y2);
            egoismRange << QPointF(x0, y3);
        }
    }

    if (allResults & 4) {
        main->worker->d->simulationResult = resultDied;
    }

    // show polygons
    polygons[0].setPolygon(altruRange);
    polygons[0].setBrush(QBrush(QColor(0xFF,0,0)));  // red
    polygons[0].setVisible(true);

    polygons[1].setPolygon(polymorphismRange);
    polygons[1].setBrush(QBrush(QColor(0,0xFF,0)));  // green
    polygons[1].setVisible(true);

    polygons[2].setPolygon(egoismRange);
    polygons[2].setBrush(QBrush(QColor(0,0,0xFF)));  // blue
    polygons[2].setVisible(true);

    // hide any unused text items
    for (i = textIndex; i < textItemsUsed; i++) {
        textItems[i].setPos(0, 0);
        textItems[i].setRotation(0);
        textItems[i].setVisible(false);    
    }
    textItemsUsed = textIndex;
    // hide any unused squares
    for (i = squaresIndex; i < squaresUsed && i < maxIslandsShown; i++) {
        squares[i].setPos(0, 0);
        squares[i].setVisible(false);    
    }
    squaresUsed = squaresIndex;
    // hide any unused polygons
    for (i = 3; i < polygonsUsed && i < maxPolygons; i++) {
        polygons[i].setPos(0, 0);
        polygons[i].setVisible(false);    
    } 
    polygonsUsed = 3;
}


void AltruistView::drawTerritories() {
    // draw map of territories
    AltruData * const d = &main->d;
    if (d == 0 || d->currentModel == 0) return;  // not initialized yet
    if (squaresUsed) clear();                    // hide unused squares and text items
    int polygonsIndex = 0;                       // index to polygons
    int textIndex = 0;                           // index to text items

    // inspect group descriptor to find offset for size and group count
    int32_t nn = 0;                              // group population
    int32_t ga = 0;                              // gene count for mutant at primary locus. determines hue or first color
    int32_t tArea = 0;                           // area of territory
    float groupProperty = 0;                     // group fitness or other group property
    int onn = -1;                                // offset of population count in group data structure
    int oga0 = 0;                                // offset of locus 0 mutant gene in group data structure
    int oGroupProperty = 0;                      // offset of groupProperty in group data structure
    int oArea = 0;                               // offset to territory area
    int varType = 0;
    int g;      // graphics field type
    int i = 0;  // index into groupFields
    do {    // search through group descriptor list
        g = d->currentModel->groupFields[i].graphics;
        switch (g) {
        case 3:   // field contains population size
            onn = d->currentModel->groupFields[i].offset;
            varType = d->currentModel->groupFields[i].varType; 
            break;
        case 6:   // field contains area of territory
            oArea = d->currentModel->groupFields[i].offset;
            break;
        case 10:  // field contains first gene count used for selecting color
            oga0 = d->currentModel->groupFields[i].offset;
            break;
        case 20:  // field contains group property
            if (d->currentModel->groupFields[i].varType == varFloat) { // float            
                oGroupProperty = d->currentModel->groupFields[i].offset;
            }
            break;
        }
    } while (d->currentModel->groupFields[i++].type != 0);  // group descriptor list ends with 0

    int graphUnit = xAxisLength / d->rowLengthTerri;       // graphics unit
    if (graphUnit < 1) graphUnit = 1;

    Habitat terrain(d);                                    // for manipulating graphics points

    int8_t * group;                                        // current group territory
    int32_t areaMedium = (d->territorySizeMax + d->territorySizeMin) / 2; // limit between big and small territories
    // loops through groups
    for (int sizeCategory = 0; sizeCategory <= 1; sizeCategory++) {
        for (int iGroup = 1; iGroup <= d->nIslands; iGroup++) {
            // get data for this group                        
            group = d->groupData + iGroup * d->groupStructureSize; // point to current group
            // extract data for this group according to structure defined in group descriptors:
            if (varType == varInt32) {
                if (onn >= 0) nn = *(int32_t*)(group + onn);
                if (oga0) ga = *(int32_t*)(group + oga0);
                if (oArea) tArea = *(int32_t*)(group + oArea);
            }
            if (tArea == 0) continue;                      // skip unused record
            // even though enclaves are rare, we want to be able to show them.
            // Therefore, we draw big territories first, then small territories so that
            // it is likely that an enclave is covered by a later drawn surrounding territory
            if (((tArea > areaMedium) ^ sizeCategory) != 0) {
                // define polygon
                QPolygonF shape;
                // find a point on border of territory
                TerriGroup * tGroup = (TerriGroup *)(group);
                TPoint borderPoint = terrain.findBorder2(tGroup); // avoid border to an enclave

                // walk around border of territory to draw a polygon
                terrain.borderWalk(tGroup->id, borderPoint, [graphUnit, &shape](TPoint&p_1) {
                    // do_us lambda function
                    shape << QPointF(p_1.x * graphUnit, p_1.y * graphUnit); // add border point to polygon
                    },
                    [](TPoint&, TPoint&) {});
                shape << QPointF(borderPoint.x * graphUnit, borderPoint.y * graphUnit); // close polygon

                // get color according to fraction of altruism gene
                QColor color(QColor(0xFF, 0xFF, 0xFF));        // white if nn = 0
                if (nn > 0) {                                  // avoid division by 0
                    int hue = (nn * 2 - ga) * 240 / (2 * nn);  // 240 = blue, 0 = red
                    color = QColor::fromHsv(hue, 0xFF, 0xFF);  // get color with hue from blue to red
                }
                if (polygonsIndex >= maxPolygons) break;
                polygons[polygonsIndex].setPolygon(shape);
                polygons[polygonsIndex].setBrush(QBrush(color));
                polygons[polygonsIndex].setVisible(true);
                polygonsIndex++;
            }
        }
    }
    if (polygonsIndex >= maxPolygons) {
        // tell that not all territories can be shown
        textItems[0].setPlainText("Too many territories. Not all can be shown.");
        textItems[0].setPos(QPointF(0, (d->numRows+1) * graphUnit));
        textItems[0].setRotation(0);
        textItems[0].setVisible(true);
        textIndex++;
    }

    // hide any unused polygons
    for (int i = polygonsIndex; i < polygonsUsed && i < maxPolygons; i++) {
        polygons[i].setPolygon(QPolygonF());
        polygons[i].setVisible(false);
    }
    polygonsUsed = polygonsIndex;
    // hide any unused text items
    for (i = textIndex; i < textItemsUsed; i++) {
        textItems[i].setPos(0, 0);
        textItems[i].setRotation(0);
        textItems[i].setVisible(false);
    }
    textItemsUsed = textIndex;
}

// mouse handlers

void AltruistView::mousePressEvent(QMouseEvent * event) {
    // handle mouse click on drawing
    if (event->button() == Qt::LeftButton) {
        QPointF scenePoint = mapToScene(event->pos());     // map to scene coordinates
        int x = (int)lround(scenePoint.x());
        int y = (int)lround(scenePoint.y());
        switch (main->d.graphicsType) {
        case graphicsIslands:
            mouseHandlerIslands(x, y);
            break;
        case graphicsTerritories:
            mouseHandlerTerritories(x, y);
            break;
        case graphicsLimits:
            mouseHandlerLimitMap(x, y);
            break;
        }
    }
}

void AltruistView::mouseHandlerIslands(int x, int y) {
    // handle mouse click on islands graphic
    AltruData * d = &main->d;                    // point to main data
    if (d == 0 || d->currentModel == 0) return;  // not initialized yet
    // number of rows and columns to show
    int32_t numColumns = d->rowLength;
    int32_t numRows = d->numRows;
    int32_t row, column;
    int32_t index = 0;
    int8_t * group;                              // point to current group
    int maxSquare = gridUnit - 4;                // max size of squares in graphics
    if (maxSquare < 2) maxSquare = 2;
    int maxDist = maxSquare * 2;                 // maximum distance from edges of drawing
    // check if point is outside drawing
    if (x < -maxDist || x > numColumns * gridUnit + maxDist) return;
    if (y < -maxDist || y > numRows * gridUnit + maxDist) return;
    // get row and column
    row = (y + gridUnit/2) / gridUnit;
    column = (x + gridUnit/2) / gridUnit;
    if (row < 0) row = 0;  if (row >= numRows) row = numRows - 1;
    if (column < 0) column = 0;  if (column >= numColumns) column = numColumns - 1;
    // get island index
    index = row * numRows + column;
    // make message box
    groupDescriptionMessageBox(index, d);
}

void AltruistView::mouseHandlerTerritories(int x, int y) {
    // handle mouse click on territories graphic
    AltruData * d = &main->d;                    // point to main data
    if (d == 0 || d->currentModel == 0) return;  // not initialized yet
    int graphUnit = xAxisLength / d->rowLengthTerri; // graphics unit
    if (graphUnit < 1) graphUnit = 1;
    // translate x,y to area units
    int xx = (x + graphUnit/2) / graphUnit;
    int yy = (y + graphUnit/2) / graphUnit;
    if (xx <= 0 || xx >= d->rowLengthTerri) return;
    if (yy <= 0 || yy >= d->numRows) return;
    Habitat terrain(d);                          // for manipulating graphics points

    int owner = terrain.owner(TPoint(xx,yy));
    if (owner <= 0 || owner > d->nIslands) return;
    // make message box
    groupDescriptionMessageBox(owner, d);
}

void groupDescriptionMessageBox(int id, AltruData * d) {
    // make a message box showing the contents of a group record
    int8_t * group = d->groupData + id * d->groupStructureSize; // point to current group

    // make text with variable values
    const int textsize = 1024;
    char text[textsize];
    text[0] = 0;
    int textlen = 0;

    // search group field descriptor for variables to show
    int offset;                                      // offset into group data structure
    int g;                                           // graphics field type
    int i = 0;                                       // index into groupFields
    int32_t xi;                                      // integer variable
    double xf;                                       // floating point variable
    do {    // search through group descriptor list
        g = d->currentModel->groupFields[i].graphics; // graphics display info
        if (g == 10 && !d->locusUsed[0]) continue;   // don't show unused locus
        if (g == 11 && !d->locusUsed[1]) continue;   // don't show unused locus
        if (g == 12 && !d->locusUsed[2]) continue;   // don't show unused locus
        if (g == 13 && !d->locusUsed[3]) continue;   // don't show unused locus

        if (g != 0) {
            textlen += snprintf(text + strlen(text), textsize - textlen, "\n%s", d->currentModel->groupFields[i].name);
            offset = d->currentModel->groupFields[i].offset;
            switch (d->currentModel->groupFields[i].varType) {
            case varInt16:
                xi = *(int16_t*)(group + offset);
                textlen += snprintf(text + strlen(text), textsize - textlen, " = %i", xi);
                break;
            case varInt32:
                xi = *(int32_t*)(group + offset);
                textlen += snprintf(text + strlen(text), textsize - textlen, " = %i", xi);
                break;
            case varFloat:
                xf = *(float*)(group + offset);
                textlen += snprintf(text + strlen(text), textsize - textlen, " = %.3G", xf);
                break;
            case varDouble:
                xf = *(double*)(group + offset);
                textlen += snprintf(text + strlen(text), textsize - textlen, " = %.3G", xf);
                break; 
            }        
        }
    } while (d->currentModel->groupFields[i++].type != 0);  // group descriptor list ends with 0 

    // make heading
    char heading[64];
    char const * title = "";
    switch (d->graphicsType) {
    case graphicsIslands:
        title = "Island";
        break;
    case graphicsTerritories:
        title = "Territory";
        break;
    default:
        title = "Group";
        break;
    }
    snprintf(heading, sizeof(heading), "%s %i", title, id);

    // make message box
    messageBox(text, heading);
}

void AltruistView::mouseHandlerLimitMap(int x, int y) {
    // handle mouse click on limits map graphic
    ParameterLoop *loops = main->worker->loops;  // parameter loops

    int maxOutside = 5;                          // max deviation outside coordinate system
    if (x < -maxOutside || x > xAxisLength + maxOutside) return;
    if (y < -maxOutside || y > yAxisLength + maxOutside) return;

    // parameter ranges for x and y
    float xStart = loops[1].startValue;
    float xEnd = loops[1].endValue;
    float yStart = loops[0].startValue;
    float yEnd = loops[0].endValue;
    bool xDescending = xStart > xEnd;             // x axis direction
    bool yDescending = yStart > yEnd;             // y axis direction
    if (xDescending) {
        float temp = xStart;                      // make sure start < end
        xStart = xEnd;
        xEnd = temp;    
    }
    if (yDescending) {
        float temp = yStart;
        yStart = yEnd;
        yEnd = temp;    
    }
    // linearized parameter ranges for x and y
    float xStartLin = loops[1].p2x(xStart);
    float xEndLin = loops[1].p2x(xEnd);
    float yStartLin = loops[0].p2x(yStart);
    float yEndLin = loops[0].p2x(yEnd);

    // convert graphics coordinates to linearized values
    float xLin = x * (xEndLin - xStartLin) / xAxisLength + xStartLin;
    float yLin = (yAxisLength - y) * (yEndLin - yStartLin) / yAxisLength + yStartLin;

    // convert to possibly logarithmic values
    float xVal = loops[1].x2p(xLin);
    float yVal = loops[0].x2p(yLin);

    // get parameter names
    int xParIndex = loops[1].parIndex;
    int yParIndex = loops[0].parIndex;
    const char * xName = sweepParameterList[xParIndex].name;
    const char * yName = sweepParameterList[yParIndex].name;

    // write text        
    const int textsize = 256;
    char text[textsize];
    snprintf(text, textsize, "%s = %.4G\n%s = %.4G", xName, xVal, yName, yVal);

    // make message box
    messageBox(text, "Parameter map");
}
