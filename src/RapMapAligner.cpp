#include <limits>
#include <iostream>
#include "RapMapAligner.hpp"

RapMapAligner::RapMapAligner(int matchIn, int misMatchIn, int gapExtendIn, int gapStartIn,
                             bool gapBefore, bool gapAfter) :
        match(matchIn),
        misMatch(misMatchIn),
        gapExtend(gapExtendIn),
        gapStart(gapStartIn),
        freeGapsBeforeRead(gapBefore),
        freeGapsAfterRead(gapAfter) {
    // init
    tm[0][0] = TraceBack::END;
    tx[0][0] = TraceBack::END;
    ty[0][0] = TraceBack::END;
    m[0][0] = 0;
    x[0][0] = 0;
    y[0][0] = 0;
    for (int i = 1; i <= 255; ++i) {
        m[i][0] = std::numeric_limits<short>::min();
        x[i][0] = std::numeric_limits<short>::min();
        if (freeGapsBeforeRead) {
            y[i][0] = 0;
            ty[i][0] = TraceBack::END;
        } else {
            y[i][0] = gapStart + (i - 1) * gapExtend;
            ty[i][0] = TraceBack::Y;
        }
    }
    for (int j = 1; j <= 255; ++j) {
        m[0][j] = std::numeric_limits<short>::min();
        x[0][j] = gapStart + (j - 1) * gapExtend;
        tx[0][j] = TraceBack::X;
        y[0][j] = std::numeric_limits<short>::min();
    }
}

int RapMapAligner::align(std::string& ref, size_t refStart, size_t refLen,
                         std::string& read, size_t readStart, size_t readLen) {
    int max, xTemp, yTemp;
    TraceBack trace;

    for (int i = 1; i <= refLen; ++i) {
        for (int j = 1; j <= readLen; ++j) {
            // Update m
            max = m[i - 1][j - 1];
            trace = TraceBack::M;
            xTemp = x[i - 1][j - 1];
            yTemp = y[i - 1][j - 1];
            if (max < xTemp) {
                max = xTemp;
                trace = TraceBack::X;
            }
            if (max < yTemp) {
                max = yTemp;
                trace = TraceBack::Y;
            }
            m[i][j] = score(ref[refStart + i - 1], read[readStart + j - 1]) + max;
            tm[i][j] = trace;
            // Update x
            max = gapStart + m[i][j - 1];
            trace = TraceBack::M;
            xTemp = gapExtend + x[i][j - 1];
            yTemp = gapStart + y[i][j - 1];
            if (max < xTemp) {
                max = xTemp;
                trace = TraceBack::X;
            }
            if (max < yTemp) {
                max = yTemp;
                trace = TraceBack::Y;
            }
            x[i][j] = max;
            tx[i][j] = trace;
            // Update y
            max = gapStart + m[i - 1][j];
            trace = TraceBack::M;
            xTemp = gapStart + x[i - 1][j];
            yTemp = gapExtend + y[i - 1][j];
            if (max < xTemp) {
                max = xTemp;
                trace = TraceBack::X;
            }
            if (max < yTemp) {
                max = yTemp;
                trace = TraceBack::Y;
            }
            y[i][j] = max;
            ty[i][j] = trace;
        }
    }

    // Find maximum total score
    size_t maxI = refLen;
    max = m[maxI][readLen];
    trace = TraceBack::M;
    xTemp = x[maxI][readLen];
    yTemp = y[maxI][readLen];
    if (max < xTemp) {
        max = xTemp;
        trace = TraceBack::X;
    }
    if (max < yTemp) {
        max = yTemp;
        trace = TraceBack::Y;
    }

    // If free gaps are allowed after the read
    if (freeGapsAfterRead) {
        for (size_t i = 1; i < refLen; ++i) {
            int mtemp = m[i][readLen];
            xTemp = x[i][readLen];
            yTemp = y[i][readLen];
            if (max < mtemp) {
                max = mtemp;
                maxI = i;
                trace = TraceBack::M;
            }
            if (max < xTemp) {
                max = xTemp;
                maxI = i;
                trace = TraceBack::X;
            }
            if (max < yTemp) {
                max = yTemp;
                maxI = i;
                trace = TraceBack::Y;
            }
        }
    }

    // Traceback
    CigarString cigar;
    bool tracing = true;
    size_t i = maxI, j = readLen;
    std::string cigarString;
    while (tracing and (i > 0 or j > 0)) {
        if (freeGapsBeforeRead and j == 0)
            break;
        switch(trace) {
            case TraceBack::M:
                cigar.push(CigarOp::M);
                trace = tm[i][j];
                --i;
                --j;
                break;
            case TraceBack::X:
                cigar.push(CigarOp::I);
                trace = tx[i][j];
                --j;
                break;
            case TraceBack::Y:
                cigar.push(CigarOp::D);
                trace = ty[i][j];
                --i;
                break;
            default:
                tracing = false;
        }
    }
    cigar.toString(cigarString);
    std::cerr << cigarString << std::endl;
    return max;
}

