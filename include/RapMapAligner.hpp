#ifndef __RAPMAP_RAPMAPALIGNER_HPP__
#define __RAPMAP_RAPMAPALIGNER_HPP__

#include <string>
#include <deque>

enum class TraceBack : uint8_t {
    END = 0,
    M,
    X,
    Y
};

enum class CigarOp : char {
    M = 'M',  // alignment match (can be a sequence match or mismatch)
    I = 'I',  // insertion to the reference
    D = 'D',  // deletion from the reference
    N = 'N',  // skipped region from the reference
    S = 'S',  // soft clipping (clipped sequences present in SEQ)
    H = 'H',  // hard clipping (clipped sequences NOT present in SEQ)
    P = 'P',  // padding (silent deletion from padded reference)
    EQ = '=', // sequence match
    X = 'X'   // sequence mismatch
};

class CigarElement {
public:
    int count;
    CigarOp op;

    CigarElement() : count(0), op(CigarOp::M) { }

    CigarElement(int count, CigarOp op) : count(count), op(op) { }
};

class CigarString {
public:
    std::deque<CigarElement> cigar;

    CigarString() : cigar() { }

    inline void push(CigarOp op) {
        if (cigar.size() > 0 and cigar.front().op == op) {
            cigar.front().count++;
        } else {
            cigar.emplace_front(1, op);
        }
    };

    inline void toString(std::string& out) {
        for (auto i = 0; i < cigar.size(); ++i) {
            out += std::to_string(cigar[i].count);
            out += static_cast<char>(cigar[i].op);
        }
    }
};


class RapMapAligner {
public:
    int m[256][256];
    int x[256][256];
    int y[256][256];
    TraceBack tm[256][256];
    TraceBack tx[256][256];
    TraceBack ty[256][256];

    int match;
    int misMatch;
    int gapStart;
    int gapExtend;
    bool freeGapsBeforeRead;
    bool freeGapsAfterRead;

    RapMapAligner(int matchIn, int misMatchIn, int gapExtendIn, int gapStartIn,
                  bool gapBefore, bool gapAfter);

    int align(std::string& ref, size_t refStart, size_t refLen,
              std::string& read, size_t readStart, size_t readLen);

    inline int score(char a, char b) {
        return a == b ? match : misMatch;
    }
};


#endif //__RAPMAP_RAPMAPALIGNER_HPP__
