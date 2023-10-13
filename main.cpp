#include <stdio.h>
#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <cassert>
#include "easygl/graphics.h"

# define ll unsigned long long
# define SOUTH 0
# define EAST 1
# define NORTH 2
# define WEST 3

using namespace std;

class PinLoc {
    public:
        char b;
        int x, y, p;
        PinLoc(char _b, int _x, int _y, int _p) {
            b = _b;
            x = _x;
            y = _y;
            p = _p;
        }
        string toStr() {
            string strB = ""; strB.push_back(b);
            return to_string(x) + "_" + to_string(y) + "_" + strB + "_" + to_string(p);
        }
};
PinLoc str2Loc(string str) {
    vector<string> tokens;
    istringstream tokenStream(str);
    string token;
    while (getline(tokenStream, token, '_')) {
        tokens.push_back(token);
    }
    char b = tokens[2][0];
    int x = stoi(tokens[0]);
    int y = stoi(tokens[1]);
    int p = stoi(tokens[3]);
    return PinLoc(b, x, y, p);
}

class Segment {
    public:
        char orient;
        int l;
        int r;
        int c;
        Segment () {
            orient = 'u';
            l = 0;
            r = 0;
            c = 0;
        };
        Segment (char _orient, int _l, int _r, int _c) {
            orient = _orient;
            l = _l;
            r = _r;
            c = _c;
        };
        bool operator==(const Segment& other) const {
            return orient == other.orient && l == other.l && r == other.r && c == other.c;
        }
        bool operator-=(const Segment& other) const {
            return r == other.r && c == other.c && orient == other.orient;
        }
        bool operator!=(const Segment& other) const {
            return !(orient == other.orient && l == other.l && r == other.r && c == other.c);
        }
};

////////////////
// Test Utils //
////////////////

class Test {
    public:
        bool isDense;
        int N;
        int W;
        vector<PinLoc> srcs, sinks;
        Test(bool _isDense, int _N, int _W, vector<PinLoc> _srcs, vector<PinLoc> _sinks) {
            isDense = _isDense;
            N = _N;
            W = _W;
            srcs = _srcs;
            sinks = _sinks;
        }
};

Test genTest (string filePath) {
    bool isDense = filePath.find("dense") != string::npos;
    int N, W;
    vector<PinLoc> srcs, sinks;
    ifstream file(filePath);
    string line;
    // get N, W
    getline(file, line);
    N = stoi(line);
    getline(file, line);
    W = stoi(line);
    // get pins & loads
    while (1) {
        if (isDense) {
            int x0, y0, p0, x1, y1, p1;
            char b0, b1;
            file >> x0;
            file >> y0;
            file >> b0;
            file >> p0;
            file >> x1;
            file >> y1;
            file >> b1;
            file >> p1;
            if (x0 == -1) {
                break;
            }
            srcs.push_back(PinLoc(b0, x0, y0, p0));
            sinks.push_back(PinLoc(b1, x1, y1, p1));
        } else {
            int x0, y0, p0, x1, y1, p1;
            file >> x0;
            file >> y0;
            file >> p0;
            file >> x1;
            file >> y1;
            file >> p1;
            if (x0 == -1) {
                break;
            }
            srcs.push_back(PinLoc('a', x0, y0, p0));
            sinks.push_back(PinLoc('a', x1, y1, p1));
        }
    }
    file.close();
    return Test(isDense, N, W, srcs, sinks);
}

//////////////////
// Gui Utils    //
//////////////////
#define TRACK_BLOCK_RATIO 3
#define PIN_BLOCK_RATIO 2.6667
#define LB_A_OFFSET_BLOCK_RATIO 0.5
#define LB_B_OFFSET_BLOCK_RATIO 0.125
class SB {
    public:
        float x1, y1; // up left
        float x2, y2; // bottom right
        SB () {
            x1 = 0;
            x2 = 0;
            y1 = 0;
            y2 = 0;
        }
        SB (float _x1, float _y1, float blockLen) {
            x1 = _x1;
            y1 = _y1;
            x2 = x1 + blockLen;
            y2 = y1 + blockLen;
        }
        void draw() {}
};

class LB {
    public:
        float x1, y1; // up left
        float x2, y2; // bottom right
        vector<float> px1, py1, px2, py2; // pin locations
        LB () {
            x1 = 0;
            x2 = 0;
            y1 = 0;
            y2 = 0;
        }
        LB (float _x1, float _y1, float blockLen, float pinLen) {
            // Left top
            x1 = _x1; 
            y1 = _y1;
            // Right bottom
            x2 = x1 + blockLen; 
            y2 = y1 + blockLen;
            // Pins
            px1.resize(4); 
            py1.resize(4); 
            px2.resize(4); 
            py2.resize(4);
            // Pin 1
            px1[SOUTH] = (x2 - x1) / 4 * 3 + x1; 
            py1[SOUTH] = y2; 
            px2[SOUTH] = px1[SOUTH]; 
            py2[SOUTH] = py1[SOUTH] + pinLen;
            // Pin 2
            px1[EAST]  = x2;
            py1[EAST]  = (y2 - y1) / 2 + y1;
            px2[EAST]  = px1[EAST] + pinLen;
            py2[EAST]  = py1[EAST];
            // Pin 3
            px1[NORTH] = (x2 - x1) / 2 + x1;
            py1[NORTH] = y1;
            px2[NORTH] = px1[NORTH];
            py2[NORTH] = py1[NORTH] - pinLen;
            // Pin 4
            px1[WEST]  = x1;
            py1[WEST]  = (y2 - y1) / 4 * 3 + y1;
            px2[WEST]  = px1[WEST] - pinLen;
            py2[WEST]  = py1[WEST];
        }
        void draw() {}
};

class Track {
    public:
        float x1, y1; // left/up
        float x2, y2; // right/down
        Track () {
            x1 = 0;
            y1 = 0;
            x2 = 0;
            y2 = 0;
        }
        Track (float _x1, float _y1, char orient, float trackLen) {
            x1 = _x1;
            y1 = _y1;
            x2 = x1 + (orient == 'h' ? trackLen : 0);
            y2 = y1 + (orient == 'v' ? trackLen : 0);
        }
        void draw() {} //TODO: 
};

class Gui {
    public:
        bool isDense;
        int N, W;
        float blockLen;
        float trackLen;
        float pinLen;
        vector<vector<SB>> switchBlocks; 
        vector<vector<map<char, LB>>> logicBlocks;
        map<char, vector<vector<Track>>> tracks;
        Gui(bool _isDense, int _N, int _W) {
            isDense = _isDense;
            N = _N;
            W = _W;
            switchBlocks = vector<vector<SB>>(N+1, vector<SB>(N+1));
            logicBlocks  = vector<vector<map<char, LB>>>(N, vector<map<char, LB>>(N));
            tracks['h']  = vector<vector<Track>>(N, vector<Track>(N+1));
            tracks['v']  = vector<vector<Track>>(N+1, vector<Track>(N));
        }
        void init(float X, float Y, float gridLen) {

            // Set length of tracks and pins.
            // Set Width of blocks
            blockLen = gridLen / (N + TRACK_BLOCK_RATIO * (N-1)); // (N * blockLen + (N-1) * trackLen = gridLen)
            trackLen = TRACK_BLOCK_RATIO * blockLen;
            pinLen   = PIN_BLOCK_RATIO * blockLen;
            
            // Set positions for each switch block on the screen
            for (int r=0; r<N+1; ++r) {
                for (int c=0; c<N+1; ++c) {
                    float x1 = X + c * (blockLen + trackLen);
                    float y1 = Y + r * (blockLen + trackLen);
                    switchBlocks[r][c] = SB(x1, y1, blockLen);
                }
            }

            // Set positions for each Logic block and their pins on the screen
            for (int r=0; r<N; ++r) {
                for (int c=0; c<N; ++c) {
                    SB sb = switchBlocks[r][c];
                    float x1 = sb.x2 + LB_A_OFFSET_BLOCK_RATIO * blockLen;
                    float y1 = sb.y2;
                    logicBlocks[r][c]['a'] = LB(x1, y1, blockLen, pinLen);
                    if (isDense) {
                        x1 += blockLen + LB_B_OFFSET_BLOCK_RATIO * blockLen;
                        y1 += blockLen;
                        logicBlocks[r][c]['b'] = LB(x1, y1, blockLen, pinLen);
                    }
                }
            }

            // Set positions for each horizontal track on the screen
            float v = blockLen / (W+1);
            for (int r=0; r<N+1; ++r) {
                for (int c=0; c<N; ++c) {
                    for (int l=0; l<W; ++l) {
                        SB sb = switchBlocks[r][c];
                        float x1 = sb.x2;
                        float y1 = sb.y1 + v*(l+1);
                        tracks['h'][r][c] = Track(x1, y1, 'h', trackLen);
                    }
                }
            }

            // Set positions for each vertical track on the screen
            float h = blockLen / (W+1);
            for (int r=0; r<N; ++r) {
                for (int c=0; c<N+1; ++c) {
                    for (int l=0; l<W; ++l) {
                        SB sb = switchBlocks[r][c];
                        float x1 = sb.x1 + h*(l+1);
                        float y1 = sb.y2;
                        tracks['v'][r][c] = Track(x1, y1, 'v', trackLen);
                    }
                }
            }

        }
};

//////////////////
// Solver Utils //
//////////////////
struct PathNode {
    Segment segment;
    ll totalCost;
    PathNode* parent; // For back tracking
    vector<PathNode*> children; // Consider multi-load paths
    bool isLocated; // Set to true when the node is located in a path through back trace
    PathNode(Segment _segment, ll _totalCost) {
        segment = _segment;
        totalCost = _totalCost;
        parent = NULL;
        children = {};
        isLocated = 0;
    }
};
class Solver {
    public:
        int N, W;
        vector<vector<vector<vector<Segment>>>> layout;
        map<char, vector<vector<vector<ll>>>> costMatrix;
        map<string, vector<PinLoc>> src2sinks;

        // Initialize the layout, compute the src2sink map
        Solver (Test test) {
            N = test.N;
            W = test.W;
            vector<PinLoc> srcs = test.srcs;
            vector<PinLoc> sinks = test.sinks;

            // Get Number of Switch Blocks per row/column
            int numSB = N + 1;

            // Initialize the layout by appending segments to each switch block
            // following the order: South -> East -> North -> West
            layout.resize(W, vector<vector<vector<Segment>>>(
                numSB, vector<vector<Segment>>(
                    numSB, vector<Segment>(4)
                )
            ));
            Segment unExistingSeg = Segment('u', -1, -1, -1);
            for (int l=0; l<W; l++) {
                for (int i=0; i<numSB; i++) {
                    for (int j=0; j<numSB; j++) {
                        // South
                        if (i < numSB - 1) {
                            layout[l][i][j][SOUTH] = Segment('v', l, i, j);
                        } else {
                            layout[l][i][j][SOUTH] = unExistingSeg;
                        }
                        // East
                        if (j < numSB - 1) {
                            layout[l][i][j][EAST] = Segment('h', l, i, j);
                        } else {
                            layout[l][i][j][EAST] = unExistingSeg;
                        }
                        // North
                        if (i > 0) {
                            layout[l][i][j][NORTH] = Segment('v', l, i-1, j);
                        } else {
                            layout[l][i][j][NORTH] = unExistingSeg;
                        }
                        // West
                        if (j > 0) {
                            layout[l][i][j][WEST] = Segment('h', l, i, j-1);
                        } else {
                            layout[l][i][j][WEST] = unExistingSeg;
                        }
                    }
                }
            }

            // Initialize the cost for each segment to 1
            costMatrix['v'] = vector<vector<vector<ll>>>(
                W, vector<vector<ll>>(
                    numSB-1, vector<ll>(
                        numSB, 1
                    )
                )
            );
            costMatrix['h'] = vector<vector<vector<ll>>>(
                W, vector<vector<ll>>(
                    numSB, vector<ll>(
                        numSB-1, 1
                    )
                )
            );

            // Condense source/sink pairs as there can be 1-to-N routings
            assert(srcs.size() == sinks.size());
            for (int i=0; i<srcs.size(); i++) {
                string srcStr = srcs[i].toStr();
                src2sinks[srcStr].push_back(sinks[i]);
            }
        }

        // Path finder
        vector<Segment>   srcPin2Segs (PinLoc output) {
            int x = output.x, y = output.y, p = output.p;
            assert(p == 4);
            int r = N - 1 - y;
            int c = x;
            int i = r;
            int j = c;
            vector<Segment> startSegs;
            for (int l=0; l<W; l++) {
                startSegs.push_back(layout[l][i][j][SOUTH]);
            }
            return startSegs;
        }
        vector<Segment>   sinkPin2Segs(PinLoc input) {
            int x = input.x, y = input.y, p = input.p;
            assert(p != 4);
            int r = N - 1 - y;
            int c = x;
            int dir = p - 1;
            int i = -1, j = -1;
            vector<Segment> endSegs;
            switch (dir) {
                case SOUTH:
                    i = r + 1;
                    j = c;
                    for (int l=0; l<W; l++) {
                        endSegs.push_back(layout[l][i][j][EAST]);
                    }
                    break;
                case EAST:
                    i = r;
                    j = c + 1;
                    for (int l=0; l<W; l++) {
                        endSegs.push_back(layout[l][i][j][SOUTH]);
                    }
                    break;
                case NORTH:
                    i = r;
                    j = c;
                    for (int l=0; l<W; l++) {
                        endSegs.push_back(layout[l][i][j][EAST]);
                    }
                    break;
                default:
                    assert(0);
            }
            return endSegs;
        }
        vector<PathNode*> genNeighbors(PathNode* cur) {
            Segment seg = cur->segment;
            ll totalCost = cur->totalCost;
            int r = seg.r, c = seg.c, l = seg.l;

            // Get the two switch blocks at the ends of the segment (r1,c1); (r2,c2)
            int r1 = r, c1 = c;
            int r2, c2;
            switch(seg.orient) {
                case 'v':
                    r2 = r+1;
                    c2 = c;
                    break;
                case 'h':
                    r2 = r;
                    c2 = c+1;
                    break;
                default:
                    assert(0);
            }

            // Get neighboring segments
            if (seg.orient == 'v') {
                Segment upEast    = layout[l][r1][c1][EAST];
                Segment upNorth   = layout[l][r1][c1][NORTH];
                Segment upWest    = layout[l][r1][c1][WEST];
                Segment downEast  = layout[l][r2][c2][EAST];
                Segment downSouth = layout[l][r2][c2][SOUTH];
                Segment downWest  = layout[l][r2][c2][WEST];
                PathNode* upEastNode    = new PathNode(upEast,    totalCost + costMatrix['h'][l][upEast.r][upEast.c]);
                PathNode* upNorthNode   = new PathNode(upNorth,   totalCost + costMatrix['v'][l][upNorth.r][upNorth.c]);
                PathNode* upWestNode    = new PathNode(upWest,    totalCost + costMatrix['h'][l][upWest.r][upWest.c]);
                PathNode* downEastNode  = new PathNode(downEast,  totalCost + costMatrix['h'][l][downEast.r][downEast.c]);
                PathNode* downSouthNode = new PathNode(downSouth, totalCost + costMatrix['v'][l][downSouth.r][downSouth.c]);
                PathNode* downWestNode  = new PathNode(downWest,  totalCost + costMatrix['h'][l][downWest.r][downWest.c]);
                return {upEastNode, upNorthNode, upWestNode, downEastNode, downSouthNode, downWestNode};
            }
            if (seg.orient == 'h') {
                Segment leftNorth  = layout[l][r1][c1][NORTH];
                Segment leftWest   = layout[l][r1][c1][WEST];
                Segment leftSouth  = layout[l][r1][c1][SOUTH];
                Segment rightSouth = layout[l][r2][c2][SOUTH];
                Segment rightEast  = layout[l][r2][c2][EAST];
                Segment rightNorth = layout[l][r2][c2][NORTH];
                PathNode* leftNorthNode  = new PathNode(leftNorth,  totalCost + costMatrix['v'][l][leftNorth.r][leftNorth.c]);
                PathNode* leftWestNode   = new PathNode(leftWest,   totalCost + costMatrix['h'][l][leftWest.r][leftWest.c]);
                PathNode* leftSouthNode  = new PathNode(leftSouth,  totalCost + costMatrix['v'][l][leftSouth.r][leftSouth.c]);
                PathNode* rightSouthNode = new PathNode(rightSouth, totalCost + costMatrix['v'][l][rightSouth.r][rightSouth.c]);
                PathNode* rightEastNode  = new PathNode(rightEast,  totalCost + costMatrix['h'][l][rightEast.r][rightEast.c]);
                PathNode* rightNorthNode = new PathNode(rightNorth, totalCost + costMatrix['v'][l][rightNorth.r][rightNorth.c]);
                return {leftNorthNode, leftWestNode, leftSouthNode, rightSouthNode, rightEastNode, rightNorthNode};
            }
            assert(0);
            return {};
        }
        vector<PathNode*> dijkstra    (PinLoc output, const vector<PinLoc>& inputs) {
            
            // Declare a min heap
            auto cmp = [&](const PathNode* n1, const PathNode* n2) {
                return n1->totalCost > n2->totalCost;
            };
            priority_queue<PathNode*, vector<PathNode*>, decltype(cmp)> pq(cmp);

            // Keep track of the minimum cost to each segment for the current route
            map<char, vector<vector<vector<ll>>>> totalCostMatrixForSingleRoute;
            int numSB = N - 1;
            totalCostMatrixForSingleRoute['v'] = vector<vector<vector<ll>>>(
                W, vector<vector<ll>>(
                    numSB-1, vector<ll>(
                        numSB, LLONG_MAX
                    )
                )
            );
            totalCostMatrixForSingleRoute['h'] = vector<vector<vector<ll>>>(
                W, vector<vector<ll>>(
                    numSB, vector<ll>(
                        numSB-1, LLONG_MAX
                    )
                )
            );

            // Get start and multi end segments
            int numEnds = inputs.size();
            vector<Segment> startSegs = srcPin2Segs(output);
            vector<Segment> multiEndSegs(numEnds);
            for (int i=0; i<numEnds; i++) {
                multiEndSegs[i] = sinkPin2Segs(inputs[i])[0]; // Don't care about which segment it is in the track
            }
            for (Segment seg : startSegs) {
                PathNode* startNode = new PathNode(seg, costMatrix[seg.orient][seg.l][seg.r][seg.c]);
                pq.push(startNode);
            }

            // Start dijkstra as the path finder algorithm
            int numEndsRemaining = numEnds;
            vector<PathNode*> endNodes(numEnds, NULL);
            while (!pq.empty()) {
                PathNode* cur = pq.top(); pq.pop();
                Segment curSeg = cur->segment;

                // Check if the current node is the end
                // TODO: Can input pins share the same segment?
                for (int i=0; i<numEnds; i++) {
                    Segment endSeg = multiEndSegs[i];
                    if ((endNodes[i] == NULL) && (endSeg -= curSeg)) {
                        endNodes[i] = cur;
                        numEndsRemaining --;
                        if (numEndsRemaining == 0) {
                            break;
                        }
                    }
                }
                if (numEndsRemaining == 0) {
                    break;
                }
                

                // Generate neighbors
                vector<PathNode*> neis = genNeighbors(cur);
                for (PathNode* nei : neis) {
                    Segment neiSeg = nei->segment;
                    ll neiTotalCost = nei->totalCost;
                    if (neiSeg.orient == 'u') continue; // The segment is out of bound. Skip it
                    if (totalCostMatrixForSingleRoute[neiSeg.orient][neiSeg.l][neiSeg.r][neiSeg.c] > neiTotalCost) {
                        // Push the neighbor node into the heap only if the new total cost is smaller
                        totalCostMatrixForSingleRoute[neiSeg.orient][neiSeg.l][neiSeg.r][neiSeg.c] = neiTotalCost;
                        nei->parent = cur;
                        pq.push(nei);
                    }
                }
            }

            // Back trace the path or tree
            vector<PathNode*> roots;
            for (PathNode* end : endNodes) {
                PathNode* cur = end;
                while ((!cur->isLocated) && !(cur->segment -= startSegs[0])) { // Don't care about which segment it is in the track
                    PathNode* par = cur->parent;
                    par->children.push_back(cur);
                    cur->isLocated = 1;
                    Segment seg = cur->segment;
                    costMatrix[seg.orient][seg.l][seg.r][seg.c] = 123456; // TODO: Update cost calculations
                    cur = par; 
                }
                if (!cur->isLocated) {
                    roots.push_back(cur);
                    cur->isLocated = 1;
                    Segment seg = cur->segment;
                    costMatrix[seg.orient][seg.l][seg.r][seg.c] = 123456; // TODO: Update cost calculations
                }
            }

            return roots;
        }

        // Route all the paths
        void dfs(vector<vector<Segment>>& paths, vector<Segment>& path, PathNode* cur) {
            path.push_back(cur->segment);
            if (cur->children.size() == 0) {
                paths.push_back(path);
            }
            for (PathNode* child : cur->children) {
                dfs(paths, path, child);
            }
            path.pop_back();
        }
        vector<vector<Segment>> routeAllPaths () {
            vector<vector<Segment>> paths;
            // TODO: implement negotiation-based approach
            for (const pair<string, vector<PinLoc>>& p : src2sinks) {
                PinLoc src = str2Loc(p.first);
                vector<PinLoc> sinks = p.second;
                vector<PathNode*> roots = dijkstra(src, sinks);
                for (PathNode* root : roots) {
                    vector<Segment> path = {};
                    dfs(paths, path, root);
                }
            }
            return paths;
        }
};

int main(int argc, char* argv[]) {
    
    // Get all tests
    vector<string> testPaths;
    for (int i=1; i<argc; ++i) {
        testPaths.push_back(argv[i]);
    }
    vector<Test> tests;
    for (string testPath : testPaths) {
        tests.push_back(genTest(testPath));
    }

    // Solve all tests
    for (Test test : tests) {

        // Initialize the gui

        // Initialize the solver

    }

    return 0;
}