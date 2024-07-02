/**
 * Very simple example strategy:
 * Search all possible positions reachable via one move,
 * and return the move leading to best position
 *
 * (c) 2006, Josef Weidendorfer
 */

#include <chrono>
#include <iostream>
#include "cstdio"
#include "search.h"
#include "board.h"
#include "eval.h"

#define max(x,y) (x>y?x:y)

/**
 * To create your own search strategy:
 * - copy this file into another one,
 * - change the class name one the name given in constructor,
 * - adjust clone() to return an instance of your class
 * - adjust last line of this file to create a global instance
 *   of your class
 * - adjust the Makefile to include your class in SEARCH_OBJS
 * - implement searchBestMove()
 *
 * Advises for implementation of searchBestMove():
 * - call foundBestMove() when finding a best move since search start
 * - call finishedNode() when finishing evaluation of a tree node
 * - Use _maxDepth for strength level (maximal level searched in tree)
 */
class AlphaBetaStrategy: public SearchStrategy
{
 public:
    // Defines the name of the strategy
    AlphaBetaStrategy(): SearchStrategy("AlphaBeta") {}

    // Factory method: just return a new instance of this class
    SearchStrategy* clone() { return new AlphaBetaStrategy(); }

 private:

    /**
     * Implementation of the strategy.
     */
    void searchBestMove();
    int alphaBeta(int depth, int alpha, int beta);

    int totalEvaluations = 0;
    std::chrono::duration<double> total_time = std::chrono::duration<double>(0);
};


void AlphaBetaStrategy::searchBestMove()
{
    _maxDepth = 2;
    int eval;
    auto start = std::chrono::steady_clock::now();
    eval = alphaBeta(0, minEvaluation(), maxEvaluation());
    auto end = std::chrono::steady_clock::now();
    total_time += end - start;
    printf("best evaluation = %d\n", eval);
    printf("Time: %f\n", totalEvaluations / total_time.count());
}


int AlphaBetaStrategy::alphaBeta(int depth, int alpha, int beta)
{
    if (depth == _maxDepth) {
        totalEvaluations++;
        return -evaluate();
    }

    Move m;
    MoveList list;
    generateMoves(list);
    int bestEval = minEvaluation();
    int eval;

    while (list.getNext(m)) {
        playMove(m);
        eval = -alphaBeta(depth + 1, -beta, -alpha);
        takeBack();
        
        if (eval > bestEval) {
            bestEval = eval;
            foundBestMove(depth, m, bestEval);
        }
        
        alpha = max(alpha, eval);
        if (alpha >= beta) {
            break;  // Beta cut-off
        }
    }

    finishedNode(depth, 0);
    return bestEval;
}

// register ourselve as a search strategy
AlphaBetaStrategy alphaBetaStrategy;