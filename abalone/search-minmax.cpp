/**
 * Very simple example strategy:
 * Search all possible positions reachable via one move,
 * and return the move leading to best position
 *
 * (c) 2006, Josef Weidendorfer
 */

#include "cstdio"
#include "search.h"
#include "board.h"
#include "eval.h"

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
class MinMaxStrategy: public SearchStrategy
{
 public:
    // Defines the name of the strategy
    MinMaxStrategy(): SearchStrategy("MinMax") {}

    // Factory method: just return a new instance of this class
    SearchStrategy* clone() { return new MinMaxStrategy(); }

 private:

    /**
     * Implementation of the strategy.
     */
    void searchBestMove();
    int minimax(int depth);
};


void MinMaxStrategy::searchBestMove()
{

    _maxDepth = 1;
    int eval;
    eval = minimax(0);
    printf("best evaluation = %d\n" , eval);
}


int MinMaxStrategy::minimax(int depth)
{
    if (depth == _maxDepth) {
        return -evaluate();
    }
    Move m;
    MoveList list;
    generateMoves(list);
    int bestEval;
    bestEval = -maxEvaluation();
    int eval;
    while(list.getNext(m)) {
        playMove(m);
        eval = -minimax(depth + 1);
        takeBack();
        printf("Self %d: %s=%d\n",depth,m.name(),eval);
        if(eval > bestEval){
            bestEval = eval;
            foundBestMove(depth,m,bestEval);
        }
    }
    printf("Self %d: best=%d\n",depth,bestEval);
    finishedNode(depth,0);
    return bestEval;
}

// register ourselve as a search strategy
MinMaxStrategy minMaxStrategy;
