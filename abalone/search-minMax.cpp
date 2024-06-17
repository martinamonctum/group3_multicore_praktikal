/**
 * Very simple example strategy:
 * Search all possible positions reachable via one move,
 * and return the move leading to best position
 *
 * (c) 2006, Josef Weidendorfer
 */
#include <iostream>
#include "search.h"
#include "board.h"
#include "eval.h"


#define max(i,j) ((i>j)?i:j)
#define min(i,j) ((i<j)?i:j)

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
class MinMaxAlgorithm: public SearchStrategy
{
 public:
    // Defines the name of the strategy
    MinMaxAlgorithm(): SearchStrategy("minMax") {}

    // Factory method: just return a new instance of this class
    SearchStrategy* clone() { return new MinMaxAlgorithm(); }

 private:

    /**
     * Implementation of the strategy.
     */
    int minimax(int depth, bool isMaximizing);
    void searchBestMove();
};


void MinMaxAlgorithm::searchBestMove()
{
    _maxDepth = 2;
    int eval;
    eval = minimax(0 ,true);
    printf("best evaluation = %d\n" , eval);
}

int MinMaxAlgorithm::minimax(int depth, bool isMaximizing)
{
    if (depth == _maxDepth) {
        return evaluate();
    }
    Move m;
    MoveList list;
    generateMoves(list);
    int bestEval;
    if (isMaximizing) {
        bestEval = minEvaluation();
        int eval;
        while(list.getNext(m)) {
            playMove(m);
            eval = -minimax(depth + 1, false);
            takeBack();
            if(eval > bestEval){
                bestEval = eval;
                foundBestMove(depth,m,bestEval);
            }
        }
        finishedNode(depth,0);
        return bestEval;
    } else {
        bestEval = maxEvaluation();
        int eval;
        while(list.getNext(m)) {
            playMove(m);
            eval = minimax(depth + 1, true);
            takeBack();
            if(eval < bestEval){
                bestEval = eval;
                foundBestMove(depth,m,bestEval);
            }
        }
        finishedNode(depth,0);
        return bestEval;
    }
}

// register ourselve as a search strategy
MinMaxAlgorithm minmaxAlgorithm;
