/**
 * Very simple example strategy:
 * Search all possible positions reachable via one move,
 * and return the move leading to best position
 *
 * (c) 2006, Josef Weidendorfer
 */
#include <iostream>
#include <chrono>
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
    int totalEvaluations = 0;
    std::chrono::duration<double> total_time;
};

void MinMaxAlgorithm::searchBestMove()
{
    _maxDepth = 3; 
    int eval;

    FILE *fp;
    fp = fopen("data.txt","a");
    
    total_time = std::chrono::duration<double>(0);

    eval = minimax(0, true);

    double evaluationsPerSecond = totalEvaluations / total_time.count();

    printf("Best evaluation = %d\n", eval);
    printf("Total evaluations = %d\n", totalEvaluations);
    printf("Time taken for evaluations = %.2f seconds\n", total_time.count());
    printf("Evaluations per second = %.2f\n", evaluationsPerSecond);
    fprintf(fp, "%d %.2f\n", _maxDepth ,evaluationsPerSecond);

}

int MinMaxAlgorithm::minimax(int depth, bool isMaximizing)
{
    if (depth == _maxDepth) {
        totalEvaluations++;
        
        int eval = evaluate();
        return eval;
    }

    Move m;
    MoveList list;

    generateMoves(list);

    int bestEval;
    if (isMaximizing) {
        auto start = std::chrono::steady_clock::now();
        bestEval = minEvaluation();
        int eval;
        while (list.getNext(m)) {
            playMove(m);
            eval = -minimax(depth + 1, false);
            takeBack();
            if (eval > bestEval) {
                bestEval = eval;
                foundBestMove(depth, m, bestEval);
            }
        }
        finishedNode(depth, 0);
        auto end = std::chrono::steady_clock::now();

        total_time += end - start;
        return bestEval;
    } else {
        auto start = std::chrono::steady_clock::now();
        bestEval = maxEvaluation();
        int eval;
        while (list.getNext(m)) {
            playMove(m);
            eval = minimax(depth + 1, true);
            takeBack();
            if (eval < bestEval) {
                bestEval = eval;
                foundBestMove(depth, m, bestEval);
            }
        }
        finishedNode(depth, 0);
        auto end = std::chrono::steady_clock::now();

        total_time += end - start;
        return bestEval;
    }
}
// register ourselve as a search strategy
MinMaxAlgorithm minmaxAlgorithm;
