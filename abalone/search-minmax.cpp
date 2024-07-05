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
#include <omp.h>

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
    int minimax(int depth, Board& board);

    int totalEvaluations = 0;
    std::chrono::duration<double> total_time = std::chrono::duration<double>(0);
};


void MinMaxStrategy::searchBestMove()
{

    _maxDepth = 1;
    int eval;
    auto start = std::chrono::steady_clock::now();
    eval = minimax(0, *_board);
    auto end = std::chrono::steady_clock::now();
    total_time += end - start;
    printf("best evaluation = %d\n" , eval);
    printf("Evals/sec: %f\n", totalEvaluations/total_time.count());
    printf("Time: %f\n", total_time.count());
}


int MinMaxStrategy::minimax(int depth, Board& board)
{

    if (depth == 0){

        Move m;
        MoveList list;
        generateMoves(list, board);
        int bestEval=-maxEvaluation();
        int eval;

        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp taskgroup task_reduction(max:bestEval)
                {
                while(list.getNext(m)) {
                    #pragma omp task firstprivate(m,board,eval) in_reduction(max:bestEval)
                    {
                        playMove(m, board);
                        /*#pragma omp critical
                        {
                            //board.print();
                        }*/
                        eval = -minimax(depth + 1, board);
                        #pragma omp critical
                        {        
                            printf("eval = %d\n", eval);
                        }
                        if (eval > bestEval){
                            bestEval = eval;
                            //bestMoveL = m;
                            //printf("%d id: bestEval: %d\n", omp_get_thread_num(), bestEval);
                        }
                        takeBack(board);
                    }
                }
                }
            }
            #pragma omp critical
            {
                printf("%d id: bestEval: %d\n", omp_get_thread_num(), bestEval);
            }
        }
        finishedNode(depth,0);
        return 0;
    }

    if (depth == _maxDepth) {
        totalEvaluations++;
        return -evaluate(board);
    }

    Move m;
    MoveList list;
    generateMoves(list, board);

    while(list.getNext(m)) {
        playMove(m, board);
        -minimax(depth + 1, board);
        takeBack(board);
    }
    finishedNode(depth,0);
    return 0;
}

// register ourselve as a search strategy
MinMaxStrategy minMaxStrategy;