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


typedef struct {
    int score;
    Move move;
} scored_move;

#pragma omp declare reduction(max_smove : scored_move :\
omp_out = omp_out.score > omp_in.score ? omp_out : omp_in)\
initializer (omp_priv=(omp_orig))


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

    int times_called = 0;
    int totalEvaluations = 0;
    std::chrono::duration<double> total_time = std::chrono::duration<double>(0);
};


void MinMaxStrategy::searchBestMove()
{

    times_called++;
    _maxDepth = 3;
    int eval;
    auto start = std::chrono::steady_clock::now();
    eval = minimax(0, *_board);
    auto end = std::chrono::steady_clock::now();
    total_time += end - start;
    //printf("best evaluation = %d\n" , eval);
    printf("Evals/sec: %f\n", totalEvaluations/total_time.count());
    printf("Time: %f\n", total_time.count()/times_called);
}


int MinMaxStrategy::minimax(int depth, Board& board)
{

    if (depth == 0){

        Move m;
        MoveList list;
        
        scored_move bestSMove;
        bestSMove.score = -maxEvaluation();
        scored_move smoveL;

        generateMoves(list, board);
        int size = list.getLength();
        Move move_array[size];
        int i = 0;

        while(list.getNext(m)){
            move_array[i] = m;
            i++;
        }

        #pragma omp parallel
        {

            #pragma omp for private(i,smoveL) reduction(max_smove:bestSMove) firstprivate(board)
            for (i = 0; i <size; i++){
                smoveL.move = move_array[i];
                playMove(move_array[i], board);
                /*#pragma omp critical
                {
                    //board.print();
                }*/
                smoveL.score = -minimax(depth + 1, board);
                /*#pragma omp critical
                {        
                    printf("eval = %d\n", smoveL.score);
                }*/
                if (smoveL.score > bestSMove.score){
                    bestSMove.score = smoveL.score;
                    bestSMove.move = smoveL.move;
                }
                takeBack(board);
            }
            /*#pragma omp critical
            {
                printf("%d id: bestEval: %d\n", omp_get_thread_num(), bestSMove.score);
                std::cout << "bestmove: " << bestSMove.move.name() << std::endl;
            }*/
        }
        foundBestMove(depth, bestSMove.move, bestSMove.score);
        finishedNode(depth,0);
        return bestSMove.score;
    }

    if (depth == _maxDepth) {
        #pragma omp atomic
            totalEvaluations++;
        return -evaluate(board);
    }
    Move m;
    MoveList list;
    generateMoves(list, board);
    int bestEval;
    bestEval = -maxEvaluation();
    int eval;

    while(list.getNext(m)) {
        playMove(m, board);
        eval = -minimax(depth + 1, board);
        takeBack(board);
        if(eval > bestEval){
            bestEval = eval;
            foundBestMove(depth,m,bestEval);
        }
    }
    finishedNode(depth,0);
    return bestEval;
}

// register ourselve as a search strategy
MinMaxStrategy minMaxStrategy;